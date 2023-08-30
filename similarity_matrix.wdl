version 1.0

workflow similarity_matrix {
    call create_sim_matrix

    output {
        File sim_matrix_final = create_sim_matrix.sim_matrix
    }
}

task create_sim_matrix {
    input {
        String output_dir # gbucket (no / at end)

        File variant_df_file # cell barcode by mutation dataframe

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 2
        Int disk_space = 128
    }

    command <<<
        set -e
        
        mkdir output_wdl

        python << CODE
        # imports
        import numpy as np
        import pandas as pd

        # load input files
        variant_df = pd.read_csv('~{variant_df_file}', index_col = 0)

        # create similarity matrix
        num_cells, num_features = variant_df.shape
        variant_array = variant_df.values
        sim_matrix = np.zeros((num_cells, num_cells))
        
        chunk_size = 50  # adjust this value based on available memory

        for start in range(0, num_cells, chunk_size):
            end = min(start + chunk_size, num_cells)
            chunk = variant_array[start:end]

            and_chunk = np.bitwise_and(chunk[:, np.newaxis, :], variant_array)
            or_chunk = np.bitwise_or(chunk[:, np.newaxis, :], variant_array)

            sum_and_chunk = np.sum(and_chunk, axis=2)
            sum_or_chunk = np.sum(or_chunk, axis=2)

            sim_matrix_chunk = sum_and_chunk / sum_or_chunk
            sim_matrix[start:end] = sim_matrix_chunk

        # save similarity matrix  
        pd.DataFrame(sim_matrix).to_csv('output_wdl/sim_matrix.csv')

        CODE

        gsutil -m cp output_wdl/sim_matrix.csv ~{output_dir}
    >>>

    output {
        File sim_matrix = 'output_wdl/sim_matrix.csv'
    }

    runtime {
        docker: docker
        memory: memory + "G"
        bootDiskSizeGb: 100
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}