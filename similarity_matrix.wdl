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

        and_matrix = np.bitwise_and(variant_array[:, np.newaxis, :], variant_array)  # Broadcast AND operation
        or_matrix = np.bitwise_or(variant_array[:, np.newaxis, :], variant_array)    # Broadcast OR operation

        sum_and = np.sum(and_matrix, axis=2)
        sum_or = np.sum(or_matrix, axis=2)

        with np.errstate(divide='ignore', invalid='ignore'):  # Handle division by zero
            sim_matrix = np.true_divide(sum_and, sum_or)

        sim_matrix[np.isnan(sim_matrix)] = 0  # Handle NaN resulting from division by zero

        # save similarity matrix  
        sim_matrix.to_csv('output_wdl/sim_matrix.csv')

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