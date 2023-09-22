version 1.0

workflow combine_coeffs {
    call combine
  
    output {
        File final_output = combine.combined_df
    }
}

task combine {
    input {
        String output_dir # gbucket (no / at end)
        String permuted_coeffs_dir # gbucket (no / at end)
        File B_labeled_file
        
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 5
        Int disk_space = 128 
    }

    command <<<
        set -e 
        
        mkdir combine_wdl_output_dir
        mkdir coeffs_dir
        gsutil -m cp ~{permuted_coeffs_dir}* coeffs_dir/

        python << CODE
        # imports
        import numpy as np
        import pandas as pd
        import os

        # load files
        B = pd.read_csv('~{B_labeled_file}', index_col=0)

        guide_list = B.columns
        gene_list = B.index 

        # initialize dict
        guide_dict = {guide: {gene: [] for gene in gene_list} for guide in guide_list}
        
        for file in os.listdir('coeffs_dir'): # for every permuted coeff df
            coeffs_df = pd.read_csv(os.path.join('coeffs_dir', file), index_col=0)
            coeffs_df.columns = guide_list
            coeffs_df.index = gene_list
        
            # separate/allocate guide columns to guide dfs
            for guide in guide_list:
                for gene in gene_list:
                    coeff = coeffs_df[guide][gene]
                    guide_dict[guide][gene].append(coeff)

        
        pd.DataFrame.from_dict(guide_dict).to_csv('combine_wdl_output_dir/all_permuted_coeffs_df.csv')
       
        CODE

        gsutil -m cp combine_wdl_output_dir/all_permuted_coeffs_df.csv ~{output_dir}
    >>>

    output {
        File combined_df = 'combine_wdl_output_dir/all_permuted_coeffs_df.csv'
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