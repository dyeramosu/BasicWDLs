version 1.0

workflow mimosca {
    input {
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 2
        Int disk_space = 128

        String output_dir # gbucket (no / at end)

        File perturb_gex_anndata_file # al_ld_073_processed_deepika.h5ad 
        File cell_by_guide_csv_file # cell_by_guide_df.csv
        
        Int num_iter # number of permutations
    }
    
    scatter (i in range(num_iter)) { 
        call run_mimosca {
            input:
                iter = i,
                cpu = cpu,
                memory = memory,
                docker = docker,
                preemptible = preemptible, 
                disk_space = disk_space,
                output_dir = output_dir,
                perturb_gex_anndata_file = perturb_gex_anndata_file,
                cell_by_guide_csv_file = cell_by_guide_csv_file
        }
    }
  
  output {
    Array[File] final_output = run_mimosca.mimosca_coeffs
  }
}

task run_mimosca {
    input {
        String output_dir # gbucket (no / at end)
        
        Int iter 
        File perturb_gex_anndata_file # al_ld_073_processed_deepika.h5ad
        File cell_by_guide_csv_file # cell_by_guide_df.csv

        Int cpu 
        Int memory 
        String docker 
        Int preemptible 
        Int disk_space 
    }

    command <<<
        set -e 
        
        mkdir tmpdir
        mkdir mimosca_output_wdl

        python << CODE
        # imports
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import sklearn
        from sklearn import linear_model
        from random import shuffle

        print('loading in anndata', flush=True)

        # load files
        adata = sc.read_h5ad('~{perturb_gex_anndata_file}') # Y = adata.X
        cell_by_guide = pd.read_csv('~{cell_by_guide_csv_file}', index_col=0) # X

        print('loaded in anndata, shuffling', flush=True)

        # shuffle rows (cell names) in X 
        shuffled_cell_names = list(cell_by_guide.index)
        shuffle(shuffled_cell_names)
        cell_by_guide.index = shuffled_cell_names
        
        # reorder gex_df to be same as shuffled cell_by_guide_df
        adata = adata[shuffled_cell_names, :].copy()

        print('starting regression', flush=True)

        # fit regression model
        lm = linear_model.Ridge(fit_intercept=True, max_iter=10000)
        lm.fit(cell_by_guide.values, adata.X.toarray())
        B = pd.DataFrame(lm.coef_) # 32659 rows (num_genes)

        print('finished regression, saving pickle object', flush=True)
        
        # save coefficients 
        #B.to_csv('mimosca_output_wdl/mimosca_coeffs_~{iter}.csv')
        #cell_by_guide.to_csv('mimosca_output_wdl/cell_by_guide_~{iter}.csv')
        B.to_pickle("mimosca_output_wdl/mimosca_coeffs_~{iter}.pkl")
       
        CODE

        gsutil -m cp mimosca_output_wdl/mimosca_coeffs_~{iter}.pkl ~{output_dir}
        #gsutil -m cp mimosca_output_wdl/cell_by_guide_~{iter}.csv ~{output_dir}
    >>>

    output {
        File mimosca_coeffs = 'mimosca_output_wdl/mimosca_coeffs_~{iter}.pkl'
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