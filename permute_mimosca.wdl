version 1.0

workflow mimosca {
    input {
        Int range_iter

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 0
        Int disk_space = 128
    }
    
    scatter (i in range(range_iter)) {
        call run_mimosca {input: iter=i }
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

        # load files
        adata = sc.read_h5ad('~{perturb_gex_anndata_file}')
        gex_df = pd.DataFrame.sparse.from_spmatrix(adata.X, columns=adata.var.index, index=adata.obs.index) # Y
        cell_by_guide = pd.read_csv('~{cell_by_guide_csv_file}', index_col=0) # X

        # shuffle

        # fit regression model
        lm = sklearn.linear_model.Ridge()
        lm.fit(cell_by_guide.values, gex_df.values)
        B = pd.DataFrame(lm.coef_) # 32659 rows (num_genes)
        
        # save coefficients 
        B.to_csv('mimosca_output_wdl/mimosca_coeffs_~{iter}.csv')
       
        CODE

        gsutil -m cp mimosca_output_wdl/mimosca_coeffs_~{iter}.csv ~{output_dir}
    >>>

    output {
        File mimosca_coeffs = 'mimosca_output_wdl/mimosca_coeffs_~{iter}.csv'
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