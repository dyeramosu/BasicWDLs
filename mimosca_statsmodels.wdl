version 1.0

workflow mimosca {
    call run_statsmodels

    output {
        File mimosca_outputs = run_statsmodels.coeffs
    }
}

task run_statsmodels {
    input {
        String output_dir # gbucket (no / at end)
        
        File perturb_gex_anndata_file # al_ld_073_processed_deepika.h5ad
        File cell_by_guide_csv_file # cell_by_guide_df.csv

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 1
        Int disk_space = 128
    }

    command <<<
        set -e 
        
        mkdir tmpdir
        mkdir mimosca_output_wdl

        python << CODE
        # imports
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import statsmodels.api as sm

        print('loading in anndata', flush=True)

        # load files
        adata = sc.read_h5ad('~{perturb_gex_anndata_file}') # Y = adata.X
        cell_by_guide = pd.read_csv('~{cell_by_guide_csv_file}', index_col=0) # X
        gex_df = pd.DataFrame.sparse.from_spmatrix(adata.X, columns=adata.var.index, index=adata.obs.index) # Y

        print('loaded in anndata, starting regression', flush=True)
        
        # fit regression model
        coeff_dict = {}
        pval_dict = {}
        
        x = cell_by_guide
        x = sm.add_constant(x)
        
        for gene in gex_df.columns:
            print(gene, flush=True)
            y = gex_df[gene]

            model = sm.OLS(y, x).fit() 

            coeff_dict[gene] = model.params
            pval_dict[gene] = model.pvalues

        print('finished regression, saving coefficients', flush=True)

        coeff_df = pd.DataFrame(coeff_dict)
        pval_df = pd.DataFrame(pval_dict)

        coeff_df.to_pickle("mimosca_output_wdl/statsmodels_coeffs.pkl")
        pval_df.to_pickle("mimosca_output_wdl/statsmodels_pvals.pkl")

        print('saved coefficients', flush=True)
       
        CODE

        gsutil -m rsync mimosca_output_wdl ~{output_dir}
        tar -zcvf mimosca_outputs.tar.gz mimosca_output_wdl
    >>>

    output {
        File coeffs = 'mimosca_outputs.tar.gz'
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