version 1.0

workflow mimosca {
    input {
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 1
        Int disk_space = 128

        String output_dir # gbucket (no / at end)

        File perturb_gex_anndata_file # 20231210_day6_processed_filtered_with-barcodes.h5ad 
        
        Int num_chunks = 10 # number of chunks
    }
    
    scatter (c in range(num_chunks)) { 
        call run_mimosca {
            input:
                chunk = c,
                cpu = cpu,
                memory = memory,
                docker = docker,
                preemptible = preemptible, 
                disk_space = disk_space,
                output_dir = output_dir,
                perturb_gex_anndata_file = perturb_gex_anndata_file
        }
    }
  
  output {
    Array[File] final_output = run_mimosca.mimosca_coeffs
  }
}

task run_mimosca {
    input {
        String output_dir # gbucket (no / at end)
        
        Int chunk 
        File perturb_gex_anndata_file 

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
        import statsmodels.api as sm

        print('loading in data', flush=True)

        # load files
        adata = sc.read_h5ad('~{perturb_gex_anndata_file}')

        # create gex df (with added grna column)
        genes = adata.var.index.to_list()
        gene_df = sc.get.obs_df(
                adata,
                keys=['guide_ID', *genes]
            )
        gene_df['guide_ID'] = gene_df.guide_ID.astype(str)

        # create gex df without grna column
        gex_df = gene_df.drop(columns=['guide_ID'])

        # create cell by guide df
        cell_by_guide = pd.get_dummies(gene_df['guide_ID'])
        cell_by_guide = cell_by_guide.fillna(0)
        
        # add cell state cov to design matrix
        cell_by_guide = cell_by_guide.join(adata.obs['states'])
        state_mapping = {'E': 0, 'M': 1, 'other': 2}
        cell_by_guide['states'] = cell_by_guide['states'].replace(state_mapping)
        cell_by_guide['states'] = cell_by_guide.states.astype(int)

        # get list of genes depending on chunk
        idxs = np.arange(0, 32109, 2919) 
        idxs = np.append(idxs, 32109)
        gene_idxs = np.arange(idxs[~{chunk}], idxs[~{chunk}+1])
        gex_df = gex_df.iloc[:, gene_idxs]

        print('loaded in data, starting regression', flush=True)

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

         # save coefficients 
        coeff_df.to_pickle("mimosca_output_wdl/statsmodels_coeffs_~{chunk}.pkl")
        pval_df.to_pickle("mimosca_output_wdl/statsmodels_pvals_~{chunk}.pkl")        
       
        CODE
        gsutil -m rsync mimosca_output_wdl ~{output_dir}
        tar -zcvf mimosca_output_wdl_~{chunk}.tar.gz mimosca_output_wdl
    >>>

    output {
        File mimosca_coeffs = 'mimosca_output_wdl_~{chunk}.tar.gz'
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