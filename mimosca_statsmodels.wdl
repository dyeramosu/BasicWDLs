version 1.0

workflow mimosca {
    input {
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/mimosca:1.0.0"
        Int preemptible = 1
        Int disk_space = 128

        String output_dir # gbucket (no / at end)

        File perturb_gex_anndata_file # al_ld_073_processed_deepika.h5ad 
        File guide_info_file # crispr-guides_al-ld-073_final.csv
        
        Int num_chunks # number of chunks (32659 genes, 11 chunks, 2969 genes in each chunk)
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
                perturb_gex_anndata_file = perturb_gex_anndata_file,
                guide_info_file = guide_info_file
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
        File perturb_gex_anndata_file # al_ld_073_processed_deepika.h5ad
        File guide_info_file # crispr-guides_al-ld-073_final.csv

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
        guide_info = pd.read_csv('~{guide_info_file}')

        # create gex df (with added grna_bam column)
        genes = adata.var.index.to_list()
        gene_df = sc.get.obs_df(
                adata,
                keys=['grna_bam', *genes]
            )
        gene_df = gene_df.reset_index()
        gene_df = pd.merge(gene_df, guide_info[['guide', 'grna']], left_on='grna_bam', right_on='guide').drop(columns=['grna_bam', 'guide'])
        gene_df.set_index('cellname', inplace=True)

        # create cell by guide df
        cell_by_guide = pd.get_dummies(gene_df['grna'])
        cell_by_guide = cell_by_guide.fillna(0)

        # add cell state cov to design matrix
        cell_by_guide = cell_by_guide.join(adata.obs['states'])
        state_mapping = {'E': 0, 'QM': 1, 'other': 2}
        cell_by_guide['states'] = cell_by_guide['states'].replace(state_mapping)

        # get list of genes depending on chunk
        idxs = np.arange(0, 32660, 2969)
        gene_idxs = np.arange(idxs[~{chunk}], idxs[~{chunk}+1])
        gex_df = gene_df.iloc[:, gene_idxs]

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