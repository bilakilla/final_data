import pandas as pd
import numpy as np
if not hasattr(np, "float"):
    np.float = float
if not hasattr(np, "int"):
    np.int = int
if not hasattr(np, "bool"):
    np.bool = bool
import diffxpy.api as de
import logging
import scipy.stats
import scipy.sparse as sp

def run_DEG_wilcoxon(
    adata,
    celltype_col = 'specific_cell_types',
    response_col = 'Response',
    min_cells_per_group = 10,
    top_n = None
):
    """
    run differential expression analysis using Wilcoxon rank-sum test to see which genes are upregulated in responders vs. non-responders.
    
    returns:
    all_deg : pd.DataFrame
        all DEGs across cell types
    """
    
    deg_results = {}
    
    for CT in adata.obs[celltype_col].dropna().unique():
        print('Running DEG on each cell type...')
        # subset adata by cell type (CT):
        adata_CT = adata[adata.obs[celltype_col] == CT].copy()
        # convert sparse to dense array for each cell type:
        adata_CT.X = adata_CT.X.toarray()
        
        response_groups = adata_CT.obs[response_col].unique()
        if len(response_groups) < 2:
            print(f"Skipping {CT} as it isn't present in both responders and non-responders!")
            continue
        
        wilcox_test = de.test.rank_test(
            data=adata_CT,
            grouping=response_col
        )
        deg_results[CT] = {'method': 'wilcoxon', 'result': wilcox_test}
        print(f'Wilcoxon rank-sum test successful for {CT}!')
        
    all_deg = []
    for CT, res in deg_results.items():
        if res is None:
            continue
        df = res['result'].summary()
        df['pval'] = df['pval'].replace(0, 1e-300)      # prevent infinite results
        df['pval'] = df['pval'].fillna(1)               # treat genes expressed in 0 cells as ns (p=1)
        df['log2fc'] = df['log2fc'].fillna(0)           # and no effect size (lfc=0)
        df['score'] = np.sign(df['log2fc']) * -np.log10(df['pval'])     # for GSEA
        df['cell_type'] = CT
        df['method'] = res['method']
        if top_n is not None:
            df = df.sort_values('qval').head(top_n)
        all_deg.append(df)
    
    all_deg = pd.concat(all_deg, ignore_index=True)
    
    return all_deg