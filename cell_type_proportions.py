import pandas as pd
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests

def run_cell_prop(
    adata,
    celltype_col,
    groupby='core_id',
    response_col='Response',
    patient_col='Melpin',
    region_col='core_name',
    covariates=None,
    transform=True,
):
    """
    fit linear mixed-effects models per cell type:
      fixed effects: response, region, response:region (+ optional covariates)
      random effects: patient ID (accounts for multiple cores per patient)

    returns:
        results_df:  wide-format fixed-effect results per cell type
        tidy_df:     long-format results (betas, SEs, p, CI, FDR, stars) for plotting
        proportions: core_id x celltype proportion matrix
        meta:        metadata (response, region) indexed by core_id
        
        echo response_col, groupby, region_col: for downstream plotting code
    """
    if covariates is None:
        covariates = []
    
    # rename specific_cell_types
    new_labels_map = {
    "Proliferating Melanoma": "Melanoma",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Classical CAF": "cCAF",
    "Inflammatory CAF": "iCAF",
    "Mast": "Mast",
    "Granulocyte":"Granulocyte",
    "Dendritic":"Dendritic",
    "M1 TAM":"M1 TAM",
    "M2 TAM":"M2 TAM",
    "Ig-expressing TAM": "Ig-TAM",
    "Plasmablast":"Plasmablast",
    "Plasma B":"Plasma",
    "TLS":"TLS",
    "CD4 T":"CD4 T",
    "CD8 T":"CD8 T"
    }
    adata.obs[celltype_col] = adata.obs[celltype_col].map(new_labels_map)
    
    # core_id x cell type counts to make proportions table
    counts = (
        adata.obs.groupby([groupby, celltype_col])
        .size()
        .unstack(fill_value=0)
    )
    proportions = counts.div(counts.sum(axis=1), axis=0)

    # metadata (response, patient ID, region) of each core
    meta = (
        adata.obs[[groupby, response_col, patient_col, region_col] + covariates]
        .drop_duplicates()
        .set_index(groupby)
    )
    # make sure response and region are categorical
    for col in [response_col, region_col]:
        if col in meta.columns and not pd.api.types.is_categorical_dtype(meta[col]):
            meta[col] = pd.Categorical(meta[col])

    # join proportions + meta tables
    base_df = proportions.join(meta, how='inner')

    results = []
    tidy_records = []

    # significance stars for downstream plotting
    def starbars(p):
        if p < 0.001: return '***'
        elif p < 0.01: return '**'
        elif p < 0.05: return '*'
        else: return ''

    for CT in proportions.columns:
        # skip cell types that are absent
        if (proportions[CT] > 0).sum() == 0:
            continue

        col = f"prop_{CT.replace(' ', '_')}"    # remove spaces in cell types for the model
        df = base_df.copy()
        df[col] = proportions[CT]

        # arcsine transform proportions to stabilise variance
        if transform:
            # clip to avoid exact 0 or 1 -> keep domain [0,1]
            p = np.clip(df[col].to_numpy(dtype=float), 0.0, 1.0)
            df[col] = np.arcsin(np.sqrt(p))

        cov_str = ' + '.join(covariates) if covariates else ''
        if df[region_col].nunique() > 1:    # global analysis (adata_filtered)
            fixed_formula = f"Q('{col}') ~ {response_col} * {region_col}" + (f' + {cov_str}' if cov_str else '')
        else:                               # region-specific analysis (adata_peritumour, adata_hightumour)
            fixed_formula = f"Q('{col}') ~ {response_col}" + (f' + {cov_str}' if cov_str else '')

        md = MixedLM.from_formula(
            fixed_formula,
            data=df,
            groups=df[patient_col],  # random intercepts for patient to account for repeated measures

            re_formula="1"
        )
        fit = md.fit(reml=False, method='lbfgs', disp=False)

        # collect fixed-effect (fe) parameters
        fe_params = fit.fe_params
        fe_se = fit.bse_fe if hasattr(fit, 'bse_fe') else fit.bse.loc[fe_params.index]
        fe_p = fit.pvalues.loc[fe_params.index]

        results_row = {'cell type': CT}
        for param in fe_params.index:
            beta = fe_params[param]
            se = fe_se[param]
            pval = fe_p[param]
            ci_low, ci_high = beta - 1.96*se, beta + 1.96*se

            mean_ref = fe_params['Intercept']   # baseline mean is the intercept of the slope

            if ':' in param:  # interaction term
                # find the main effects involved
                parts = param.split(':')
                mean_alt = mean_ref + beta
                for part in parts:
                    if part in fe_params.index:
                        mean_alt += fe_params[part] # for interaction, add to main effects
            elif param == 'Intercept':  
                mean_alt = mean_ref
            else:  # main effect
                mean_alt = mean_ref + beta  # for main effects (region or response)
                
            # back-transform to get proportions again for intuitive plots
            p_ref = np.sin(mean_ref)**2
            p_alt = np.sin(mean_alt)**2
            beta_prop = p_alt - p_ref
            p_low = np.sin(mean_ref + ci_low)**2 - p_ref
            p_high = np.sin(mean_ref + ci_high)**2 - p_ref

            results_row[param] = beta
            results_row[param + '_se'] = se
            results_row[param + '_pval'] = pval

            tidy_records.append({
                'cell_type': CT,
                'term': param,
                'beta': beta,
                'se': se,
                'pval': pval,
                'ci_low': ci_low,
                'ci_high': ci_high,
                'beta_prop': beta_prop,     # use for forest plots, box plots
                'ci_low_prop': p_low,       # use for forest plots
                'ci_high_prop': p_high      # use for forest plots
            })

        results.append(results_row)

    results_df = pd.DataFrame(results)
    tidy_df = pd.DataFrame(tidy_records)

    # FDR correction for any term that involves the response (main effect or interactions)
    if not tidy_df.empty:
        response_terms = [t for t in tidy_df['term'].unique() if response_col in t]
        tidy_df['fdr'] = np.nan
        for term in response_terms:
            mask = tidy_df['term'] == term
            if mask.any():
                tidy_df.loc[mask, 'fdr'] = multipletests(
                    tidy_df.loc[mask, 'pval'].astype(float), method='fdr_bh'
                )[1]
        # stars only for terms where pval is defined
    tidy_df['sig'] = tidy_df['pval'].apply(lambda p: starbars(p) if pd.notnull(p) else '')

    return results_df, tidy_df, proportions, meta, response_col, groupby, region_col