from collections import defaultdict
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from pathlib import Path
from glob import glob
import scanpy as sc
import seaborn as sns

def filter_adata(
    adata,
    clinical,
    clinical_patient='Melpin',
    adata_patient='melpin',
    core='core_id',
    cluster_col='X_scANVI_predicted_2'
):
    """
    merge adata with clinical file to gain metadata like 'response'
    clean up the core names and filter out cores for downstream analysis

    args:
        adata (_type_): _description_
        clinical (_type_): _description_
        response_col (str, optional): _description_. Defaults to 'Response'.
        patient_col (str, optional): _description_. Defaults to 'Melpin'.
        
    returns:
        adata_filtered: both peritumour and high tumour cores
        adata_hightumour: only high tumour cores
        adata_peritumour: only peritumour cores
    """

    # convert patient_col to a string and standardise casing
    clinical[clinical_patient] = clinical[clinical_patient].astype(str)
    adata.obs[adata_patient] = adata.obs[adata_patient].astype(str)
    adata.obs.rename(columns={adata_patient: clinical_patient}, inplace=True)
    adata.obs.rename(columns={'Response': 'response_adata'}, inplace=True)

    # merge both on clinical_patient column
    adata.obs = adata.obs.join(
        clinical.set_index(clinical_patient),
        on=clinical_patient
    )
    adata.obs.rename(columns={'Response': 'Response_RECIST'}, inplace=True)
    adata.obs.rename(columns={'RESPONSE': 'Response'}, inplace=True)
    adata.obs['Response'] = adata.obs['Response'].astype('category')
    adata.obs['Response'] = adata.obs['Response'].cat.rename_categories(
        {"Non-responder": "Non-Responder"})
    print(f'adata and clinical merged, rename the columns if you want')

    # rename all the cores to contain region of interest
    rename_core_dict = {
        '0029039_Region_1_4-A': 'High Tumour: 0029039_Region_1_4-A',
        '0029039_Region_1_4-B': 'Peritumour: 0029039_Region_1_4-B',
        '0029039_Region_1_4-C': 'High Tumour: 0029039_Region_1_4-C',
        '0029039_Region_1_4-D': 'Peritumour: 0029039_Region_1_4-D',
        '0029039_Region_1_4-E': 'High Tumour: 0029039_Region_1_4-E',
        '0029039_Region_1_4-F': 'Peritumour: 0029039_Region_1_4-F',
        '0029039_Region_1_5-A': 'High Tumour: 0029039_Region_1_5-A',
        '0029039_Region_1_5-B': 'Peritumour: 0029039_Region_1_5-B',
        '0029039_Region_1_5-C': 'High Tumour: 0029039_Region_1_5-C',
        '0029039_Region_1_5-D': 'Peritumour: 0029039_Region_1_5-D',
        '0029039_Region_1_5-E': 'High Tumour: 0029039_Region_1_5-E',
        '0029039_Region_1_5-F': 'Peritumour: 0029039_Region_1_5-F',
        '0029039_Region_1_6-A': 'High Tumour: 0029039_Region_1_6-A',
        '0029039_Region_1_6-B': 'Peritumour: 0029039_Region_1_6-B',
        '0029039_Region_1_6-C': 'High Tumour: 0029039_Region_1_6-C',
        '0029039_Region_1_6-D': 'Peritumour: 0029039_Region_1_6-D',
        '0029039_Region_1_6-E': 'High Tumour: 0029039_Region_1_6-E',
        '0029039_Region_1_6-F': 'Peritumour: 0029039_Region_1_6-F',
        '0029039_Region_2_7-C': 'High Tumour: 0029039_Region_2_7-C',
        '0029039_Region_2_7-D': 'Peritumour: 0029039_Region_2_7-D',
        '0029039_Region_2_7-E': 'High Tumour: 0029039_Region_2_7-E',
        '0029039_Region_2_7-F': 'Peritumour: 0029039_Region_2_7-F',
        '0029039_Region_2_8-D': 'Peritumour: 0029039_Region_2_8-D',
        '0029039_Region_2_8-E': 'High Tumour: 0029039_Region_2_8-E',
        '0029039_Region_2_8-F': 'Peritumour: 0029039_Region_2_8-F',
        '0029039_Region_2_9-D': 'Peritumour: 0029039_Region_2_9-D',
        '0029039_Region_2_9-E': 'High Tumour: 0029039_Region_2_9-E',
        '0029039_Region_2_9-F': 'Peritumour: 0029039_Region_2_9-F',
        '0029039_Region_2_10-D': 'Peritumour: 0029039_Region_2_10-D',
        '0029039_Region_2_10-E': 'High Tumour: 0029039_Region_2_10-E',
        '0029039_Region_2_11-D': 'High Tumour: 0029039_Region_2_11-D',
        '0029039_Region_2_11-E': 'High Tumour: 0029039_Region_2_11-E',
        '0029039_Region_3_11-A': 'High Tumour: 0029039_Region_3_11-A',
        '0029039_Region_3_11-B': 'Peritumour: 0029039_Region_3_11-B',
        '0029039_Region_3_7-A': 'High Tumour: 0029039_Region_3_7-A',
        '0029039_Region_3_7-B': 'Peritumour: 0029039_Region_3_7-B',
        '0029039_Region_3_8-A': 'High Tumour: 0029039_Region_3_8-A',
        '0029039_Region_3_8-B': 'Peritumour: 0029039_Region_3_8-B',
        '0029039_Region_3_8-C': 'High Tumour: 0029039_Region_3_8-C',
        '0029039_Region_3_9-A': 'High Tumour: 0029039_Region_3_9-A',
        '0029039_Region_3_9-B': 'Peritumour: 0029039_Region_3_9-B',
        '0029039_Region_3_9-C': 'High Tumour: 0029039_Region_3_9-C',
        '0029039_Region_3_10-A': 'High Tumour: 0029039_Region_3_10-A',
        '0029039_Region_3_10-B': 'Peritumour: 0029039_Region_3_10-B',
        '0029039_Region_3_10-C': 'High Tumour: 0029039_Region_3_10-C',
        '0040204_Region_1_14-C': 'High TIL: 0040204_Region_1_14_C',
        '0040204_Region_1_14-D': 'High TIL: 0040204_Region_1_14-D',
        '0040207_Region_4_14-A': 'High Tumour: 0040207_Region_4_14-A',
        '0040207_Region_4_14-B': 'High Tumour: 0040207_Region_4_14-B',
        '0040207_Region_4_14-C': 'High Tumour: 0040207_Region_4_14-C',
        '0040207_Region_4_14-D': 'High Tumour: 0040207_Region_4_14-D',
        '0052306_Region_4_14-A': 'Peritumour: 0052306_Region_4_14-A',
        '0052306_Region_4_14-B': 'Peritumour: 0052306_Region_4_14-B',
        '0052306_Region_4_14-C': 'Peritumour: 0052306_Region_4_14-C',
        '0052306_Region_4_14-D': 'Peritumour: 0052306_Region_4_14-D'
    }
    adata.obs[core] = adata.obs[core].cat.rename_categories(rename_core_dict)

    # dictionaries for specific and broad cell type annotation
    specific_cell_types = {
    'mk_0': 'Proliferating Melanoma', # More dedifferentiated than proliferative: SLC7A11, CCND1, SLC1A5, KIT, CD44
    'mk_1': 'Endothelial', # CLEC14A, PLVAP, SPARCL1
    'mk_2': 'CD4 T', # Helper CD4 T cell: CD3E, CD96, IL7R, CD6, CD5
    'mk_3': 'Inflammatory CAF', # IFN-stimulated CAF enriched in ISGs
    'mk_5': 'Classical CAF', # Classical ECM-remodeling CAFs: DCN, LUM, VCAN, PDGFRA
    'mk_7': 'Granulocyte', # CEACAM8, CCL16, CCL28
    'mk_8': 'Proliferating Melanoma', # Proliferating melanoma cells: CENPF, CDK1
    'mk_9': 'Proliferating Melanoma', # Proliferating melanoma cells: S100B, SOX9, MAGEA12, MET, ERBB3, ERBB4, CDK4, CDK6, EIF4EBP1, PSAT1, LDHB, STMN1
    'mk_11': 'Proliferating Melanoma', # Proliferating melanoma cells: CENPF, SOX9
    'mk_12': 'Proliferating Melanoma', # Proliferating melanoma cells: ERBB3, CCND1, MAT2A, CTNNB1, NOTCH2
    'mk_13': 'Epithelial', # Infiltrating mucosal epithelial cells: EPCAM, MUC5AC, TFF3, REG4, CDX2
    'mk_14': 'Proliferating Melanoma', # Proliferating melanoma cells: UBE2C, CDK1, MKI67, CENPF, ORC6, CDKN2C, STMN1, TUBA1B
    'mk_15': 'Inflammatory CAF', # IFN-stimulated CAF enriched in ISGs
    'mk_17': 'Dendritic', # Monocyte-derived DCs: ITGAX, FCER1A, CD80, CD86, CLEC10A, IRF8, BATF3
    'mk_18': 'CD8 T', # Exhausted CD8 T cell
    'mk_20': 'CD8 T', # Exhausted CD8 T cell
    'mk_21': 'M2 TAM', # M2 TAMs: CHIT1, CD163, MARCO, TREM2, VSIG4, LILRB4, LILRB2, GPR34
    'mk_22': 'Epithelial', # Infiltrating mucosal epithelial cells: TFF3, MUC5AC, REG4, CDX2, EPCAM, ACE2
    'mk_24': 'Plasmablast', # Plasmablast: MS4A1, CD37, BANK1, CD79A/B, CD19, 
    'mk_25': 'Plasma B', # Plasma cell: TNFRSF17, IGHG...
    'mk_28': 'M1 TAM', # M1 TAMs: CXCL9, CXCL10, CXCL11
    'mk_29': 'Plasma B', # Plasma cell: IGHG4, IGHG1, JCHAIN, IGKC, IGHG3, LUM, IGHGP, DCN, IGLC3, IGHM, IGHG2, CCL21, CXCL12, CCL19, MZB1, CXCL14, CD79A
    'mk_30': 'Ig-expressing TAM', # IG-expressing TAMs: IGHG1, IGKC, JCHAIN, IGHG3, C1QB, CSF1R, IGHG4, C1QA, CD163, IGHGP, CD14
    'mk_31': 'M2 TAM', # M2 TAMs: CD163, CHIT1, GPR34, TREM2, VSIG4, MARCO, CX3CR1, CSF1R, CD14, CD68
    'mk_32': 'TLS', # Tertiary lymphoid structure: B cells, plasma cells, CD8⁺ cytotoxic T cells, Tfh cells, and regulatory T cells
    'mk_34': 'Proliferating Melanoma', # Proliferating melanoma cells: ORC6, ERBB3, MKI67, CCND1
    'mk_36': 'CD8 T', # IFN-activated CD8 T cells: CD8A, CD8B, PRF1, GZMB, GZMK, GZMH, GZMA, GNLY, NKG7, CTSW
    'mk_37': 'Mast' # Tumour-associated mast cells: CPA3, MS4A2, HPGDS, IL1R1
    }
    broad_cell_types = {
    'mk_0': 'Melanoma', # More dedifferentiated than proliferative: SLC7A11, CCND1, SLC1A5, KIT, CD44
    'mk_1': 'Endothelial', # CLEC14A, PLVAP, SPARCL1
    'mk_2': 'T', # Helper CD4 T cell: CD3E, CD96, IL7R, CD6, CD5
    'mk_3': 'Inflammatory CAF', # IFN-stimulated CAF enriched in ISGs
    'mk_5': 'Classical CAF', # Classical ECM-remodeling CAFs: DCN, LUM, VCAN, PDGFRA
    'mk_7': 'Neutrophil', # CEACAM8, CCL16, CCL28
    'mk_8': 'Melanoma', # Proliferating melanoma cells: CENPF, CDK1
    'mk_9': 'Melanoma', # Proliferating melanoma cells: S100B, SOX9, MAGEA12, MET, ERBB3, ERBB4, CDK4, CDK6, EIF4EBP1, PSAT1, LDHB, STMN1
    'mk_11': 'Melanoma', # Proliferating melanoma cells: CENPF, SOX9
    'mk_12': 'Melanoma', # Proliferating melanoma cells: ERBB3, CCND1, MAT2A, CTNNB1, NOTCH2
    'mk_13': 'Epithelial', # Infiltrating mucosal epithelial cells: EPCAM, MUC5AC, TFF3, REG4, CDX2
    'mk_14': 'Melanoma', # Proliferating melanoma cells: UBE2C, CDK1, MKI67, CENPF, ORC6, CDKN2C, STMN1, TUBA1B
    'mk_15': 'Inflammatory CAF', # IFN-stimulated CAF enriched in ISGs
    'mk_17': 'Myeloid', # Monocyte-derived DCs: CD1C, FCER1A, CD80, CD86, CLEC10A, IRF8, BATF3
    'mk_18': 'T', # Exhausted CD8 T cell
    'mk_20': 'T', # Exhausted CD8 T cell
    'mk_21': 'Myeloid', # M2 TAMs: CHIT1, CD163, MARCO, TREM2, VSIG4, LILRB4, LILRB2, GPR34
    'mk_22': 'Epithelial', # Infiltrating mucosal epithelial cells: TFF3, MUC5AC, REG4, CDX2, EPCAM, ACE2
    'mk_24': 'B', # Plasmablast: MS4A1, CD37, BANK1, CD79A/B, CD19, 
    'mk_25': 'B', # Plasma cell: TNFRSF17, IGHG...
    'mk_28': 'Myeloid', # M1 TAMs: CXCL9, CXCL10, CXCL11
    'mk_29': 'B', # Plasma cell: IGHG4, IGHG1, JCHAIN, IGKC, IGHG3, LUM, IGHGP, DCN, IGLC3, IGHM, IGHG2, CCL21, CXCL12, CCL19, MZB1, CXCL14, CD79A
    'mk_30': 'Myeloid', # IG-expressing TAMs: IGHG1, IGKC, JCHAIN, IGHG3, C1QB, CSF1R, IGHG4, C1QA, CD163, IGHGP, CD14
    'mk_31': 'Myeloid', # M2 TAMs: CD163, CHIT1, GPR34, TREM2, VSIG4, MARCO, CX3CR1, CSF1R, CD14, CD68
    'mk_32': 'TLS', # Tertiary lymphoid structure: B cells, plasma cells, CD8⁺ cytotoxic T cells, Tfh cells, and regulatory T cells
    'mk_34': 'Epithelial', # Infiltrating mucosal epithelial cells: MUC5AC, REG4, CDX2, EPCAM, ACE2, PGA5
    'mk_36': 'T', # IFN-activated CD8 T cells: CD8A, CD8B, PRF1, GZMB, GZMK, GZMH, GZMA, GNLY, NKG7, CTSW
    'mk_37': 'Mast' # Tumour-associated mast cells: CPA3, MS4A2, HPGDS, IL1R1
    }

    # order cluster_col
    sorted_cluster_labels = sorted(adata.obs[cluster_col].unique(), key=lambda x: int(x.split("_")[1]))
    adata.obs[cluster_col] = pd.Categorical(
        adata.obs[cluster_col],
        categories=sorted_cluster_labels,
        ordered=True
    )

    # map both to the adata
    adata.obs['specific_cell_types'] = adata.obs[cluster_col].map(specific_cell_types)
    unique_cell_types = adata.obs['specific_cell_types'].astype('category').cat.categories
    adata.uns['specific_cell_types_colors'] = sns.color_palette('tab20', n_colors=len(unique_cell_types)).as_hex()
    adata.obs['broad_cell_types'] = adata.obs[cluster_col].map(broad_cell_types)
    print(f'clusters are now annotated with specific and broad cell types')

    # filter out cores with <100 melanoma cells or <1000 total cells, as well as high TIL cores
    melanoma_clusters = adata.obs['broad_cell_types'] == 'Melanoma'

    # count melanoma cells per core
    melanoma_counts = (
        adata.obs.loc[melanoma_clusters]
        .groupby('core_id')
        .size()
        .rename('melanoma_cell_count')
    )
    adata.obs = adata.obs.join(
        melanoma_counts.rename("melanoma_cell_count"),
        on="core_id",
        how="left"
    )
    adata.obs['melanoma_cell_count'] = adata.obs['melanoma_cell_count'].fillna(0)

    high_tumour_ids = adata.obs.loc[adata.obs['core_name'] == 'High tumour', 'core_id'].unique()
    peritumour_ids  = adata.obs.loc[adata.obs['core_name'] == 'Peritumour', 'core_id'].unique()

    rank_high_tumour = melanoma_counts.loc[high_tumour_ids].sort_values(ascending=False)
    rank_peri        = melanoma_counts.loc[peritumour_ids].sort_values(ascending=False)

    total_counts_per_core = (
        adata.obs
        .groupby('core_id')
        .size()
        .rename('total_cell_count')
    )
    adata.obs['total_cell_count'] = adata.obs['core_id'].map(total_counts_per_core)

    adata.obs['melanoma_cell_count'] = pd.to_numeric(adata.obs['melanoma_cell_count'], errors='coerce')
    adata.obs['total_cell_count']    = pd.to_numeric(adata.obs['total_cell_count'], errors='coerce')

    adata.obs['melanoma_percentage'] = (
        adata.obs['melanoma_cell_count'] / adata.obs['total_cell_count'] * 100
    )

    melanoma_summary_counts = (
        adata.obs[['core_id', 'melanoma_cell_count', 'total_cell_count', 'melanoma_percentage']]
        .drop_duplicates(subset='core_id')
        .sort_values('melanoma_percentage', ascending=True)
    )

    #  remove cores with <100 melanoma cells or <1000 total cells
    melanoma_summary_counts.loc[(melanoma_summary_counts['melanoma_cell_count'] < 100) | (melanoma_summary_counts['total_cell_count'] < 1000)]

    adata_filtered = adata.copy()
    adata_filtered = adata_filtered[(adata_filtered.obs.melanoma_cell_count >= 100) & (adata_filtered.obs.total_cell_count >= 1000)]
    print(f'low quality cores with <1000 total cells or <100 melanoma cells have been removed')

    adata_filtered = adata_filtered[(adata_filtered.obs.core_name == 'Peritumour') | (adata_filtered.obs.core_name == 'High tumour')]
    print(f'adata_filtered (both high tumour and peritumour cores) has been created')

    adata_hightumour = adata_filtered[adata_filtered.obs.core_name=='High tumour']
    print(f'adata_hightumour (only high tumour cores) has been created')

    adata_peritumour = adata_filtered[adata_filtered.obs.core_name=='Peritumour']
    print(f'adata_peritumour (only peritumour cores) has been created')

    print(f'ready for downstream analysis, you can choose adata_filtered for global analysis and adata_peritumour/adata_hightumour for region-specific analyses')

    return adata_filtered, adata_hightumour, adata_peritumour

adata_path = '/Users/nabilazulkapeli/Documents/Honours Thesis 2025/nabs_data/res_2/final_scanvi_adata.h5ad'
adata = sc.read_h5ad(adata_path)
clinical = pd.read_excel("lag3_clinical.xlsx")

adata_filtered, adata_hightumour, adata_peritumour = filter_adata(
    adata=adata, 
    clinical=clinical
    )