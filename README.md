## Code for analysis presented in the thesis "Identifying biomarkers of response and resistance to anti-LAG-3-based therapy in advanced melanoma." by Nabila Zulkapeli.

1. <B>filter_adata.py</b>: merge AnnData with metadata and filter out cores with low melanoma cells, low total cells, remove high TIL cores, and other QC to aid downstream analysis
2. <B>cell_type_proportions.py</b>: fit a linear mixed effects model per cell type with fixed effects (response (R/NR) and region (high tumour, peritumour)), random effects (patient ID) and optional covariates (e.g. batch, slide ID, site of biopsy)
3. <b>differential_expression.py</b>: built using <B>diffxpy</b>, run differential expression gene (DEG) analysis using Wilcoxon rank-sum test to compile top DEG lists in R/NR
4. <b>nhood_per_core</b>: using <B>Squidpy</b> functions, build spatial neighbours per core instead of across all cores and merge and return results as an AnnData object

### References:
1. diffxpy: https://github.com/theislab/diffxpy <br>
2. Fang, Z., Liu, X., & Peltz, G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. <i>Bioinformatics</i>, <b>39</b>(1), btac757, (2023). https://doi.org/10.1093/bioinformatics/btac757
3. Troulé, K., Petryszak, R., Cakir, B. <i>et al.</i> CellPhoneDB v5: inferring cell–cell communication from single-cell multiomics data. <i>Nat Protoc</i> (2025). https://doi.org/10.1038/s41596-024-01137-1
4. Palla, G., Spitzer, H., Klein, M. <i>et al.</i> Squidpy: a scalable framework for spatial omics analysis. <i>Nat Methods</i> <b>19</b>, 171–178 (2022). https://doi.org/10.1038/s41592-021-01358-2

