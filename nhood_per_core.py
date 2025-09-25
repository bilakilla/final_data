import pandas as pd
import numpy as np
import squidpy as sq
import scipy.sparse as sp

def make_spatial_neighbors_per_core(
    adata, 
    group_key="core_id", 
    coord_type="generic", 
    delaunay=False, 
    radius=None
):
    """
    build spatial neighbors per core and combine results into the full AnnData object.
    
    parameters
    ----------
    adata : AnnData
        AnnData with obsm['spatial']
    group_key : str
        column in `adata.obs` defining groups (e.g. core_id)
    coord_type : str
        coordinate type for squidpy (default "generic")
    delaunay : t/f
        if true, build a spatial graph by connecting points based on a Delaunay triangulation of spatial coordinates
    radius : float, optional
        spatial radius cutoff (µm or pixels, depending on coords)
    """

    n = adata.n_obs
    connectivities = sp.lil_matrix((n, n))
    distances = sp.lil_matrix((n, n))

    for group in adata.obs[group_key].unique():
        idx = np.where(adata.obs[group_key] == group)[0]
        sub = adata[idx, :].copy()

        # run squidpy neighbors (choose n_neigh OR radius)
        sq.gr.spatial_neighbors(
            sub, 
            coord_type=coord_type,
            delaunay=delaunay, 
            radius=radius
        )

        # store into global matrices
        connectivities[np.ix_(idx, idx)] = sub.obsp["spatial_connectivities"]
        distances[np.ix_(idx, idx)] = sub.obsp["spatial_distances"]

    # save into AnnData
    adata.obsp["spatial_connectivities"] = connectivities.tocsr()
    adata.obsp["spatial_distances"] = distances.tocsr()
    
    print(f"spatial neighbors built per {group_key} (n={adata.obs[group_key].nunique()} groups)")
    return adata
