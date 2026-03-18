import numpy as np
from .single_mesh_superquadrics import single_mesh_superquadrics
from .mesh2tri import mesh2tri


def mesh_superquadrics(x, arclength=0.02, taper=False):
    """
    Convert multiple superquadric parameters to a merged triangulated mesh.

    Parameters
    ----------
    x : ndarray (K, 11) or (K, 13)
        Superquadric parameters, one per row.
    arclength : float
        Arclength sampling resolution.
    taper : bool
        Whether tapering is enabled.

    Returns
    -------
    mesh : dict
        Dictionary with 'v' (vertices, Nx3) and 'f' (faces, Mx3) keys.
    """
    n = x.shape[0]
    meshes = []

    v_offset = 0
    all_v = []
    all_f = []

    for i in range(n):
        x_mesh, y_mesh, z_mesh = single_mesh_superquadrics(
            x[i], arclength=arclength, taper=taper
        )
        f, v = mesh2tri(x_mesh, y_mesh, z_mesh, 'f')

        # Remove mesh sections (matching MATLAB: idx = s/8 : 3s/8, 5s/8 : 7s/8)
        s = f.shape[0]
        idx1 = np.arange(s // 8, 3 * s // 8)
        idx2 = np.arange(5 * s // 8, 7 * s // 8)
        idx = np.concatenate([idx1, idx2])
        f = f[idx]

        # Offset face indices
        f = f + v_offset
        v_offset += v.shape[0]

        all_v.append(v)
        all_f.append(f)

    mesh = {
        'v': np.vstack(all_v),
        'f': np.vstack(all_f)
    }
    return mesh
