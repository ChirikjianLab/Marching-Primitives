import numpy as np


def mesh2tri(X, Y, Z, tri_type='f'):
    """
    Convert a regular mesh defined by X, Y, Z into a triangulation.

    Parameters
    ----------
    X, Y, Z : ndarray
        2D arrays defining the mesh (from meshgrid/ndgrid).
    tri_type : str
        'f' for forward slash, 'b' for back slash, 'x' for cross division.

    Returns
    -------
    F : ndarray (M, 3)
        Face indices (0-based).
    V : ndarray (N, 3)
        Vertex coordinates.
    """
    nrows, ncols = X.shape
    # Create index grids (0-based)
    I, J = np.meshgrid(np.arange(nrows - 1), np.arange(ncols - 1), indexing='ij')
    I = I.ravel()
    J = J.ravel()

    def sub2ind(i, j):
        return i * ncols + j

    if tri_type == 'f':
        # Forward slash division
        tri1 = np.column_stack([
            sub2ind(I, J),
            sub2ind(I + 1, J + 1),
            sub2ind(I + 1, J)
        ])
        tri2 = np.column_stack([
            sub2ind(I, J),
            sub2ind(I, J + 1),
            sub2ind(I + 1, J + 1)
        ])
        F = np.vstack([tri1, tri2])
        V = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    elif tri_type == 'b':
        # Back slash division
        tri1 = np.column_stack([
            sub2ind(I, J + 1),
            sub2ind(I + 1, J),
            sub2ind(I, J)
        ])
        tri2 = np.column_stack([
            sub2ind(I + 1, J + 1),
            sub2ind(I + 1, J),
            sub2ind(I, J + 1)
        ])
        F = np.vstack([tri1, tri2])
        V = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    elif tri_type == 'x':
        # Cross division with extra center points
        n_quads = len(I)
        n_orig = X.size
        center_idx = np.arange(n_orig, n_orig + n_quads)

        tri1 = np.column_stack([sub2ind(I + 1, J), sub2ind(I, J), center_idx])
        tri2 = np.column_stack([sub2ind(I + 1, J + 1), sub2ind(I + 1, J), center_idx])
        tri3 = np.column_stack([sub2ind(I, J + 1), sub2ind(I + 1, J + 1), center_idx])
        tri4 = np.column_stack([sub2ind(I, J), sub2ind(I, J + 1), center_idx])
        F = np.vstack([tri1, tri2, tri3, tri4])

        # Compute center points
        corners = np.column_stack([
            sub2ind(I, J),
            sub2ind(I + 1, J),
            sub2ind(I + 1, J + 1),
            sub2ind(I, J + 1)
        ])
        Xf = X.ravel()
        Yf = Y.ravel()
        Zf = Z.ravel()
        Xe = np.mean(Xf[corners], axis=1)
        Ye = np.mean(Yf[corners], axis=1)
        Ze = np.mean(Zf[corners], axis=1)

        V = np.column_stack([
            np.concatenate([Xf, Xe]),
            np.concatenate([Yf, Ye]),
            np.concatenate([Zf, Ze])
        ])
    else:
        raise ValueError(f"Unknown tri_type '{tri_type}'. Use 'f', 'b', or 'x'.")

    return F, V
