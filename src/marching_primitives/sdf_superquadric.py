import numpy as np


def eul2rotm(eul):
    """
    Convert ZYX Euler angles to rotation matrix.

    Parameters
    ----------
    eul : array-like, shape (3,)
        Euler angles [z, y, x] in radians.

    Returns
    -------
    R : ndarray (3, 3)
        Rotation matrix.
    """
    ct = np.cos(eul)
    st = np.sin(eul)
    R = np.array([
        [ct[1] * ct[0],
         st[2] * st[1] * ct[0] - ct[2] * st[0],
         ct[2] * st[1] * ct[0] + st[2] * st[0]],
        [ct[1] * st[0],
         st[2] * st[1] * st[0] + ct[2] * ct[0],
         ct[2] * st[1] * st[0] - st[2] * ct[0]],
        [-st[1],
         st[2] * ct[1],
         ct[2] * ct[1]]
    ])
    return R


def rotm2eul(R):
    """
    Convert rotation matrix to ZYX Euler angles.

    Parameters
    ----------
    R : ndarray (3, 3)
        Rotation matrix.

    Returns
    -------
    eul : ndarray (3,)
        Euler angles [z, y, x] in radians.
    """
    sy = -R[2, 0]
    cy = np.sqrt(R[0, 0] ** 2 + R[1, 0] ** 2)

    if cy > 1e-6:
        z = np.arctan2(R[1, 0], R[0, 0])
        y = np.arctan2(-R[2, 0], cy)
        x = np.arctan2(R[2, 1], R[2, 2])
    else:
        z = np.arctan2(-R[0, 1], R[1, 1])
        y = np.arctan2(-R[2, 0], cy)
        x = 0.0

    return np.array([z, y, x])


def sdf_superquadric(para, points, truncation=0):
    """
    Compute signed distance from points to a single superquadric.

    Parameters
    ----------
    para : array-like, shape (11,)
        [eps1, eps2, ax, ay, az, eul_z, eul_y, eul_x, tx, ty, tz]
    points : ndarray (3, N)
        3D point coordinates.
    truncation : float
        Truncation value. 0 means no truncation.

    Returns
    -------
    sdf : ndarray (N,)
        Signed distance values.
    """
    R = eul2rotm(para[5:8])
    t = para[8:11]
    X = R.T @ points - (R.T @ t[:, np.newaxis])

    r0 = np.linalg.norm(X, axis=0)

    eps1, eps2 = max(para[0], 1e-10), max(para[1], 1e-10)
    ax, ay, az = max(para[2], 1e-10), max(para[3], 1e-10), max(para[4], 1e-10)

    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        term_xy = np.power((np.abs(X[0]) / ax) ** 2, 1.0 / eps2) + \
                  np.power((np.abs(X[1]) / ay) ** 2, 1.0 / eps2)
        term_z = np.power((np.abs(X[2]) / az) ** 2, 1.0 / eps1)
        scale = np.power(
            np.power(term_xy, eps2 / eps1) + term_z,
            -eps1 / 2.0
        )

    sdf = r0 * (1.0 - np.nan_to_num(scale, nan=0.0, posinf=1e30, neginf=-1e30))

    if truncation != 0:
        sdf = np.clip(sdf, -truncation, truncation)

    return sdf


def sdf_multi_superquadrics(x, points, truncation=0):
    """
    Compute signed distance to the union of multiple superquadrics.

    Parameters
    ----------
    x : ndarray (K, 11)
        Superquadric parameters, one per row.
    points : ndarray (3, N)
        3D point coordinates.
    truncation : float
        Truncation value. 0 means no truncation.

    Returns
    -------
    sdf_u : ndarray (N,)
        Minimum SDF across all superquadrics.
    """
    sdf_u = sdf_superquadric(x[0], points, truncation)
    for i in range(1, x.shape[0]):
        sdf_u = np.minimum(sdf_u, sdf_superquadric(x[i], points, truncation))
    return sdf_u
