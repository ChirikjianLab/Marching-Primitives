import numpy as np
from .sdf_superquadric import eul2rotm


def _dtheta(theta, arclength, threshold, scale, sigma):
    """Compute angular step for uniform arclength sampling."""
    if theta < threshold:
        dt = abs((arclength / scale[1] + theta ** sigma) ** (1.0 / sigma) - theta)
    else:
        ct = np.cos(theta)
        st = np.sin(theta)
        dt = arclength / sigma * np.sqrt(
            (ct ** 2 * st ** 2) /
            (scale[0] ** 2 * ct ** (2 * sigma) * st ** 4 +
             scale[1] ** 2 * st ** (2 * sigma) * ct ** 4)
        )
    return dt


def _angle2points(theta, scale, sigma):
    """Convert angles to superellipse points."""
    ct = np.cos(theta)
    st = np.sin(theta)
    points = np.zeros((2, len(theta)))
    points[0] = scale[0] * np.sign(ct) * np.abs(ct) ** sigma
    points[1] = scale[1] * np.sign(st) * np.abs(st) ** sigma
    return points


def _uniform_sampled_superellipse(epsilon, scale, arclength):
    """
    Uniformly sample a superellipse by arclength.

    Parameters
    ----------
    epsilon : float
        Shape parameter.
    scale : array-like, shape (2,)
        Scale factors [a, b].
    arclength : float
        Target arclength between samples.

    Returns
    -------
    point : ndarray (2, M)
        Sampled superellipse points.
    """
    threshold = 1e-2
    num_limit = 10000
    theta = np.zeros(num_limit)

    # Forward pass
    critical = num_limit
    for i in range(1, num_limit):
        dt = _dtheta(theta[i - 1], arclength, threshold, scale, epsilon)
        theta_temp = theta[i - 1] + dt
        if theta_temp > np.pi / 4:
            critical = i
            break
        if i < num_limit - 1:
            theta[i] = theta_temp
        else:
            raise RuntimeError(
                f"Sampled points exceed limit of {num_limit * 4}. "
                "Increase arclength or raise the limit."
            )

    # Backward pass with flipped scale
    scale_flip = scale[::-1]
    num_pt = critical
    for j in range(critical, num_limit):
        dt = _dtheta(theta[j - 1], arclength, threshold, scale_flip, epsilon)
        theta_temp = theta[j - 1] + dt
        if theta_temp > np.pi / 4:
            num_pt = j
            break
        if j < num_limit - 1:
            theta[j] = theta_temp
        else:
            raise RuntimeError(
                f"Sampled points exceed limit of {num_limit * 4}. "
                "Increase arclength or raise the limit."
            )

    theta = theta[:num_pt]

    # Forward points
    points_fw = _angle2points(theta[:critical], scale, epsilon)
    # Backward points (with flipped scale, then flip back)
    points_bw_raw = _angle2points(theta[critical - 1:], scale_flip, epsilon)
    points_bw = points_bw_raw[:, ::-1]  # flip columns
    points_bw = np.array([points_bw[1], points_bw[0]])  # swap rows

    point = np.hstack([points_fw, points_bw])

    # Build full superellipse (4 quadrants)
    n = num_pt - 1
    q2 = np.array([-point[0, :n], point[1, :n]])[:, ::-1]
    q3 = np.array([-point[0, 1:], -point[1, 1:]])
    q4 = np.array([point[0, :n], -point[1, :n]])[:, ::-1]

    point = np.hstack([point, q2, q3, q4])
    return point


def single_mesh_superquadrics(x, arclength=0.02, taper=False):
    """
    Generate a 3D mesh for a single superquadric.

    Parameters
    ----------
    x : array-like, shape (11,) or (13,)
        Superquadric parameters: [eps1, eps2, ax, ay, az, eul_z, eul_y, eul_x, tx, ty, tz]
        With taper: additional [taper_x, taper_y].
    arclength : float
        Arclength sampling resolution.
    taper : bool
        Whether tapering is enabled.

    Returns
    -------
    x_mesh, y_mesh, z_mesh : ndarray
        2D arrays of mesh coordinates.
    """
    x = np.array(x, dtype=float)

    if taper:
        if len(x) != 13:
            raise ValueError("Input parameters should have dimension 13 for tapered SQ.")
    else:
        if len(x) != 11:
            raise ValueError("Input parameters should have dimension 11 for SQ.")

    # Clamp small shape parameters for numerical stability
    if x[0] < 0.01:
        x[0] = 0.01
    if x[1] < 0.01:
        x[1] = 0.01

    R = eul2rotm(x[5:8])
    t = x[8:11]

    point_eta = _uniform_sampled_superellipse(x[0], np.array([1.0, x[4]]), arclength)
    point_omega = _uniform_sampled_superellipse(x[1], np.array([x[2], x[3]]), arclength)

    n_omega = point_omega.shape[1]
    n_eta = point_eta.shape[1]

    x_mesh = np.zeros((n_omega, n_eta))
    y_mesh = np.zeros((n_omega, n_eta))
    z_mesh = np.zeros((n_omega, n_eta))

    for m in range(n_omega):
        for n in range(n_eta):
            pt = np.array([
                point_omega[0, m] * point_eta[0, n],
                point_omega[1, m] * point_eta[0, n],
                point_eta[1, n]
            ])

            if taper:
                fx = x[11] * pt[2] / x[4] + 1.0
                fy = x[12] * pt[2] / x[4] + 1.0
                pt[0] *= fx
                pt[1] *= fy

            pt = R @ pt + t
            x_mesh[m, n] = pt[0]
            y_mesh[m, n] = pt[1]
            z_mesh[m, n] = pt[2]

    return x_mesh, y_mesh, z_mesh
