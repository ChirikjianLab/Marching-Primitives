import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .sdf_superquadric import eul2rotm
from .single_mesh_superquadrics import single_mesh_superquadrics


def show_tsdf(sdf, grid, ax=None):
    """
    Visualize truncated signed distance field as 3D scatter plot.

    Parameters
    ----------
    sdf : ndarray
        Flattened SDF values.
    grid : dict
        Grid metadata with 'x', 'y', 'z', 'size', 'disp_range' keys.
    ax : matplotlib 3D axis, optional
        Axis to plot on.
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    sdf_display = sdf.copy()
    mask = (sdf_display > grid['disp_range'][0]) & (sdf_display < grid['disp_range'][1])

    # Build 3D coordinates
    gx, gy, gz = np.meshgrid(grid['x'], grid['y'], grid['z'], indexing='ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    pts = coords[mask]
    vals = sdf_display[mask]

    sc = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=vals, alpha=0.2, s=1)
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax)
    return ax


def show_superquadrics(x, ax=None, taper=False, color='r', view_axis=(0, 0),
                       cam_roll=0, show_axis=False, arclength=0.02,
                       face_alpha=1.0, face_lighting='flat', lighting=False):
    """
    Visualize a single superquadric as a 3D surface plot.

    Parameters
    ----------
    x : array-like, shape (11,) or (13,)
        Superquadric parameters.
    ax : matplotlib 3D axis, optional
    taper : bool
    color : str or tuple
    view_axis : tuple (elev, azim)
    cam_roll : float
    show_axis : bool
    arclength : float
    face_alpha : float
    face_lighting : str
    lighting : bool
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    x_mesh, y_mesh, z_mesh = single_mesh_superquadrics(
        x, arclength=arclength, taper=taper
    )

    ax.plot_surface(x_mesh, y_mesh, z_mesh, color=color, alpha=face_alpha,
                    shade=(face_lighting != 'none'))

    ax.set_aspect('equal')
    ax.view_init(elev=view_axis[1], azim=view_axis[0])

    if not show_axis:
        ax.set_axis_off()

    return ax


def show_multi_superquadrics(x, color='r', arclength=0.02, ax=None):
    """
    Visualize multiple superquadrics in a single figure.

    Parameters
    ----------
    x : ndarray (K, 11)
        Superquadric parameters, one per row.
    color : str or tuple
    arclength : float
    ax : matplotlib 3D axis, optional
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    for i in range(x.shape[0]):
        show_superquadrics(x[i], ax=ax, color=color,
                           face_alpha=1.0, show_axis=True,
                           lighting=False, arclength=arclength)
    return ax
