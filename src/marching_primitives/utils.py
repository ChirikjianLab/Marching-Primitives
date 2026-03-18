"""
Utility functions for Marching-Primitives.

Provides SDF loading, mesh decimation, and STL export.
"""

from struct import pack

import numpy as np


def load_sdf_csv(csv_file):
    """
    Load an SDF from a CSV file and build the voxel grid.

    The CSV format encodes: grid_size, range (6 values), then flattened SDF values.

    Parameters
    ----------
    csv_file : str
        Path to the SDF .csv file.

    Returns
    -------
    sdf : ndarray (N,)
        Truncated SDF values.
    voxel_grid : dict
        Voxel grid with keys: 'size', 'range', 'x', 'y', 'z', 'points',
        'interval', 'truncation', 'disp_range'.
    """
    raw = np.loadtxt(csv_file, delimiter=',').flatten()

    grid_size = int(raw[0])
    voxel_grid = {
        'size': np.array([grid_size, grid_size, grid_size]),
        'range': raw[1:7],
    }
    sdf = raw[7:]

    voxel_grid['x'] = np.linspace(voxel_grid['range'][0], voxel_grid['range'][1], grid_size)
    voxel_grid['y'] = np.linspace(voxel_grid['range'][2], voxel_grid['range'][3], grid_size)
    voxel_grid['z'] = np.linspace(voxel_grid['range'][4], voxel_grid['range'][5], grid_size)

    gx, gy, gz = np.meshgrid(voxel_grid['x'], voxel_grid['y'], voxel_grid['z'], indexing='ij')
    voxel_grid['points'] = np.column_stack([
        gx.ravel(order='F'), gy.ravel(order='F'), gz.ravel(order='F')
    ]).T  # (3, N)

    voxel_grid['interval'] = (
        (voxel_grid['range'][1] - voxel_grid['range'][0]) / (grid_size - 1)
    )
    voxel_grid['truncation'] = 1.2 * voxel_grid['interval']
    voxel_grid['disp_range'] = [-np.inf, voxel_grid['truncation']]

    sdf = np.clip(sdf, -voxel_grid['truncation'], voxel_grid['truncation'])

    return sdf, voxel_grid


def reduce_mesh(faces, vertices, ratio=0.5):
    """
    Mesh decimation. Tries PyMeshLab, then Open3D, then random face removal.

    Parameters
    ----------
    faces : ndarray (M, 3)
    vertices : ndarray (N, 3)
    ratio : float
        Target ratio of faces to keep (0 to 1).

    Returns
    -------
    faces_reduced : ndarray
    vertices : ndarray
    """
    try:
        import pymeshlab
        ms = pymeshlab.MeshSet()
        m = pymeshlab.Mesh(vertices, faces)
        ms.add_mesh(m)
        target = int(faces.shape[0] * ratio)
        ms.meshing_decimation_quadric_edge_collapse(targetfacenum=target)
        reduced = ms.current_mesh()
        return reduced.face_matrix(), reduced.vertex_matrix()
    except ImportError:
        pass

    try:
        import open3d as o3d
        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(vertices)
        mesh.triangles = o3d.utility.Vector3iVector(faces)
        target = int(faces.shape[0] * ratio)
        mesh_simplified = mesh.simplify_quadric_decimation(target)
        return (np.asarray(mesh_simplified.triangles),
                np.asarray(mesh_simplified.vertices))
    except ImportError:
        pass

    n_keep = max(1, int(faces.shape[0] * ratio))
    idx = np.random.choice(faces.shape[0], n_keep, replace=False)
    return faces[idx], vertices


def save_stl(filename, faces, vertices):
    """
    Save mesh as binary STL.

    Parameters
    ----------
    filename : str
        Output file path.
    faces : ndarray (M, 3)
        Triangle face indices.
    vertices : ndarray (N, 3)
        Vertex positions.
    """
    normals = np.zeros((faces.shape[0], 3))
    for i, face in enumerate(faces):
        v0, v1, v2 = vertices[face]
        n = np.cross(v1 - v0, v2 - v0)
        norm = np.linalg.norm(n)
        if norm > 0:
            n /= norm
        normals[i] = n

    with open(filename, 'wb') as f:
        f.write(b'\0' * 80)  # header
        f.write(pack('<I', faces.shape[0]))
        for i in range(faces.shape[0]):
            f.write(pack('<3f', *normals[i]))
            for j in range(3):
                f.write(pack('<3f', *vertices[faces[i, j]]))
            f.write(pack('<H', 0))
