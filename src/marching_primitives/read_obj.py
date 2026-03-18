import numpy as np


def read_obj(fname):
    """
    Parse a Wavefront OBJ file.

    Parameters
    ----------
    fname : str
        Path to the .obj file.

    Returns
    -------
    obj : dict
        Dictionary with keys:
        - 'v': vertices (N, 3)
        - 'vt': texture coordinates (M, 2 or 3)
        - 'vn': normal coordinates (L, 3)
        - 'f': dict with 'v', 'vt', 'vn' face index arrays
    """
    v = []
    vt = []
    vn = []
    f_v = []
    f_vt = []
    f_vn = []

    with open(fname, 'r') as fid:
        for line in fid:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            token = parts[0]

            if token == 'v':
                v.append([float(x) for x in parts[1:4]])
            elif token == 'vt':
                vt.append([float(x) for x in parts[1:]])
            elif token == 'vn':
                vn.append([float(x) for x in parts[1:4]])
            elif token == 'f':
                face_v = []
                face_vt = []
                face_vn = []
                for vert in parts[1:]:
                    indices = vert.split('/')
                    face_v.append(int(indices[0]))
                    if len(indices) > 1 and indices[1]:
                        face_vt.append(int(indices[1]))
                    if len(indices) > 2 and indices[2]:
                        face_vn.append(int(indices[2]))
                f_v.append(face_v)
                if face_vt:
                    f_vt.append(face_vt)
                if face_vn:
                    f_vn.append(face_vn)

    obj = {
        'v': np.array(v) if v else np.array([]).reshape(0, 3),
        'vt': np.array(vt) if vt else np.array([]),
        'vn': np.array(vn) if vn else np.array([]).reshape(0, 3),
        'f': {
            'v': np.array(f_v, dtype=int) if f_v else np.array([], dtype=int),
            'vt': np.array(f_vt, dtype=int) if f_vt else np.array([], dtype=int),
            'vn': np.array(f_vn, dtype=int) if f_vn else np.array([], dtype=int),
        }
    }
    return obj
