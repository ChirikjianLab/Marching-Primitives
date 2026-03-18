import numpy as np
import struct


def plyread(path, mode=None):
    """
    Read a PLY 3D data file.

    Parameters
    ----------
    path : str
        Path to PLY file.
    mode : str, optional
        If 'tri', return (triangles, points) for mesh display.

    Returns
    -------
    If mode == 'tri':
        tri : ndarray (M, 3) - triangular face indices (0-based)
        pts : ndarray (N, 3) - vertex coordinates
    Else:
        elements : dict - element data keyed by element name
    """
    with open(path, 'rb') as f:
        # Read header
        line = f.readline().decode('ascii').strip()
        if line != 'ply':
            raise ValueError("Not a PLY file.")

        fmt = None
        elements_info = []  # list of (name, count, properties)
        current_element = None
        comments = []

        while True:
            line = f.readline().decode('ascii').strip()
            if line.startswith('format'):
                parts = line.split()
                fmt = parts[1]
            elif line.startswith('comment'):
                comments.append(line[8:])
            elif line.startswith('element'):
                parts = line.split()
                current_element = {
                    'name': parts[1],
                    'count': int(parts[2]),
                    'properties': []
                }
                elements_info.append(current_element)
            elif line.startswith('property'):
                parts = line.split()
                if parts[1] == 'list':
                    current_element['properties'].append({
                        'type': 'list',
                        'count_type': parts[2],
                        'data_type': parts[3],
                        'name': parts[4]
                    })
                else:
                    current_element['properties'].append({
                        'type': parts[1],
                        'name': parts[2]
                    })
            elif line.startswith('end_header'):
                break

        # Type mappings
        ply_to_numpy = {
            'char': 'i1', 'uchar': 'u1', 'short': 'i2', 'ushort': 'u2',
            'int': 'i4', 'uint': 'u4', 'float': 'f4', 'double': 'f8',
            'char8': 'i1', 'uchar8': 'u1', 'short16': 'i2', 'ushort16': 'u2',
            'int32': 'i4', 'uint32': 'u4', 'float32': 'f4', 'double64': 'f8',
            'int8': 'i1', 'uint8': 'u1', 'int16': 'i2', 'uint16': 'u2',
            'float64': 'f8',
        }

        ply_to_struct = {
            'char': 'b', 'uchar': 'B', 'short': 'h', 'ushort': 'H',
            'int': 'i', 'uint': 'I', 'float': 'f', 'double': 'd',
            'char8': 'b', 'uchar8': 'B', 'short16': 'h', 'ushort16': 'H',
            'int32': 'i', 'uint32': 'I', 'float32': 'f', 'double64': 'd',
            'int8': 'b', 'uint8': 'B', 'int16': 'h', 'uint16': 'H',
            'float64': 'd',
        }

        ply_type_size = {
            'char': 1, 'uchar': 1, 'short': 2, 'ushort': 2,
            'int': 4, 'uint': 4, 'float': 4, 'double': 8,
            'char8': 1, 'uchar8': 1, 'short16': 2, 'ushort16': 2,
            'int32': 4, 'uint32': 4, 'float32': 4, 'double64': 8,
            'int8': 1, 'uint8': 1, 'int16': 2, 'uint16': 2,
            'float64': 8,
        }

        elements = {}

        if fmt == 'ascii':
            for elem_info in elements_info:
                elem_data = {p['name']: [] for p in elem_info['properties']}
                for _ in range(elem_info['count']):
                    line = f.readline().decode('ascii').strip().split()
                    idx = 0
                    for prop in elem_info['properties']:
                        if prop['type'] == 'list':
                            count = int(line[idx])
                            idx += 1
                            vals = [float(line[idx + k]) for k in range(count)]
                            idx += count
                            elem_data[prop['name']].append(vals)
                        else:
                            elem_data[prop['name']].append(float(line[idx]))
                            idx += 1

                # Convert to numpy arrays
                for prop in elem_info['properties']:
                    if prop['type'] != 'list':
                        elem_data[prop['name']] = np.array(elem_data[prop['name']])

                elements[elem_info['name']] = elem_data

        else:
            # Binary format
            endian = '<' if fmt == 'binary_little_endian' else '>'

            for elem_info in elements_info:
                has_list = any(p['type'] == 'list' for p in elem_info['properties'])

                if not has_list:
                    # Fast path: all fixed-size properties
                    dtype_list = []
                    for prop in elem_info['properties']:
                        np_type = ply_to_numpy.get(prop['type'])
                        if np_type is None:
                            raise ValueError(f"Unknown type: {prop['type']}")
                        bo = '<' if fmt == 'binary_little_endian' else '>'
                        dtype_list.append((prop['name'], bo + np_type))

                    dt = np.dtype(dtype_list)
                    raw = np.frombuffer(f.read(dt.itemsize * elem_info['count']),
                                        dtype=dt, count=elem_info['count'])
                    elem_data = {}
                    for prop in elem_info['properties']:
                        elem_data[prop['name']] = raw[prop['name']].astype(float)
                    elements[elem_info['name']] = elem_data
                else:
                    # Slow path: has list properties
                    elem_data = {p['name']: [] for p in elem_info['properties']}
                    for _ in range(elem_info['count']):
                        for prop in elem_info['properties']:
                            if prop['type'] == 'list':
                                ct = prop['count_type']
                                dt_name = prop['data_type']
                                sz = ply_type_size[ct]
                                fmt_char = endian + ply_to_struct[ct]
                                count = struct.unpack(fmt_char, f.read(sz))[0]
                                sz2 = ply_type_size[dt_name]
                                fmt_char2 = endian + str(count) + ply_to_struct[dt_name]
                                vals = list(struct.unpack(fmt_char2, f.read(sz2 * count)))
                                elem_data[prop['name']].append(vals)
                            else:
                                sz = ply_type_size[prop['type']]
                                fmt_char = endian + ply_to_struct[prop['type']]
                                val = struct.unpack(fmt_char, f.read(sz))[0]
                                elem_data[prop['name']].append(val)

                    for prop in elem_info['properties']:
                        if prop['type'] != 'list':
                            elem_data[prop['name']] = np.array(elem_data[prop['name']])
                    elements[elem_info['name']] = elem_data

    # Handle 'tri' mode
    if mode and mode.lower() == 'tri':
        # Find vertex element
        pts = None
        for name in ['vertex', 'Vertex', 'point', 'Point', 'pts', 'Pts']:
            if name in elements:
                ve = elements[name]
                if 'x' in ve and 'y' in ve and 'z' in ve:
                    pts = np.column_stack([ve['x'], ve['y'], ve['z']])
                break
        if pts is None:
            pts = np.zeros((1, 3))

        # Find face element
        tri = np.zeros((0, 3), dtype=int)
        for name in ['face', 'Face', 'poly', 'Poly', 'tri', 'Tri']:
            if name in elements:
                fe = elements[name]
                for pname in ['vertex_indices', 'vertex_indexes', 'vertex_index',
                              'indices', 'indexes']:
                    if pname in fe:
                        face_lists = fe[pname]
                        tris = []
                        for face in face_lists:
                            face = [int(v) for v in face]
                            # Triangulate polygon
                            for j in range(1, len(face) - 1):
                                tris.append([face[0], face[j], face[j + 1]])
                        tri = np.array(tris, dtype=int)
                        break
                break

        return tri, pts

    return elements
