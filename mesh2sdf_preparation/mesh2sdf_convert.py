import numpy as np
import os
import argparse
import csv
import trimesh
import mesh2sdf
import sys

def main(argv):
    parser = argparse.ArgumentParser(
        description='Transforming the meshes of ShapeNet in to Signed Distance Fields.')

    parser.add_argument(
        'path_to_data',
        help='Path to the directory containing the test data.'
    )

    parser.add_argument(
        '--grid_resolution', type=int, default=64,
        help='Set the resolution of the voxel grids in the order of x, y, z, e.g. 64 means 64^3.'
    )

    parser.add_argument(   
        '--normalize', action='store_true'
    )

    parser.add_argument(
        '--level', type=float, default=2,
        help='Set watertighting thicken level. By default 2'
    )

    parser.set_defaults(normalize=False)

    args = parser.parse_args(argv)

    assert os.path.isfile(args.path_to_data), {'The input file does not exist. Please Verify the path ' + args.path_to_data + '.'}
    
    dir_name = os.path.dirname(args.path_to_data)
    mesh_file = os.path.basename(args.path_to_data)
    mesh_file_woext = os.path.splitext(mesh_file)

    if not args.normalize:
        ply = os.path.join(dir_name, mesh_file_woext[0] + '_watertight.ply')
        stl = os.path.join(dir_name, mesh_file_woext[0] + '_watertight.stl')
        csv_path = os.path.join(dir_name, mesh_file_woext[0] + '.csv')
    else:
        ply = os.path.join(dir_name, mesh_file_woext[0] + '_normalized_watertight.ply')
        stl = os.path.join(dir_name, mesh_file_woext[0] + '_normalized_watertight.stl')
        csv_path = os.path.join(dir_name, mesh_file_woext[0] + '_normalized.csv')

    mesh = trimesh.load(os.path.join(dir_name, mesh_file), force='mesh')
    print('Original mesh loaded from ' + args.path_to_data + '.')

    mesh_scale = 0.8
    size = args.grid_resolution
    level = args.level / size


    # normalize mesh
    vertices = mesh.vertices
    bbmin = vertices.min(0)
    bbmax = vertices.max(0)
    center = (bbmin + bbmax) * 0.5
    scale = 2.0 * mesh_scale / (bbmax - bbmin).max()
    vertices = (vertices - center) * scale
    
    # generate watertight mesh and sdf
    print('Converting to watertight mesh...')
    sdf, mesh = mesh2sdf.compute(
        vertices, mesh.faces, size, fix=True, level=level, return_mesh=True)

    # output mesh and sdf
    if args.normalize:
        # output mesh
        mesh.export(ply)
        mesh.export(stl)
        print('Watertight mesh generated and saved to ' + ply + '.')

        # output sdf
        resolution = args.grid_resolution
        grid_config = np.array(
            [[resolution], 
            [-1], [1], 
            [-1], [1], 
            [-1], [1]]
        )
        writevoxel = np.reshape(np.swapaxes(sdf, 0, 2), (resolution**3, 1))
        writevoxel = np.append(grid_config, writevoxel).reshape(-1, 1)
    
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(writevoxel)
        print('SDF saved to ' + csv_path + '.')
    else:
        # output mesh
        mesh.vertices = mesh.vertices / scale + center
        mesh.export(ply)
        print('Watertight mesh generated and saved to ' + ply + '.')

        # output sdf
            # output sdf
        resolution = args.grid_resolution
        grid_config = np.array(
            [[resolution], 
            [-1 / scale + center[0]], [1 / scale + center[0]], 
            [-1 / scale  + center[1]], [1 / scale + center[1]], 
            [-1 / scale  + center[2]], [1 / scale + center[2]]]
        )
        writevoxel = np.reshape(np.swapaxes(sdf / scale, 0, 2), (resolution**3, 1))
        writevoxel = np.append(grid_config, writevoxel).reshape(-1, 1)
    
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(writevoxel)
        print('SDF saved to ' + csv_path + '.')   

if __name__ == "__main__":
    main(sys.argv[1:])
    