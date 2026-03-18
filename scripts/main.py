"""
Concise Marching-Primitives pipeline.

Usage:
    python main.py <path_to_sdf.csv> [--ply <path.ply>] [--no-save] [--no-display]
"""

import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt

from marching_primitives import (
    MPS, mesh_superquadrics, single_mesh_superquadrics, mesh2tri, plyread,
    load_sdf_csv, reduce_mesh, save_stl,
)


def main():
    parser = argparse.ArgumentParser(description='Marching-Primitives')
    parser.add_argument('csv_file', help='Path to the SDF .csv file')
    parser.add_argument('--ply', help='Path to watertight .ply for comparison', default=None)
    parser.add_argument('--no-save', action='store_true', help='Do not save output files')
    parser.add_argument('--no-display', action='store_true', help='Do not display plots')
    args = parser.parse_args()

    # Load SDF
    sdf, voxel_grid = load_sdf_csv(args.csv_file)
    visualize_arclength = 0.01 * np.sqrt(voxel_grid['range'][1] - voxel_grid['range'][0])

    # Run MPS
    t0 = time.time()
    x = MPS(sdf, voxel_grid)
    print(f"MPS: {time.time() - t0:.2f}s, {x.shape[0]} superquadrics")

    # Generate and compress mesh
    mesh_orig = mesh_superquadrics(x, arclength=visualize_arclength)
    mesh_faces, mesh_vertices = reduce_mesh(mesh_orig['f'], mesh_orig['v'], ratio=0.1)

    # Save
    if not args.no_save:
        pathname = os.path.dirname(args.csv_file)
        name = os.path.splitext(os.path.basename(args.csv_file))[0]

        np.savez(os.path.join(pathname, f"{name}_sq_py.npz"), x=x.astype(np.float32))
        save_stl(os.path.join(pathname, f"{name}_sq_py.stl"), mesh_faces, mesh_vertices)
        np.savetxt(
            os.path.join(pathname, f"{name}_sq_py.csv"), x, delimiter=',',
            header='eps1,eps2,ax,ay,az,eul_z,eul_y,eul_x,tx,ty,tz', comments='',
        )

    if args.no_display:
        return

    # Visualize
    grid_range = voxel_grid['range']

    def setup_axes(ax, title):
        ax.set_xlim(grid_range[0], grid_range[1])
        ax.set_ylim(grid_range[2], grid_range[3])
        ax.set_zlim(grid_range[4], grid_range[5])
        ax.set_aspect('equal')
        ax.view_init(elev=-40, azim=151)
        ax.set_axis_off()
        ax.set_title(title)

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(8, 6))
    cmap = plt.get_cmap('tab20', x.shape[0])
    for i in range(x.shape[0]):
        xm, ym, zm = single_mesh_superquadrics(x[i], arclength=visualize_arclength, taper=False)
        f_i, v_i = mesh2tri(xm, ym, zm, 'f')
        s = f_i.shape[0]
        f_i = f_i[np.concatenate([np.arange(s//8, 3*s//8), np.arange(5*s//8, 7*s//8)])]
        ax.plot_trisurf(v_i[:, 0], v_i[:, 1], v_i[:, 2],
                        triangles=f_i, color=cmap(i), edgecolor='none', alpha=1.0)
    setup_axes(ax, 'Superquadrics from Marching Primitives')

    if args.ply and os.path.exists(args.ply):
        tri_gt, pts_gt = plyread(args.ply, 'tri')
        color = np.array([145, 163, 176]) / 255.0

        fig1, ax1 = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(8, 6))
        ax1.plot_trisurf(pts_gt[:, 0], pts_gt[:, 1], pts_gt[:, 2],
                         triangles=tri_gt, color=color, edgecolor='none', alpha=1.0)
        setup_axes(ax1, 'Ground truth')

        fig2, ax2 = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(8, 6))
        ax2.plot_trisurf(pts_gt[:, 0], pts_gt[:, 1], pts_gt[:, 2],
                         triangles=tri_gt, color='g', edgecolor='none', alpha=0.5)
        ax2.plot_trisurf(mesh_vertices[:, 0], mesh_vertices[:, 1], mesh_vertices[:, 2],
                         triangles=mesh_faces, color=color, edgecolor='none', alpha=1.0)
        setup_axes(ax2, 'Overlap')

    plt.show()


if __name__ == '__main__':
    main()
