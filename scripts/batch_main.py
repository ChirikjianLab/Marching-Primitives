"""
Batch-run Marching-Primitives on all models in a data directory.

Usage:
    python scripts/batch_main.py data/
    python scripts/batch_main.py data/ --names chair1 chair2
    python scripts/batch_main.py data/ --no-save
"""

import os
import glob
import time
import argparse
import numpy as np

from marching_primitives import (
    MPS, mesh_superquadrics, load_sdf_csv, reduce_mesh, save_stl,
)


def process_one(csv_path, save=True):
    """Run MPS on a single SDF CSV and optionally save results."""
    name = os.path.splitext(os.path.basename(csv_path))[0]
    pathname = os.path.dirname(csv_path)

    sdf, voxel_grid = load_sdf_csv(csv_path)
    visualize_arclength = 0.01 * np.sqrt(voxel_grid['range'][1] - voxel_grid['range'][0])

    t0 = time.time()
    x = MPS(sdf, voxel_grid)
    elapsed = time.time() - t0

    n_sq = x.shape[0]
    print(f"  {name}: {n_sq} superquadrics in {elapsed:.2f}s")

    if save:
        mesh_orig = mesh_superquadrics(x, arclength=visualize_arclength)
        mesh_faces, mesh_vertices = reduce_mesh(mesh_orig['f'], mesh_orig['v'], ratio=0.1)

        np.savez(os.path.join(pathname, f"{name}_sq_py.npz"), x=x.astype(np.float32))
        save_stl(os.path.join(pathname, f"{name}_sq_py.stl"), mesh_faces, mesh_vertices)
        np.savetxt(
            os.path.join(pathname, f"{name}_sq_py.csv"), x, delimiter=',',
            header='eps1,eps2,ax,ay,az,eul_z,eul_y,eul_x,tx,ty,tz', comments='',
        )

    return {'name': name, 'n_sq': n_sq, 'time': elapsed}


def main():
    parser = argparse.ArgumentParser(description='Batch Marching-Primitives')
    parser.add_argument('data_dir', help='Path to data directory')
    parser.add_argument('--names', nargs='+', default=None,
                        help='Object names to process (e.g., chair1 chair2). Default: all.')
    parser.add_argument('--no-save', action='store_true', help='Do not save output files')
    args = parser.parse_args()

    if args.names is not None:
        csv_files = []
        for name in args.names:
            pattern = os.path.join(args.data_dir, name, f'{name}_normalized.csv')
            matches = glob.glob(pattern)
            if matches:
                csv_files.extend(matches)
            else:
                print(f"Warning: no CSV found for '{name}' at {pattern}")
    else:
        pattern = os.path.join(args.data_dir, '*', '*_normalized.csv')
        csv_files = sorted(glob.glob(pattern))

    if not csv_files:
        print("No SDF CSV files found.")
        return

    print(f"Processing {len(csv_files)} model(s)...")
    t_total = time.time()

    results = []
    for csv_path in csv_files:
        results.append(process_one(csv_path, save=not args.no_save))

    elapsed_total = time.time() - t_total

    # Summary
    print(f"\n{'=' * 50}")
    print(f"  {'Name':20s} | {'#SQ':>5s} | {'Time':>8s}")
    print(f"  {'-' * 44}")
    for r in results:
        print(f"  {r['name']:20s} | {r['n_sq']:5d} | {r['time']:7.2f}s")
    print(f"  {'-' * 44}")
    print(f"  {'TOTAL':20s} | {sum(r['n_sq'] for r in results):5d} | {elapsed_total:7.2f}s")


if __name__ == '__main__':
    main()
