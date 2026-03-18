"""
Compare MATLAB and Python superquadric fitting results against ground-truth SDF.

Usage:
    python scripts/compare_matlab_python.py data/              # all objects in data/
    python scripts/compare_matlab_python.py data/ --names chair2      # single object
    python scripts/compare_matlab_python.py data/ --names chair2 chair5 --plot
"""

import argparse
import glob
import os
import numpy as np
import scipy.io

from marching_primitives.sdf_superquadric import sdf_superquadric, sdf_multi_superquadrics


def load_ground_truth_sdf(csv_path):
    """Load ground-truth SDF from CSV. Returns grid points (3, N) and SDF values (N,)."""
    data = np.loadtxt(csv_path)
    res = int(data[0])
    bounds = data[1:7]  # xmin, xmax, ymin, ymax, zmin, zmax
    sdf_vals = data[7:]

    x = np.linspace(bounds[0], bounds[1], res)
    y = np.linspace(bounds[2], bounds[3], res)
    z = np.linspace(bounds[4], bounds[5], res)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

    return points, sdf_vals, res


def load_matlab_params(mat_path):
    """Load superquadric params from .mat file. Returns (K, 11) array."""
    d = scipy.io.loadmat(mat_path)
    keys = [k for k in d if not k.startswith('_')]
    return d[keys[0]]


def load_python_params(csv_path):
    """Load superquadric params from Python CSV. Returns (K, 11) array."""
    return np.loadtxt(csv_path, delimiter=',', skiprows=1)


def compute_sdf_metrics(sdf_pred, sdf_gt):
    """Compute comparison metrics between predicted and ground-truth SDF."""
    diff = sdf_pred - sdf_gt
    abs_diff = np.abs(diff)

    # Only evaluate near the surface (where GT SDF is small)
    near_surface = np.abs(sdf_gt) < 0.1
    n_near = near_surface.sum()

    metrics = {
        'mae': np.mean(abs_diff),
        'rmse': np.sqrt(np.mean(diff ** 2)),
        'max_err': np.max(abs_diff),
        'median_err': np.median(abs_diff),
    }
    if n_near > 0:
        metrics['mae_near_surface'] = np.mean(abs_diff[near_surface])
        metrics['rmse_near_surface'] = np.sqrt(np.mean(diff[near_surface] ** 2))
    else:
        metrics['mae_near_surface'] = float('nan')
        metrics['rmse_near_surface'] = float('nan')

    # IoU approximation: fraction of grid points where both agree on sign
    sign_agree = np.mean(np.sign(sdf_pred) == np.sign(sdf_gt))
    metrics['sign_agreement'] = sign_agree

    # Volume IoU: treat negative SDF as interior
    interior_pred = sdf_pred <= 0
    interior_gt = sdf_gt <= 0
    intersection = np.sum(interior_pred & interior_gt)
    union = np.sum(interior_pred | interior_gt)
    metrics['volume_iou'] = intersection / union if union > 0 else 1.0

    return metrics


def compare_params(matlab_params, python_params):
    """Compare parameter arrays directly (number of primitives, per-param stats)."""
    info = {
        'n_matlab': matlab_params.shape[0],
        'n_python': python_params.shape[0],
    }

    # If same number of primitives, compute per-primitive differences
    # But primitives may be in different order, so we skip direct row comparison
    # Instead, we compare aggregate statistics
    col_names = ['eps1', 'eps2', 'ax', 'ay', 'az', 'eul_z', 'eul_y', 'eul_x', 'tx', 'ty', 'tz']
    info['matlab_means'] = {col_names[i]: np.mean(matlab_params[:, i]) for i in range(11)}
    info['python_means'] = {col_names[i]: np.mean(python_params[:, i]) for i in range(11)}

    return info


def print_metrics(name, metrics, indent=2):
    prefix = ' ' * indent
    print(f"{prefix}{name}:")
    for k, v in metrics.items():
        if isinstance(v, float):
            print(f"{prefix}  {k:25s} = {v:.6f}")
        else:
            print(f"{prefix}  {k:25s} = {v}")


def compare_one_object(obj_name, data_dir, do_plot=False):
    """Run full comparison for one object."""
    base = os.path.join(data_dir, obj_name, f'{obj_name}_normalized')
    gt_csv = base + '.csv'
    py_csv = base + '_sq_py.csv'
    mat_file = base + '_sq.mat'

    missing = []
    for f, label in [(gt_csv, 'ground truth'), (py_csv, 'python result'), (mat_file, 'matlab result')]:
        if not os.path.exists(f):
            missing.append(f'{label} ({f})')
    if missing:
        print(f"\n=== {obj_name} === SKIPPED (missing: {', '.join(missing)})")
        return None

    print(f"\n{'=' * 60}")
    print(f"  {obj_name}")
    print(f"{'=' * 60}")

    # Load data
    print("  Loading ground truth SDF...")
    points, sdf_gt, res = load_ground_truth_sdf(gt_csv)

    matlab_params = load_matlab_params(mat_file)
    python_params = load_python_params(py_csv)

    # Parameter comparison
    pinfo = compare_params(matlab_params, python_params)
    print(f"  Primitives: MATLAB={pinfo['n_matlab']}, Python={pinfo['n_python']}")

    # Compute SDF from both (truncate to avoid overflow far from surface)
    truncation = 1.0
    print("  Computing MATLAB SDF on grid...")
    sdf_matlab = sdf_multi_superquadrics(matlab_params, points, truncation)
    print("  Computing Python SDF on grid...")
    sdf_python = sdf_multi_superquadrics(python_params, points, truncation)
    # Truncate GT to same range for fair comparison
    sdf_gt = np.clip(sdf_gt, -truncation, truncation)

    # Metrics: each vs GT
    m_matlab_gt = compute_sdf_metrics(sdf_matlab, sdf_gt)
    m_python_gt = compute_sdf_metrics(sdf_python, sdf_gt)
    # Metrics: MATLAB vs Python
    m_mat_py = compute_sdf_metrics(sdf_matlab, sdf_python)

    print_metrics("MATLAB vs Ground Truth", m_matlab_gt)
    print_metrics("Python vs Ground Truth", m_python_gt)
    print_metrics("MATLAB vs Python", m_mat_py)

    # Print parameter means side-by-side
    col_names = ['eps1', 'eps2', 'ax', 'ay', 'az', 'eul_z', 'eul_y', 'eul_x', 'tx', 'ty', 'tz']
    print(f"\n  Parameter means (MATLAB | Python):")
    for c in col_names:
        mv = pinfo['matlab_means'][c]
        pv = pinfo['python_means'][c]
        diff = abs(mv - pv)
        print(f"    {c:8s}: {mv:10.5f} | {pv:10.5f}  (diff={diff:.5f})")

    if do_plot:
        plot_comparison(obj_name, sdf_gt, sdf_matlab, sdf_python, res, points, data_dir)

    return {
        'name': obj_name,
        'n_matlab': pinfo['n_matlab'],
        'n_python': pinfo['n_python'],
        'matlab_vs_gt': m_matlab_gt,
        'python_vs_gt': m_python_gt,
        'matlab_vs_python': m_mat_py,
    }


def plot_comparison(obj_name, sdf_gt, sdf_matlab, sdf_python, res, points, data_dir):
    """Generate comparison plots."""
    import matplotlib.pyplot as plt

    # Take a 2D slice at z=0 (middle of grid)
    mid = res // 2
    gt_3d = sdf_gt.reshape(res, res, res)
    mat_3d = sdf_matlab.reshape(res, res, res)
    py_3d = sdf_python.reshape(res, res, res)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'{obj_name} - SDF comparison (z=0 slice)', fontsize=14)

    vmin, vmax = -0.3, 0.3

    # Row 1: SDF slices
    axes[0, 0].imshow(gt_3d[:, :, mid].T, origin='lower', cmap='RdBu', vmin=vmin, vmax=vmax)
    axes[0, 0].set_title('Ground Truth')
    axes[0, 1].imshow(mat_3d[:, :, mid].T, origin='lower', cmap='RdBu', vmin=vmin, vmax=vmax)
    axes[0, 1].set_title('MATLAB')
    axes[0, 2].imshow(py_3d[:, :, mid].T, origin='lower', cmap='RdBu', vmin=vmin, vmax=vmax)
    axes[0, 2].set_title('Python')

    # Row 2: Error maps
    diff_mat = np.abs(mat_3d[:, :, mid] - gt_3d[:, :, mid])
    diff_py = np.abs(py_3d[:, :, mid] - gt_3d[:, :, mid])
    diff_mp = np.abs(mat_3d[:, :, mid] - py_3d[:, :, mid])
    err_max = max(diff_mat.max(), diff_py.max(), diff_mp.max(), 0.01)

    axes[1, 0].imshow(diff_mat.T, origin='lower', cmap='hot', vmin=0, vmax=err_max)
    axes[1, 0].set_title('|MATLAB - GT|')
    axes[1, 1].imshow(diff_py.T, origin='lower', cmap='hot', vmin=0, vmax=err_max)
    axes[1, 1].set_title('|Python - GT|')
    im = axes[1, 2].imshow(diff_mp.T, origin='lower', cmap='hot', vmin=0, vmax=err_max)
    axes[1, 2].set_title('|MATLAB - Python|')

    fig.colorbar(im, ax=axes[1, :], shrink=0.6, label='Absolute Error')
    plt.tight_layout()

    out_path = os.path.join(data_dir, obj_name, f'{obj_name}_comparison.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"  Plot saved to {out_path}")
    plt.close()


def print_summary_table(results):
    """Print a summary table across all objects."""
    results = [r for r in results if r is not None]
    if not results:
        return

    print(f"\n{'=' * 80}")
    print("  SUMMARY TABLE")
    print(f"{'=' * 80}")
    header = (f"  {'Name':12s} | {'#Mat':>4s} {'#Py':>4s} | "
              f"{'MAT RMSE':>9s} {'PY RMSE':>9s} | "
              f"{'MAT IoU':>8s} {'PY IoU':>8s} | "
              f"{'M-P RMSE':>9s}")
    print(header)
    print(f"  {'-' * (len(header) - 2)}")

    for r in results:
        mg = r['matlab_vs_gt']
        pg = r['python_vs_gt']
        mp = r['matlab_vs_python']
        print(f"  {r['name']:12s} | {r['n_matlab']:4d} {r['n_python']:4d} | "
              f"{mg['rmse']:9.6f} {pg['rmse']:9.6f} | "
              f"{mg['volume_iou']:8.4f} {pg['volume_iou']:8.4f} | "
              f"{mp['rmse']:9.6f}")

    # Average
    avg_mat_rmse = np.mean([r['matlab_vs_gt']['rmse'] for r in results])
    avg_py_rmse = np.mean([r['python_vs_gt']['rmse'] for r in results])
    avg_mat_iou = np.mean([r['matlab_vs_gt']['volume_iou'] for r in results])
    avg_py_iou = np.mean([r['python_vs_gt']['volume_iou'] for r in results])
    avg_mp_rmse = np.mean([r['matlab_vs_python']['rmse'] for r in results])
    print(f"  {'-' * (len(header) - 2)}")
    print(f"  {'AVERAGE':12s} | {'':4s} {'':4s} | "
          f"{avg_mat_rmse:9.6f} {avg_py_rmse:9.6f} | "
          f"{avg_mat_iou:8.4f} {avg_py_iou:8.4f} | "
          f"{avg_mp_rmse:9.6f}")


def main():
    parser = argparse.ArgumentParser(description='Compare MATLAB vs Python superquadric results')
    parser.add_argument('data_dir', help='Path to data directory')
    parser.add_argument('--names', nargs='+', default=None,
                        help='Object names to compare (e.g., chair2 table1). Default: all.')
    parser.add_argument('--plot', action='store_true', help='Generate comparison plots')
    args = parser.parse_args()

    if args.names is None:
        # Auto-detect all objects that have both _sq.mat and _sq.csv
        pattern = os.path.join(args.data_dir, '*', '*_normalized_sq.mat')
        mat_files = sorted(glob.glob(pattern))
        names = [os.path.basename(os.path.dirname(f)) for f in mat_files]
    else:
        names = args.names

    if not names:
        print("No objects found to compare.")
        return

    print(f"Comparing {len(names)} object(s): {', '.join(names)}")

    results = []
    for name in names:
        r = compare_one_object(name, args.data_dir, do_plot=args.plot)
        results.append(r)

    print_summary_table(results)


if __name__ == '__main__':
    main()
