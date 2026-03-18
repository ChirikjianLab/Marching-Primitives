import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Batch process mesh files in a data folder with mesh2sdf_convert.py."
    )
    parser.add_argument(
        "--data_dir",
        help="Directory containing mesh files (default: data).",
    )
    parser.add_argument(
        "--grid_resolution",
        type=int,
        default=100,
        help="Voxel grid resolution (default: 100).",
    )
    parser.add_argument(
        "--level",
        type=float,
        default=2.0,
        help="Watertight thicken level (default: 2.0).",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Enable normalize mode for mesh2sdf_convert.py.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    converter_script = script_dir / "mesh2sdf_convert.py"
    data_dir = Path(args.data_dir).resolve()

    mesh_extensions = ("*.obj", "*.stl", "*.ply", "*.off", "*.glb", "*.gltf", "*.dae")
    mesh_files = sorted(
        f for ext in mesh_extensions for f in data_dir.glob(ext)
    )
    if not mesh_files:
        print(f"No mesh files found in {data_dir}")
        return 0

    failed = []
    for obj_file in mesh_files:
        cmd = [
            sys.executable,
            str(converter_script),
            str(obj_file),
            "--grid_resolution",
            str(args.grid_resolution),
            "--level",
            str(args.level),
        ]
        if args.normalize:
            cmd.append("--normalize")

        print(f"Processing: {obj_file.name}")
        print(" ".join(cmd))
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            failed.append(obj_file)

    if failed:
        print("\nFailed files:")
        for f in failed:
            print(f"- {f}")
        return 1

    print(f"\nDone. Processed {len(mesh_files)} mesh files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
