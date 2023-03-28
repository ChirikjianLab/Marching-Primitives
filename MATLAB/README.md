## MATLAB Implementation
### Marching-Primitives algorithm
The source code of the algorithm is in `/MATLAB/src/MPS.m`
```
x = MPS(sdf, grid)
```
The algorithm depends on the Image Processing Toolbox of MATLAB.
The algorithm requires a Signed Distance Function discretized on a voxel grid as input. More specifically, for a grid of size $(x,y,z):M*N*W$, `grid.size = [M, N, W]` is the size of the voxel grid; `grid.range = [x_min, x_max, y_min, y_max, z_min, z_max]` stores the range of the voxel grid, and `sdf` is a 1-D array flattened from the 3-D array storing the signed distance of points in the voxel grid. 


The output of the function is a 2D array of size $K*11$, where each row stores the parameter of a superquadric $[\epsilon_1, \epsilon_2, a_x, a_y, a_z, euler_z, euler_y, euler_x, t_x, t_y, t_z]$.

### Preparing SDF from meshes
If you do not have SDF files but want to test the algorithm, we have a simple script [mesh2sdf_preparation](/ mesh2sdf_preparation) to generate SDF from meshes. The script is based on this [package](https://github.com/wang-ps/mesh2sdf).
The mesh file will be first transformed to be watertight so that a valid SDF can be extracted.
To run the script, first run
```
pip install mesh2sdf
```
Then simply run
```
python3 mesh2sdf_convertion.py $location of the mesh file$ --normalize --grid_resolution 100
```
where the mesh and sdf will be normalized within $[-1, 1]$ if you add `–-normalized`; and `--grid_resolution` specifies the resolution of the sdf (default $64$).
The script accepts mesh forms: .stl, .off, .ply, .collada, .json, .dict, .glb, .dict64, .msgpack, .obj.
The SDF (*.csv) and preprocessed watertight mesh (*_watertight.stl) will be saved at the same folder of the input mesh.
A few .obj meshes from ShapeNets are prepared in the [data](/MATLAB/data) at `/MATALB/data`.

### Demo script
After obtaining the SDF file, run the [demo script](/MATLAB/demo_script.m) at `/MATLAB/demo_script.m`.
The script conducts the shape abstraction and visualize the results.
The recovered superquadric representation is saved as `*.mat` at the same location of the input SDF.
For visualization, mesh file of the superquadric representation is also saved as `*_sq.stl`.

## Citation
If you find this repo useful, please cite:
```
@Inproceedings{Liu2023CVPR,
     title = {Marching-Primitives: Shape Abstraction from Signed Distance Function},
     author = {Liu, Weixiao and Wu, Yuwei and Ruan, Sipu and Chirikjian, Gregory},
     booktitle = {Proceedings IEEE Conf. on Computer Vision and Pattern Recognition (CVPR)},
     year = {2023}
}
```
