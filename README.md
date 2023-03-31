# Marching-Primitives: Shape Abstraction from Signed Distance Function
Project|Paper|Supplementary|[Arxiv](https://arxiv.org/abs/2303.13190)|[3D-Demos](/examples)|[Data](/MATLAB/data)

<img src="/examples/example.jpg" alt="example" width="600"/>

This repo provides the source code for the CVPR2023 paper:
> [**Marching-Primitives: Shape Abstraction from Signed Distance Function**](https://arxiv.org/abs/2303.13190 "ArXiv version of the paper.")  
> Weixiao Liu<sup>1,2</sup>, Yuwei Wu<sup>2</sup>, [Sipu Ruan](https://ruansp.github.io/)<sup>2</sup>, [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)<sup>2</sup>  
> <sup>1</sup> National University of Singapore, <sup>2</sup> Johns Hopkins University

## Update
 - **March 27th, 2023** - V0.1 basic version is online, including MATLAB implementation of the algorithm, data(SDF) preprocess script, and visualization tools.
Python implementation is planned in the coming updates.
 - **March 28th, 2023** - implementation details has been updated.

## Abstraction
Representing complex objects with basic geometric primitives has long been a topic in computer vision. Primitive-based representations have the merits of compactness and computational efficiency in higher-level tasks such as physics simulation, collision checking, and robotic manipulation. Unlike previous works which extract polygonal meshes from a signed distance function (SDF), in this paper, we present a novel method, named Marching-Primitives, to obtain a primitive-based abstraction directly from an SDF. Our method grows geometric primitives (such as superquadrics) iteratively by analyzing the connectivity of voxels while marching at different levels of signed distance. For each valid connected volume of interest, we march on the scope of voxels from which a primitive is able to be extracted in a probabilistic sense and simultaneously solve for the parameters of the primitive to capture the underlying local geometry. We evaluate the performance of our method on both synthetic and real-world datasets. The results show that the proposed method outperforms the state-of-the-art in terms of accuracy, and is directly generalizable among different categories and scales.

## Implementation
### Marching-Primitives algorithm
The source code of the algorithm is in `/MATLAB/src/MPS.m`
```
x = MPS(sdf, grid)
```
The algorithm depends on the Image Processing Toolbox of MATLAB.
The algorithm requires a Signed Distance Function discretized on a voxel grid as input. More specifically, for a grid of size $(x,y,z):M\times N\times W$, `grid.size = [M, N, W]` is the size of the voxel grid; `grid.range = [x_min, x_max, y_min, y_max, z_min, z_max]` stores the range of the voxel grid, and `sdf` is a 1-D array flattened from the 3-D array storing the signed distance of points in the voxel grid. 


The output of the function is a 2D array of size $K*11$, where each row stores the parameter of a superquadric $[\epsilon_1, \epsilon_2, a_x, a_y, a_z, euler_z, euler_y, euler_x, t_x, t_y, t_z]$.

### Preparing SDF from meshes
If you do not have SDF files but want to test the algorithm, we have a simple script [mesh2sdf_preparation](/mesh2sdf_preparation) to generate SDF from meshes. The script is based on this [package](https://github.com/wang-ps/mesh2sdf).
The mesh file will be first transformed to be watertight so that a valid SDF can be extracted.
To run the script, first run
```
pip install mesh2sdf
```
Then simply run
```
python3 mesh2sdf_convert.py $location of the mesh file$ --normalize --grid_resolution 100
```
where the mesh and sdf will be normalized within $[-1, 1]$ if you add `–-normalized`; and `--grid_resolution` specifies the resolution of the sdf (default $100$).
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
