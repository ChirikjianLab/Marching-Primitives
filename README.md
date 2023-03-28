# Marching-Primitives: Shape Abstraction from Signed Distance Function
Project|Paper|Supplementary|[Arxiv](https://arxiv.org/abs/2303.13190)|[3D-Demos](/MATLAB/examples)|[Data](/MATLAB/data)

<img src="/examples/example.jpg" alt="example" width="600"/>

This repo provides the source code for the CVPR2023 paper:
> [**Marching-Primitives: Shape Abstraction from Signed Distance Function**](https://arxiv.org/abs/2303.13190 "ArXiv version of the paper.")  
> Weixiao Liu<sup>1,2</sup>, Yuwei Wu<sup>2</sup>, [Sipu Ruan](https://ruansp.github.io/)<sup>2</sup>, [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)<sup>2</sup>  
> <sup>1</sup> National University of Singapore, <sup>2</sup> Johns Hopkins University

## Update
**March 27th, 2023** - V0.1 basic version is online, including MATLAB implementation of the algorithm, data(SDF) preprocess script, and visualization tools.
Python implementation is planned in the coming updates.

## Abstraction
Representing complex objects with basic geometric primitives has long been a topic in computer vision. Primitive-based representations have the merits of compactness and computational efficiency in higher-level tasks such as physics simulation, collision checking, and robotic manipulation. Unlike previous works which extract polygonal meshes from a signed distance function (SDF), in this paper, we present a novel method, named Marching-Primitives, to obtain a primitive-based abstraction directly from an SDF. Our method grows geometric primitives (such as superquadrics) iteratively by analyzing the connectivity of voxels while marching at different levels of signed distance. For each valid connected volume of interest, we march on the scope of voxels from which a primitive is able to be extracted in a probabilistic sense and simultaneously solve for the parameters of the primitive to capture the underlying local geometry. We evaluate the performance of our method on both synthetic and real-world datasets. The results show that the proposed method outperforms the state-of-the-art in terms of accuracy, and is directly generalizable among different categories and scales.

## Implementation
Updating...

## Acknowledgement

## Citation
