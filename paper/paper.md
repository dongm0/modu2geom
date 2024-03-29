# Getting Optimal Modular from node Connection using minimal-generarion-tree and harmonic-deformation

## Abstract

在本文中，我们提出了一种从六面体模板的连接方式得到高单元质量模板的方法。我们的方法基于单元堆积的思想。我们在模板上遵照他的连接方式逐步加入六面体单元，并在每次加入之后优化，这样我们就能得到理想的模板。我们的方法分为两步，首先我们列举了模板加入单元的所有形式，并以此构建了生成生成树（这里怎么写再想想）。然后我们用一种xx的deformation来确保每次加入新单元的时候很顺利，而且不会存在几何上的不合法性。

In this paper, we proposed a method getting optimal mesh from merely its node connection. Our method is based on the idea of element stacking. We add one hex element to a modular following its connection step by step and optimize the modular in each iteration, in this way we get the ideal modular. Out method contains two main steps. First, we enumerates different cases of adding one hex to modular, and we construct a generation tree, which ensures the deformation in the latter step less complicated. Then we use flipless harmonic-deformation to make adding new element convenient and without tangle. This method can be applied in mesh optimization, such as untangle. Using this method, we construct the 17-18 transformation modular.

## Key Words

## Introduction

At present, in the field of simulating calculation, hexahedral/quad mesh has superiority to tetrahedral/tri mesh. Comparing to Tetrahedral mesh, hex mesh requires much less space to store, and they usually convergence faster than tet mesh. Owing to its properties, people put forward many methods to get hex mesh of high quality. In some cases, wo only know the connections of the mesh, sometimes may have boundary positions, and we want to acquire node positions that ensures mesh quality as high as possible. For example, ()

Many existing methods can partly solve these questions, such as smoothing, interpolating, and energy-based optimizing methods. Unfortunately, when vertex connections are of low quality, for example, when the input mesh is highly-unstructured, containing lots of singularties, it's difficult to get high-quality and untangled mesh. Some distortion energy may have barrier layers to prevent tangling, but those barriers need an initial untangle mesh. To solve those problems, we put forward a procedure based on thought of mesh accumulation.

Our main approach is to add hexahedral cells step by step, during which we make sure the mesh is untangled and with high quality. We take two steps to accomplish this. First we classified several different situations of adding one cell, and we propose a method to estimate the relationship between cell valence and deform in every step. Using this method, we construct a most suitable sequence to generate the whole mesh. 

## background

## our method

### 3.1 initial untangled mesh



### 3.2 