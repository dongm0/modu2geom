# Getting Optimal Modular from node Connection using minimal-generarion-tree and harmonic-deformation

## Abstract

在本文中，我们提出了一种从六面体模板的连接方式得到高单元质量模板的方法。我们的方法基于单元堆积的思想。我们在模板上遵照他的连接方式逐步加入六面体单元，并在每次加入之后优化，这样我们就能得到理想的模板。我们的方法分为两步，首先我们列举了模板加入单元的所有形式，并以此构建了生成生成树（这里怎么写再想想）。然后我们用一种xx的deformation来确保每次加入新单元的时候很顺利，而且不会存在几何上的不合法性。

In this paper, we proposed a method getting optimal mesh from merely its node connection. Our method is based on the idea of element stacking. We add one hex element to a modular following its connection step by step and optimize the modular in each iteration, in this way we get the ideal modular. Out method contains two main steps. First, we enumerates different cases of adding one hex to modular, and we construct a generation tree, which ensures the deformation in the latter step less complicated. Then we use flipless harmonic-deformation to make adding new element convenient and without tangle. This method can be applied in mesh optimization, such as untangle. Using this method, we construct the 17-18 transformation modular.

hispapershowsthatconstraintprogrammingtechniquescansuccessfullybeusedtosolvechallenginghex-meshingproblems.Schneiders'pyramidisasquare-basedpyramidwhosefacetsaresubdividedintothreeorfourquadranglesbyaddingverticesatedgemidpointsandfacetcentroids.Inthispaper,weprovethatSchneiders'pyramidhasnohexahedralmesheswithfewerthan18interiorverticesand17hexahedra,andintroduceavalidmeshwith44hexahedra.Wealsoconstructthesmallestknownmeshoftheoctagonalspindle,with40hexahedraand42interiorvertices.Theseresultswereobtainedthroughageneralpurposealgorithmthatcomputesthehexahedralmeshesconformaltoagivenquadrilateralsurfaceboundary.ThelowerboundforSchneiders'pyramidisobtainedbyexhaustivelylistingthehexahedralmeshesawithupto17interiorverticesandwhichhavethesameboundaryasthepyramid.Our44-elementmeshisobtainedbymodifyingapriorsolutionwith88hexahedra.Thenumberofelementswasreducedusinganalgorithmwhichlocallysimplifiesgroupsofhexahedra.Giventheboundaryofsuchagroup,ouralgorithmisusedtofindameshofitsinteriorthathasfewerelementsthantheinitialsubdivision.Theresultingmeshisuntangledtoobtainavalidhexahedralmesh.

## Key Words

## Introduction

Triangle meshes remain a predominant representation of 3D sur-faces. When the complexity of a given mesh exceeds computationalresources, we rely on mesh simplification methods to remove ver-tices, edges, and faces. In rendering, efficient simplification methodscan dramatically reduce the complexity of a mesh without affectingits appearance. It is tempting to repurpose appearance-preservingsimplification methods for other geometry processing tasks.

Unfortunately, appearance-based methods do not preserve thespectral properties of the important differential operators upon whichmuch of modern geometry processing is built (see Figure1). As aresult, solutions computed on such a coarse mesh can be incorrector misleading. Alternatively, previous coarsening methods that dopreserve spectral properties work purely algebraically on the oper-ator matrices and do not produce a geometric mesh (see Figure2).The lack of a mesh limits the use of coarsening in many downstreamgeometry processing tasks.We present the first mesh simplification method intentionallydesigned to preserve spectral properties. 

We propose adapting thestandard greedy edge-collapse mesh-simplification algorithm with anovel cost function that measures spectral preservation of a given op-erator (e.g., the cotangent Laplacian). Unlike algebraic methods thatdirectly output a reduced operator (i.e., matrix), our method outputsa manifold triangle mesh with 3D vertex positions. Reconstructingthe operator on the output mesh will preserve both the eigenvaluesand eigenvectors of the operator on the input mesh.

Confirmed by a series of experiments, our method preserves spec-tral properties nearly as well as purely algebraic methods, whilestill outputting an embedded mesh like standard simplification algo-rithms. We demonstrate our approach’s effectiveness for geodesicdistance approximation and functional maps correspondence

在当下，六面体网格在仿真计算中仍具有相当优势性的地位，人们用各种各样的方法去得到一个高质量的六面体网格，利用预先计算好的网格模板是其中一种。模板可以用以确定划分时候的奇异graph，也可以在局部优化的时候做替换（ALM，模板（比利时人））。这些任务中都需要从拓扑到几何的转换，通常我们有很多方法来得到优化的网格，例如（能量函数法，插值法，光滑法等等），但是对于拓扑质量低的模板，这些方法往往很难得到高质量无扭的网格（怎么很难奏效的说清除），在本文中我们提出了一套流程，用以解决这个问题。（转化一下，看看第一个做xx问题的人是怎么做的）

写网格优化，很多情况下，我们只知道拓扑链接，希望得到高质量的几何位置，
st  
在当下，六面体网格在仿真计算中仍具有相当优势性的地位，人们用各种各样的方法去得到一个高质量的六面体网格。

然而，很多情况下我们只知道拓扑连接，需要从拓扑链接中得到好的几何位置。已有的方法可以解决部分问题，如基于能量函数的方法，插值法，光滑法等等

我们的流程基于单元堆积的思想。在一个已有部分单元集合的模板上，每次生长一个单元，在这个过程中保证它高质量以及合法性（untangle），逐步的完成整个模板的生成。整个过程中需要

