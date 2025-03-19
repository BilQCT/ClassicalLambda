# ClassicalLambda

## Overview
The central focus of this repository is the jupyter notebook [ClassicalLambda](ClassicalLambda.ipynb) which demonstrates the double description (DD) method described in arXiv:2312.10734 to obtain vertices of $\Lambda_2$, the 2-qubit $\Lambda$ polytope. In particular, we characterize vertices of $\Lambda_2$ that can be described by a noncontextual hidden variable model. For this we take advantage of the fact that for all $n\in\mathbb{N}$ we have that $\Lambda_n\subset \Lambda_n^{(l)}$, where $\Lambda_n^{(l)} \subset \mathbb{R}^{4^n-1}$ can be identified with the nonsignaling polytope for the $n$-partite Bell scenario with three binary outcome observables per party. Certain vertices of $\Lambda_n^{\ell}$, the so-called deterministic vertices, provide a noncontextual hidden variable model (ncHVM) and the corresponding classical (or Bell) polytope is defined as the convex hull of deterministic vertices. In this repository we study the vertices $A_\alpha \in \Lambda_2$ that admit an ncHVM in the sense that $A_\alpha $ is in the classical polytope of the bipartite $(2,3,2)$ Bell scenario.

## Usage

The necessary packages are:

Polymake

GAP

Combinatorics

Nemo v0.22

You can choose to run the script [init.jl](./init.jl) to generate a new Julia project, which will automatically install the dependencies.

## Contents

### Simplicial distributions

There are several polytopes involved in our analysis. Many of them, such as $\text{NS}$, $\Lambda_2$, and $\text{MP}$, can be described using the framework of [simplicial distributions](https://quantum-journal.org/papers/q-2023-05-22-1009/). From the simplicial point of view, measurement scenarios are encoded as topological spaces using combinatorial objects known as simplicial sets. The script [twodim.jl](./lib/twodim.jl) is dedicated to tools for working with the two-dimensional simplicial sets that are used here. (A more comprehensive repository can be found at [TwoDim](https://github.com/BilQCT/TwoDim).)

For a simplicial scenario $(X,N\mathbb{Z}_d)$ we have that the set of simplicial outcome maps $r:X\to N\mathbb{Z}_d$ are determined by the $2$-skeleton. These outcome maps can be computed in GAP which can then in turn be converted into deterministic vertices. The latter can be used as an input to compute the facets of the corresponding Bell polytope. This is implemented in the script [bell.jl](./lib/bell.jl).

### Double description method

For highly degenerate polytopes the DD method (see e.g., [Fukuda-Prodon, 1995](https://link.springer.com/chapter/10.1007/3-540-61576-8_77)) is an effective tool for description conversion of a polyhedral cone from its $H$-representation to its $V$-representation, or vice versa. Our notebook [ClassicalLambda.ipynb](./ClassicalLambda.ipynb) gives a walkthrough of how we use the main proposition of the DD method to obtain vertices of a polytope $\overline{\text{MP}}$. Functions needed for this implementation can be found in [doubledescription.jl](./lib/doubledescription.jl).

### Symmetries

We deal with highly symmetric polytopes, thus we can reduce by group action the number of vertices to a small number of representatives of each orbit. Currently we use the Julia interface to [Polymake](https://polymake.org/doku.php/start) via the [Oscar](https://www.oscar-system.org/) project to handle this part, which is in the script [symmetries.jl](./lib/symmetries.jl) However, eventually consideration of symmetries will be implemented directly via [GAP](https://www.gap-system.org/).






**Acknowledgments:** This work is supported by the Digital Horizon Europe project FoQaCiA, GA no.101070558. The authors also want to acknowledge support from the US Air Force Office of Scientific Research
under award number FA9550-21-1-0002.
