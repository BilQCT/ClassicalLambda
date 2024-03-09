# ClassicalLambda

## Overview
The central focus of this repository is the jupyter notebook [ClassicalLambda](ClassicalLambda.ipynb) which demonstrates the double description (DD) method described in arXiv:2312.10734 to obtain vertices of $\Lambda_2$, the 2-qubit $\Lambda$ polytope.

## Usage

The necessary packages are:

Polymake
GAP

You can choose to run the script [init.jl](./init.jl) to generate a new Julia project, which will automatically install the dependencies.

## Contents

### Simplicial distributions

There are several polytopes involved in our analysis. Many of them, such as $\text{NS}$, $\Lambda_2$, and $\text{MP}$, can be described using the framework of [simplicial distributions](https://quantum-journal.org/papers/q-2023-05-22-1009/). From the simplicial point of view, measurement scenarios are encoded as topological spaces using combinatorial objects known as simplicial sets. The script [twodim.jl](./lib/twodim.jl) is dedicated to tools for working with the two-dimensional simplicial sets that are used here. (A more comprehensive repository can be found at [SimpDist](https://github.com/okaygit/SimpDist).)

### Double description method

For highly degenerate polytopes the DD method (see e.g., [Fukuda-Prodon, 1995](https://link.springer.com/chapter/10.1007/3-540-61576-8_77)) is an effective tool for description conversion of a polyhedral cone from its $H$-representation to its $V$-representation, or vice versa. Our Julia implementation of the DD algorithm can be found in [doubledescription.jl](./lib/doubledescription.jl).

### Symmetries

We deal with highly symmetric polytopes, thus we can reduce by group action the number of vertices to a small number of representatives of each orbit. Currently we use the Julia interface to [Polymake](https://polymake.org/doku.php/start) via the [Oscar](https://www.oscar-system.org/) project to handle this part, which is in the script [symmetries.jl](./lib/symmetries.jl) However, eventually consideration of symmetries will be implemented directly via [GAP](https://www.gap-system.org/).






**Acknowledgments:** This work is supported by the Digital Horizon Europe project FoQaCiA, GA no.101070558. The authors also want to acknowledge support from the US Air Force Office of Scientific Research
under award number FA9550-21-1-0002.
