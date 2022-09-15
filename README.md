# OscoNet.jl: inferring oscillatory gene networks

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://burtonjosh.github.io/OscoNet.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://burtonjosh.github.io/OscoNet.jl/dev/)
[![Build Status](https://github.com/burtonjosh/OscoNet.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/burtonjosh/OscoNet.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/burtonjosh/OscoNet.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/burtonjosh/OscoNet.jl)

A Julia implementation of OscoNet[^Cutillo2020], a method for detecting oscillatory gene networks from snapshot single cell data. For more details, see the [paper](https://doi.org/10.1186/s12859-020-03561-y).

# Installation

Open the Julia REPL by typing the following in the command line:
```bash
$ julia
```

Enter the Pkg REPL by typing `]` from the Julia REPL.
```
julia> ]
```

To install OscoNet, type the following into the Pkg REPL

```
(@v1.6) pkg> add "https://github.com/burtonjosh/OscoNet.jl.git"
```

To exit the Pkg REPL either press backspace or ^C (CTRL+C).

You can now start using the OscoNet package,

```
julia> using OscoNet
```

# Usage

Basic usage of the package:

```julia
using OscoNet
data, _, _ = simulate_data();
n_permutations = 100;

adjacency_matrix, _, cost = bootstrap_hypothesis_test(data, n_permutations);

gene_names = ["gene_$i" for i = 1:20];

edge_network = create_edge_network(adjacency_matrix, 1 ./ cost, gene_names)
```

For more details and examples, refer to the package [documentation](https://burtonjosh.github.io/OscoNet.jl).

[^Cutillo2020]: Luisa Cutillo, Alexis Boukouvalas, Elli Marinopoulou, Nancy Papalopulu & Magnus Rattray (2020).
    OscoNet: inferring oscillatory gene networks.
    BMC Bioinformatics: [10.1186/s12859-020-03561-y](https://doi.org/10.1186/s12859-020-03561-y).
    [Code](https://github.com/alexisboukouvalas/OscoNet)