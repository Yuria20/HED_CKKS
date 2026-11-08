# HED-CKKS: Efficient Homomorphic Edit Distance Computation

This repository contains a proof-of-concept implementation of **HED-CKKS**, accompanying the paper *"HED-CKKS: Efficient Homomorphic Edit Distance Computation"*. It is built on [Lattigo v6](https://github.com/tuneinsight/lattigo).

> **Note:** This code is a proof-of-concept and is intended for research and reproducibility purposes only. It is not production-ready.

## Repository Structure

```
banded_HED/
└── lattigo/
    └── examples/
        └── singleparty/
            └── ckks_bootstrapping/
                └── slim/    # HED-CKKS implementation and experiment code
```

## Requirements

- Go 1.23.3 (linux/amd64)
- Memory: **128 GB RAM**
- Lattigo v6 (included in this repository under `lattigo/`)

## Getting Started

```bash
git clone https://github.com/Yuria20/banded_HED.git
cd banded_HED
```

## Reproducing the Experiments

The code for all experiments reported in the paper is located in `lattigo/examples/singleparty/ckks_bootstrapping/slim`. Running `go run .` in that directory reproduces the results.

```bash
cd lattigo/examples/singleparty/ckks_bootstrapping/slim
go run .
```

Experiments are run single-threaded for string lengths m = 8, 128, and 256. Timings depend on the machine, so absolute numbers may differ from those reported in the paper.
