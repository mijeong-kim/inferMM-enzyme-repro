# inferMM-enzyme-repro

Reproduction materials for *Variance-Aware Estimation and Inference for
Michaelis–Menten Regression with Heteroscedastic and Clustered Enzyme-Kinetic
Data*.

## Revised-manuscript materials

The verified materials are in [`revision_20261125/`](revision_20261125/).
They contain the analysis scripts, raw inputs, replication-level reference
results, software-version records, table-generation code and verification
commands supplied with the revised manuscript.

This directory is an unchanged extraction of
`reproducibility_9078092_20261123_final.zip`, the numerical archive included
with the 25 November 2026 submission files. The archive date identifies the
reproduction-materials version; subsequent manuscript wording changes do not
alter these numerical materials.

The repository's earlier `scripts/`, `data/` and `results/` directories describe
the original analysis workflow. Use the versioned directory above for the
revised-manuscript calculations.

## Verification and execution

Start in the versioned directory and verify the 108 recorded file hashes:

```sh
cd revision_20261125
python3 verify_manifest.py
```

Follow [`README.txt`](revision_20261125/README.txt) for the reference R
environment and external dependencies. The supplied source package is
`inferMM` 0.1.2, which is available on CRAN. Install the supplied archive to
use the exact package build accompanying these materials:

```sh
Rscript -e "dir.create('library012', showWarnings=FALSE)"
R CMD INSTALL --library=library012 packages/inferMM_0.1.2.tar.gz
Rscript run.R all
Rscript verify.R
```

Individual analysis modules can also be run separately, as listed in
`README.txt`. Full simulations and bootstrap analyses can take substantial
time. Generated files are written separately from `expected/`.

The original simulation execution and timing records are retained, including
their original package-version metadata. The README documents the supplied
0.1.2 build and compatibility checks. Data sources, preprocessing, simulation
designs, inference assumptions and numerical-convergence qualifications are
also documented in the versioned materials.

## Package repository

The reusable R package is maintained at
[mijeong-kim/inferMM](https://github.com/mijeong-kim/inferMM).

## Archive identity

The unchanged reproduction ZIP has SHA-256:

```text
e15eed58a6f32b19b6e668e741145e86f80cb2730229c07aa28a8ba1f9a2a448
```

Manuscript PDFs, reviewer responses and internal review checklists are
maintained separately from this numerical reproduction repository.
