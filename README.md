# Open Problem Toolkit

### Test how secure tomorrow's encryption really is.

Quantum computers will eventually break today's encryption. **Lattice-based cryptography** is the leading replacement — but how well does it actually perform? This toolkit lets you build, test, and benchmark the core algorithms so you can find out.

> 量子コンピュータは現在の暗号を破る可能性があります。**格子暗号**はその有力な代替手段です。本ツールキットでは、格子暗号のコアアルゴリズムを構築・テスト・ベンチマークできます。

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)]()
[![Open Source](https://badgen.net/badge/open/source/green)]()

---

## Packages

| Package | What it does | Status |
|---------|-------------|--------|
| **LatticeBasedCryptography.jl** | Encrypt & decrypt using lattice schemes | 🟡 In progress |
| **LatticeReductionAlgorithms.jl** | Break lattices with reduction algorithms (LLL, KZ) | 🟡 In progress |
| **HomomorphicCryptography.jl** | Compute on encrypted data without decrypting | 🟡 In progress |
| **LibOQS.jl** | Julia bindings for [liboqs](https://openquantumsafe.org/) (NIST PQC algorithms) | 🟡 In progress |
| **PQCValidator.jl** | Validate post-quantum parameter security | 🔵 Early |
| **ZKPValidator.jl** | Validate zero-knowledge proof constructions | 🔵 Early |

---

## Quick Start

```julia
using Pkg
Pkg.add(url="https://github.com/Kiara-01-Lab/open-problem-toolkit")
```

```julia
# Example: Build a ring-based lattice and reduce its basis
using LatticeReductionAlgorithms

L = random_lattice(dim=64)
reduced = lll_reduce(L)
println("Reduction ratio: ", reduction_quality(reduced))
```

> ⚠️ API is under active development. Expect breaking changes before v1.0.

---

## Why This Exists

Post-quantum crypto is mostly studied in theory. Real-world performance data is scarce.

This toolkit provides **empirical evidence**: how fast these algorithms run, where they break down, and which parameter choices actually matter — benchmarked against the [Lattice Challenge](https://www.latticechallenge.org/).

Results are published as open preprints on [arXiv](https://arxiv.org/) with full datasets.

---

## Roadmap

| Phase | Milestone | Target |
|-------|-----------|--------|
| 1 | Lattice construction + LLL reduction | Sep–Oct |
| 2 | KZ reduction + first preprint | Nov |
| 3 | Modular solver framework v2 | Dec–Jan |
| 4 | Large-scale benchmarks + scalability analysis | Feb–Mar |

**Long-term goal:** A reference platform where anyone can build, test, and break post-quantum encryption with full transparency.

---

## Contributing

1. Fork → branch → commit → PR.
2. See [`COMMUNITIES.md`](COMMUNITIES.md) for related open-science communities.

All skill levels welcome — math, code, docs, or ideas.

---

## License

Apache 2.0 — free to use, modify, and distribute.
