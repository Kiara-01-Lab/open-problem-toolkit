# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

LatticeReductionAlgorithms.jl implements lattice basis reduction algorithms in Julia, including:
- **GSO (Gram-Schmidt Orthogonalization)**: Core data structure and operations
- **LLL (Lenstra-Lenstra-Lovász)**: Primary lattice reduction algorithm
- **MLLL (Modified LLL)**: Handles zero vectors and edge cases
- **BKZ (Block Korkine-Zolotarev)**: Advanced block-wise reduction using enumeration
- **ENUM_reduce**: Shortest vector enumeration for BKZ

The implementation uses exact integer arithmetic for basis vectors while performing floating-point GSO computations.

## Development Commands

### Environment Setup
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Enter REPL with project activated
julia --project=.
```

### Running Code
```bash
# Run the demo (main.jl)
julia --project=. main.jl

# Interactive development with Revise
julia --project=.
# Then: using Revise, LatticeReductionAlgorithms
```

### Testing
```bash
# Run tests (when test/ directory exists)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run with coverage
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
```

### Code Formatting
```bash
# Format code (requires JuliaFormatter)
julia -e 'using JuliaFormatter; format(".")'
```

## Architecture

### Core Data Structure: GSOData

Located in `src/LatticeReductionAlgorithms.jl:9-40`, `GSOData` maintains the Gram-Schmidt orthogonalization state:
- `B::Matrix{I}`: Integer basis matrix (column vectors)
- `B⃗::Vector{F}`: Squared norms of GS vectors (||b*ᵢ||²)
- `Q::Matrix{F}`: Orthogonalized basis (GS vectors scaled by original R diagonals)
- `R::Matrix{F}`: Size-reduction coefficients μᵢⱼ (upper triangular)

The constructor computes QR factorization and carefully handles zero/near-zero vectors by setting their row to represent no coupling.

### Algorithm Flow

1. **LLL_reduce** (`src/LatticeReductionAlgorithms.jl:171-192`):
   - Main reduction loop: size-reduce each basis vector, check Lovász condition
   - If condition fails: swap basis vectors via `gsoupdate!`, backtrack
   - Parameter δ must satisfy 0.25 < δ < 1

2. **MLLL_reduce!** (`src/LatticeReductionAlgorithms.jl:194-331`):
   - In-place variant handling zero vectors
   - Incrementally builds GSO during reduction
   - Swaps zero vectors to end, reduces dimension z
   - Called iteratively until no zero vectors remain

3. **BKZ_reduction!** (`src/LatticeReductionAlgorithms.jl:348-385`):
   - Initial LLL reduction
   - For each block [k, k+β-1]: find shortest vector via `find_svp_by_enum`
   - Insert found vector and re-reduce with MLLL
   - Terminates after n-1 rounds without improvement

4. **ENUM_reduce** (`src/LatticeReductionAlgorithms.jl:57-107`):
   - Pruned enumeration for shortest lattice vector
   - Uses branch-and-bound with radius constraints R²
   - Returns coefficient vector and success flag
   - Core of BKZ's SVP oracle

### Key Implementation Details

- **Mutating conventions**: Functions ending in `!` modify their arguments
- **Unicode identifiers**: μ (mu), ν (nu), B⃗ (B-vector) used throughout
- **Type parameters**: `GSOData{I<:Integer, F<:Real}` separates exact/floating types
- **Zero vector handling**: `iszerovec` checks approximate equality; MLLL swaps to end
- **Size reduction**: `partial_size_reduce!` subtracts integer multiple to minimize μᵢⱼ

### Example Usage (main.jl)

The demo constructs a 40×40 basis from an SVP challenge, extracts a 4×4 submatrix, and runs BKZ with block size β=3 and δ=0.75.

## Julia Conventions

- **Naming**:
  - Types/modules: `CamelCase` (e.g., `GSOData`)
  - Functions: `snake_case` (e.g., `size_reduce!`)
  - Mutating functions: suffix `!`
- **Indexing**: 1-based (Julia standard)
- **Broadcasting**: Use `.=` and `.+` for element-wise operations
- **Views**: `@view` macro avoids allocation for slices
- **Macros**: `@assert`, `@inbounds`, `@info` for assertions, optimization, logging

## Testing Approach

When creating tests (in `test/runtests.jl`):
- Use small deterministic matrices (e.g., 5×5 random integer)
- Verify structural properties: basis dimensions, orthogonality bounds
- Test edge cases: zero columns, ill-conditioned bases
- Example: `@test size(LLL_reduce(rand(Int, 5, 5), 0.75).B) == (5,5)`

## Performance Considerations

- Large integer arithmetic (BigInt in ENUM_reduce) is costly
- BKZ complexity grows exponentially with block size β
- `main.jl` processes 40-dimensional SVP challenge lattices; full BKZ runs can be very slow
- Consider smaller examples or lower β for quick iteration

## Dependencies

- **LinearAlgebra**: QR factorization, dot products, norms
- **OffsetArrays**: 0-indexed arrays in ENUM_reduce (r::OffsetVector)

See `Project.toml` for version constraints.
