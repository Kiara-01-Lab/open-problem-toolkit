# Repository Guidelines

## Project Structure & Module Organization
- `src/LatticeReductionAlgorithms.jl`: Main module implementing GSO, LLL, MLLL, BKZ routines.
- `Project.toml` / `Manifest.toml`: Julia environment; keep synced with `--project`.
- `main.jl`: Minimal example invoking `BKZ_reduction!`. Treat as a runnable demo.
- Suggested: add tests under `test/` (see below).

## Build, Test, and Development Commands
- Install deps: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run demo: `julia --project=. main.jl`
- REPL dev: `julia --project=.`, then `using Revise, LatticeReductionAlgorithms` (if you use Revise).
- Run tests (when `test/` exists): `julia --project=. -e 'using Pkg; Pkg.test()'`
- Optional formatting: `julia -e 'using JuliaFormatter; format(".")'` (install with `Pkg.add("JuliaFormatter")`).

## Coding Style & Naming Conventions
- Use standard Julia style (4‑space indent, no hard tabs).
- Names: modules/types in `CamelCase` (e.g., `GSOData`); functions in `snake_case` (e.g., `size_reduce!`).
- Mutating functions end with `!`.
- Unicode identifiers are used (e.g., `μ`, `ν`, `B⃗`); be consistent and consider ASCII aliases in comments where helpful.
- Keep public API minimal and documented with docstrings `""" ... """`.

## Testing Guidelines
- Framework: Julia `Test` stdlib. Create `test/runtests.jl` and focused files like `test/test_lll.jl`.
- Example: `@test size(LLL_reduce(rand(Int, 5, 5), 0.75).B) == (5,5)`.
- Run with coverage: `julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'`.
- Prefer small, deterministic matrices; include edge cases (zero columns, ill‑conditioned bases).

## Commit & Pull Request Guidelines
- History shows short messages (e.g., “WIP”). Please switch to concise, imperative titles with context, or Conventional Commits (e.g., `feat(bkz): add block size param`).
- PRs should include: purpose, brief design notes, benchmarks if performance‑impacting, and checklists for tests/docs.
- Link related issues and annotate breaking changes clearly.

## Security & Configuration Tips
- Use project activation (`--project=.`) to avoid global env drift.
- Large integer/float operations can be costly; gate heavy examples behind `@info` and keep demos reproducible.
- Add `[compat]` bounds for dependencies and Julia version in `Project.toml` to ensure resolver stability.

