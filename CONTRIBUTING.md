# Contributing

Some advice for first-time contributors to GraphsOptim.jl:

- Before pushing, format the code with `JuliaFormatter.format(GraphsOptim)`, otherwise the tests will fail.
- Whenever possible, try to write a version of the function that mutates an existing `JuMP.Model` and one that creates it, following the template given in `src/flow.jl`.
- Don't use any LaTeX in the docstrings, because it will not display nicely in the REPL. Instead, use Unicode symbols, and put the mathematical details in `docs/src/algorithms.md`.