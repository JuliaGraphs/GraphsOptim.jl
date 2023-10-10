# Contributing

Some advice for first-time contributors to GraphsOptim.jl.

## Content

- For each algorithm, write a version of the function that mutates an existing `JuMP.Model` and another version that creates a new one.
- Add tests comparing the output to expected values

## Documentation

- Put the mathematical details in `docs/src/algorithms.md` (you can use LaTeX).
- Put the description of the arguments in the function docstring (you can use Unicode but no LaTeX)
- Build the docs locally to check that it works
- You can also preview the docs build at <https://juliagraphs.org/GraphsOptim.jl/previews/PR##> where "##" is the number of your PR

## Formalities

- Before pushing, format the code with `JuliaFormatter.format(GraphsOptim)`, otherwise the tests will fail.
