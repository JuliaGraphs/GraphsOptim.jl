# Benchmarks

Benchmarks of GraphsOptim.jl against other Julia and Python libraries (Graphs.jl, graspologic, NetworkX, SciPy), in terms of running time and solution quality.
The results are shown on the [Benchmarks](https://juliagraphs.org/GraphsOptim.jl/dev/benchmarks/) page of the documentation, which also describes the methodology.

| File | Content |
| :--- | :--- |
| `run.jl` | entry point, runs the suites and draws the figures |
| `suites/*.jl` | one suite per family of problems: instances, methods and metrics |
| `common.jl` | shared tools: timing, calls to Python, CSV files |
| `plots.jl` | figures, saved in `docs/src/assets/benchmarks/` |
| `results/*.csv` | raw results |

To update the results and figures, run from the root of the repository:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/run.jl                    # all suites
julia --project=benchmark benchmark/run.jl cliques            # a single suite
julia --project=benchmark benchmark/run.jl --plots-only       # only redraw the figures
```

The Python libraries are installed in `benchmark/.CondaPkg` by [CondaPkg.jl](https://github.com/JuliaPy/CondaPkg.jl), and called through [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl).

To add a suite, create `suites/<name>.jl` with a function `run_<name>()` writing its CSV files, add a `plot_<name>()` function to `plots.jl`, register both in `SUITES` in `run.jl`, and add a section to `docs/src/benchmarks.md`.
