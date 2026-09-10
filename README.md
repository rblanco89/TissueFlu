# TissueFlu

Simulation of influenza spread in the lungs of mice: a tissue of individual
epithelial cells on the GPU, exchanging virions and interferon (IFN) with their
neighbours. One time step is one minute.

## Build

```
make                       # needs the CUDA toolkit and a supported gcc
make ARCH=sm_75            # for a different GPU (sm_86 = RTX 30xx / A2000)
```

The host compiler follows the CUDA version, so one Makefile serves both
machines: CUDA 13 or newer uses `gcc-15`, older CUDA uses the newest of
`gcc-14`/`gcc-13`/`gcc-12`/`gcc` it finds. `make HOST_COMPILER=...` overrides.

## Run

```
./tissueFlu --config lung.conf --structure cells_normalized_z0.csv --output results
```

* `--config` — parameter file, one `name = value` per line (`#` starts a comment)
* `--structure` — cell positions, a CSV with a header line then `x,y,z` per cell
* `--output` — directory for the result files (created if missing)

The neighbour list is read from `<structure>_neighbors_r<R>.csv` next to the
structure file. If that file is missing it is computed and saved there, so the
search happens only once per structure and radius. Delete the file to force a
rebuild (for instance after changing `neighRadius`).

## Output

| file | contents |
|---|---|
| `tissueState_<n>.csv` | tissue-wide time course of replicate `n` (always written, one per replicate) |
| `auc_results.csv` | area under the viral-load and IFN curves, per replicate |
| `cellState.csv` | time course of the last seeded cell (single replicate) |
| `snapshots.xyz` | extended-XYZ frames for Ovito, if `printSnapshots` is 1 or 2 |

## Source layout

| file | contents |
|---|---|
| `src/model.cu` | **the biological model**: cell state machine, transport, CTL response |
| `src/main.cu` | command line, initial condition, the time-step loop |
| `src/tissue.cu` | reading the structure, neighbour lists, GPU allocation |
| `src/params.cu` | the configuration file and its defaults |
| `src/output.cu` | measurements and result files |
| `src/headers.h` | data types shared by all of the above |
| `src/ranNumbers.h` | host random-number generators (Numerical Recipes) |

Averages over replicates are left to the analysis scripts.

## Random seeds and reproducibility

* `ranSeed` — replicate `n` uses `ranSeed + n - 1` for everything that varies
  between replicates: where the inoculum lands, infected-cell lifetimes, and
  every stochastic event during the run.
* `tissueSeed` (optional) — if set, the non-permissive cells are drawn from this
  seed alone, so the same cells are non-permissive in every replicate and in
  every run that uses the same `tissueSeed` and structure file. If it is absent,
  a new set of non-permissive cells is drawn for each replicate from `ranSeed`.

A run is fully determined by its seeds: each cell draws its random numbers from
(seed, cell, time step) rather than from a shared stream, so repeating a run
gives byte-identical results, and the outcome does not depend on how the cells
happen to be ordered internally or on how the GPU schedules its threads.
