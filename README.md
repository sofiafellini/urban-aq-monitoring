# Urban Air Quality Monitoring Network Optimization

This repository contains the MATLAB code associated with the paper:

> **A model-independent matrix method for optimizing urban air quality monitoring networks**
> Beatrice Carlini, Pietro Salizzoni, Luca Ridolfi, Sofia Fellini
> *Under review*, 2026

---

## Overview

The code implements a matrix-based methodology for optimizing the placement of air quality sensors in urban environments. The method operates on a source-receptor transport matrix **A**, which encodes pollutant propagation from potential emission sources to candidate sensor locations. It identifies optimal sensor configurations by maximizing two complementary metrics: the spatial extent of the detection basin and the precision of source identification.

The methodology is model-independent: the transport matrix **A** can be derived from any dispersion model (simplified parametric, CFD, statistical) or from observational data. In this repository, **A** is constructed using the network-based dispersion model of Fellini et al. (2019, 2020) as an illustrative example.

---

## Repository Structure

```
repository/
│
├── step01_receptor_pairs.m          ← main analysis: optimal receptor pairs
├── step02_receptor_pairs_basin.m    ← detection basin node lists for pairs
├── step03_receptor_triplets.m       ← extension to receptor triplets
│
├── FUNCTIONS/
│   ├── distances_threshold.m        ← shortest-path decay matrix
│   └── dijkstra_threshold.m         ← weighted Dijkstra algorithm
│
├── INPUT/
│   └── 5_<dir>_C0_Cth_10_<city>_25_03.mat   ← urban network data
│
└── ELABORATION/
    ├── step04_plot_basin.m          ← plot detection basin (Figures 3, 4, 8)
    ├── step05_plot_receptor_effectiveness.m  ← plot R_a map (Figures 5, 7, 9)
    └── PLOT_FUNCTIONS/
        ├── tight_subplot.m
        └── customcolormap_preset.m
```

---

## Pipeline

The scripts must be run in the following order:

**1. `step01_receptor_pairs.m`**
Computes the detection basin index B_(ab,n), the precision index P_ab, and the combined performance index R_ab for all receptor pairs. Identifies the optimal pair and classifies nodes in the detection basin by precision level (Pi = 0, Pi > 0, Pi = 1). Also computes the individual node performance score R_a.

**2. `step02_receptor_pairs_basin.m`**
Builds the full detection basin node list for every receptor pair. This output is required as input for step 3.

**3. `step03_receptor_triplets.m`** *(optional)*
Extends the analysis to receptor triplets. Computes R_abc for all triplet combinations and identifies the optimal triplet. Note: this step is computationally intensive for large urban networks.

**4. `step04_plot_basin.m`** *(run from ELABORATION/)*
Visualises the detection basin of the best receptor pair (Figures 3 and 4) or triplet (Figure 8). Set `receptor_mode = 'pairs'` or `receptor_mode = 'triplets'`.

**5. `step05_plot_receptor_effectiveness.m`** *(run from ELABORATION/)*
Visualises the spatial distribution of R_a across the urban network (Figures 5, 7, 9). Set `receptor_mode = 'pairs'` or `receptor_mode = 'triplets'`.

---

## Input Data

The `INPUT/` folder contains urban network data for two example cities:

- **Lyon** (France) — 748 nodes
- **Florence** (Italy) — 706 nodes

Each `.mat` file corresponds to one city and one wind direction, and contains the following variables:

| Variable | Description |
|---|---|
| `G_str` | Directed graph of the urban street network (MATLAB graph object) |
| `X_node_a`, `Y_node_a` | Node coordinates |
| `L_str_ord_a` | Link lengths [m] |
| `Ud_str_ord_a` | Mean wind speed projected on each link [m/s] |
| `H_str_ord_a` | Mean building height along each link [m] |
| `U_str_ord_a` | Mean street-canyon wind speed [m/s] |

To apply the method to a new city, construct a transport matrix **A** of size (n x n) where A(i,j) > 0 if source j influences receptor i, and replace the transport matrix construction block in `step01_receptor_pairs.m`. No other changes are required.

---

## Output Files

All output files are written as `.txt` in the working directory (for steps 1–3).

| File | Description |
|---|---|
| `M_<city>_<dir>_delta0<d>_original.txt` | Full pair metric matrix (unsorted) |
| `M_bestpair_<city>_<dir>_delta0<d>.txt` | Pair metric matrix sorted by R_ab |
| `v_certain_nodes_<city>_<dir>_delta0<d>.txt` | Basin nodes with Pi = 1 |
| `v_uncertain_<city>_<dir>_delta0<d>.txt` | Basin nodes with Pi > 0 |
| `v_zero_nodes_<city>_<dir>_delta0<d>.txt` | Basin nodes with Pi = 0 |
| `single_node_<city>_<dir>_delta0<d>.txt` | Individual node scores R_a (pairs) |
| `pairs_<city>_<dir>_basin_nodes.txt` | Basin node lists for all pairs |
| `M_triplets_<city>_<dir>_delta0<d>.txt` | Triplet metric matrix sorted by R_abc |
| `best_triplet_basin_<city>_<dir>_delta0<d>.txt` | Sub-basins of the best triplet |
| `triplets_single_node_<city>_<dir>_delta0<d>.txt` | Individual node scores R_a (triplets) |

---

## Requirements

- MATLAB R2019b or later
- MATLAB Graph and Network Algorithms (built-in)
- No additional toolboxes required

---

## Key Parameters

| Parameter | Description | Default |
|---|---|---|
| `delta` | Relative uncertainty in the transport matrix A | 0.1 |
| `city_names` | List of city labels | `{'Firenze', 'Parigi', 'New_York', 'Lione'}` |
| `wind_dirs` | Wind directions to analyse [deg] | `[360, 45, 90, 135, 180, 225, 270, 315]` |
| `receptor_mode` | `'pairs'` or `'triplets'` (plot scripts only) | `'pairs'` |

---

## License

The code is released under the [MIT License](LICENSE).
The input data in `INPUT/` are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Citation

If you use this code, please cite:

```
Carlini, B., Salizzoni, P., Ridolfi, L., Fellini, S. (2026).
A model-independent matrix method for optimizing urban air quality monitoring networks.
(Under review)
doi: [to be added upon publication]
```

---

## References

Fellini, S., Salizzoni, P., Soulhac, L., Ridolfi, L. (2019).
Propagation of toxic substances in the urban atmosphere: A complex network perspective.
*Atmospheric Environment*, 198, 291–301.

Fellini, S., Salizzoni, P., Ridolfi, L. (2020).
Centrality metric for the vulnerability of urban networks to toxic releases.
*Physical Review E*, 101, 032312.
