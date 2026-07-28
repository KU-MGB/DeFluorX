# FAcDs / PFAS-27 - Example run (controls only)

This is a **test / example** Boltz-2 run kept in the repository so anyone can see how a
completed pipeline run is organised on disk. **It is not a full production run.**

## Why only a slice of the data is here

The real production run screens **2,150 FAcD variants x 27 PFAS = 58,050 complexes** and
grows to tens of gigabytes (the best-complex CIFs alone are ~12 GB, the master/ranked CSVs
~350 MB each). That far exceeds GitHub's file and repository limits, so a full run cannot be
committed. To still show the folder layout with real, inspectable data, this example keeps
**only the control cases**:

- **3R3U** (the RPA1163 / 3R3U crystal fluoroacetate dehalogenase) and
- **DeHa4** (the functional positive control)

each paired with **fluoroacetate (FA)**, **difluoroacetate (DFA)** and **trifluoroacetate
(TFA)** - the six `0000000_*` control jobs.

## What contains data

Data is present **only** in these two directories (kept small, well under GitHub limits):

- **`1_Boltz2_Production/`**
  - `1_Input_Data/` - the input FASTA + PFAS SMILES roster
  - `2_Boltz2_YAML_Configs/` - the 6 control Boltz-2 job configs
  - `3_Sequence_Reference_Data/Active_Site_Alignments/` - the 3R3U + DeHa4 control alignments
  - `4_Prediction_Jobs/` - the **6 control prediction jobs** (co-folded structures, best complex,
    interaction + summary files; heavy regenerable binaries such as MSA search caches and the
    `.npz` confidence tensors are omitted to keep the example light)
  - `5_Boltz2_FAcDs_Master_*.csv` and `6_Boltz2_FAcDs_Ranked_*.csv` - the master and ranked
    CSVs **filtered to the control rows only**
- **`2_Best_Complexes_CIFs/`** - the 6 best control complex CIF structures

## What is empty

Every other phase folder is present as an **empty skeleton** so the directory structure is
visible, but carries no data (its real outputs would exceed GitHub limits):

- `0_FAcDs_Pipeline_Logs/`
- `3_Validation_Figures/` (analysis data + the 7 figure folders)
- `4_Dendrogram/` (per-tier interactive apps)
- `5_TopN_and_Preparation/` (converted / prepared PDBs, comparative analysis, ESP charges)
- `6_Physics_Validation/` (prepared proteins, ESP complexes, WaterMaps, system builder,
  MD simulations, analysis + Prime-MMGBSA / defluorination)
- `7_MD_Thermodynamics_Results/`

(`.gitkeep` files hold the empty folders in git.)

To reproduce a real run, use the pipeline in the repository root (`00_00_run_pipeline_FAcDs.sh`
through `07_MD_QMMM_Defluorination_FAcDs.py`).
