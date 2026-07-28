# Boltz-2_Run_EXAMPLE_Controls

An example run directory, retained in the repository so that the on-disk layout of a completed
FAcDs / PFAS-27 pipeline run can be inspected directly. It is a controls-only illustration, not a
production run: it carries real, viewable data for the control cases while every other stage is
represented by its empty folder skeleton.

## Purpose

A production run screens 2,150 FAcD variants against 27 PFAS ligands (58,050 complexes) and reaches
tens of gigabytes. The best-complex CIF structures alone occupy roughly 12 GB, and the master and
ranked tables are approximately 350 MB each, so a complete run cannot be committed within GitHub's
file and repository limits. This directory therefore preserves the folder structure with a small,
representative slice of genuine output.

## Control cases

Six control complexes are included, formed from two reference proteins and three short-chain
substrates:

- Proteins: 3R3U (the RPA1163 / 3R3U crystal fluoroacetate dehalogenase) and DeHa4 (the functional
  positive control).
- Ligands: fluoroacetate (FA), difluoroacetate (DFA) and trifluoroacetate (TFA).

These correspond to the six `0000000_*` prediction jobs.

## Directories that contain data

Data is present only where it remains well within GitHub limits:

- `1_Boltz2_Production/`
  - `1_Input_Data/` - input protein FASTA and PFAS SMILES roster.
  - `2_Boltz2_YAML_Configs/` - the six control Boltz-2 job configurations.
  - `3_Sequence_Reference_Data/Active_Site_Alignments/` - the 3R3U and DeHa4 control alignments.
  - `4_Prediction_Jobs/` - the six control prediction jobs, each with its co-folded models, best
    complex and interaction and summary files.
  - `5_Boltz2_FAcDs_Master_*.csv` and `6_Boltz2_FAcDs_Ranked_*.csv` - the master and ranked tables,
    restricted to the six control rows.
- `2_Best_Complexes_CIFs/` - the six best control complex CIF structures.

Heavy, regenerable intermediates (multiple-sequence-alignment search caches and the `.npz`
confidence tensors) are omitted from the prediction jobs to keep the example compact.

## Directories retained as empty skeletons

Every remaining stage is present with its folder structure but without data, so the overall layout
is visible without exceeding size limits:

- `0_FAcDs_Pipeline_Logs/`
- `3_Validation_Figures/` - analysis data and the seven thematic figure folders.
- `4_Dendrogram/` - the per-tier interactive phylogeny applications.
- `5_TopN_and_Preparation/` - converted and prepared PDBs, comparative analysis, ligand ESP charges.
- `6_Physics_Validation/` - prepared proteins, ESP-charged complexes, WaterMaps, system builder,
  MD simulations, and the analysis folder (Prime MM-GBSA and defluorination).
- `7_MD_Thermodynamics_Results/`

Empty folders are held in version control by `.gitkeep` placeholders.

## Reproducing a full run

The complete pipeline resides in the repository root and runs in sequence from
`00_00_run_pipeline_FAcDs.sh` through `07_MD_QMMM_Defluorination_FAcDs.py`.
