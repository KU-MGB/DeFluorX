# Graph Report - .  (2026-07-28)

## Corpus Check
- 14 files · ~475,823 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1141 nodes · 2018 edges · 12 communities
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 47 edges (avg confidence: 0.53)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_00 00 run pipeline|00 00 run pipeline]]
- [[_COMMUNITY_00 01 Project Config|00 01 Project Config]]
- [[_COMMUNITY_00 02 Project Utils|00 02 Project Utils]]
- [[_COMMUNITY_00 03 Environment|00 03 Environment]]
- [[_COMMUNITY_01 Merge|01 Merge]]
- [[_COMMUNITY_02 Production|02 Production]]
- [[_COMMUNITY_03 Validation Figures|03 Validation Figures]]
- [[_COMMUNITY_04 Dendrogram|04 Dendrogram]]
- [[_COMMUNITY_05 TopN and PDB Preparation|05 TopN and PDB Preparation]]
- [[_COMMUNITY_06 Physics Validation|06 Physics Validation]]
- [[_COMMUNITY_07 MD QMMM Defluorination|07 MD QMMM Defluorination]]
- [[_COMMUNITY_Docs (README + deps)|Docs (README + deps)]]

## God Nodes (most connected - your core abstractions)
1. `main()` - 31 edges
2. `_echo()` - 31 edges
3. `main()` - 31 edges
4. `_generate_comprehensive_figures_impl()` - 28 edges
5. `process_single_job()` - 27 edges
6. `CFG` - 26 edges
7. `console_info()` - 26 edges
8. `prep_and_convert_phase()` - 24 edges
9. `run_mmgbsa_phase()` - 20 edges
10. `process_single_job()` - 19 edges

## Surprising Connections (you probably didn't know these)
- `CFG` --documents--> `Feasibility-weighted mechanistic tier ladder (Tier_1A..Tier_3/Decoy)`  [EXTRACTED]
  00_01_Project_Config_FAcDs.py → README.md

## Import Cycles
- None detected.

## Communities (12 total, 0 thin omitted)

### Community 0 - "00 00 run pipeline"
Cohesion: 0.22
Nodes (13): _cancel_schrodinger_jobs(), _ensure_jobserver(), _print_step_menu(), _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_FAcDs.sh script, _sync_staging_log() (+5 more)

### Community 1 - "00 01 Project Config"
Cohesion: 0.09
Nodes (8): 00_01_Project_Config_FAcDs.py (CFG SSOT), CFG, ===============================================================================, Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, Central Configuration Repository - FAcDs Pipeline.      Attribute prefix → secti, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant

### Community 2 - "00 02 Project Utils"
Cohesion: 0.03
Nodes (85): apply_figure_style(), atomic_write_csv(), auto_label_colour(), _box_diag(), calculate_angle(), calculate_burgi_dunitz(), calculate_dihedral(), calculate_flippin_lodge() (+77 more)

### Community 3 - "00 03 Environment"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 4 - "01 Merge"
Cohesion: 0.17
Nodes (15): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), main(), process_and_write(), Initialises a dual-handler logger (File + Console).     Log is saved beside the (+7 more)

### Community 5 - "02 Production"
Cohesion: 0.03
Nodes (130): Active_Site_Alignments/alignment stats CSV (residue mapping), 4_Prediction_Jobs/ (Boltz CIF models), 02_Production_FAcDs.py (upstream, authoritative rank), analyse_candidate_structure(), analyse_model_task(), analyse_pi_interactions(), analysis_worker_loop(), append_rows_to_csv() (+122 more)

### Community 6 - "03 Validation Figures"
Cohesion: 0.02
Nodes (256): analyse_conflicts (confidence-vs-tier conflict, hidden gems), calculate_pareto_fronts (non-dominated sorting), _xn__ensure_multimodel_variance_csv (per-model variance), 05_Figure_Descriptions.txt, 03_Figure_Enriched_Dataset.csv (figure columns, NOT rank source), generate_comprehensive_figures (figure suite driver), generate_ramachandran_figures (control backbone geometry), load_and_prep_data (reads ranked CSV, renames cols) (+248 more)

### Community 7 - "04 Dendrogram"
Cohesion: 0.16
Nodes (16): 04_Dendrogram_FAcDs.py (downstream consumer), clean_id(), console_info(), console_separator(), generate_phylogenies(), generate_upgma_newick(), get_kmer_counts(), main() (+8 more)

### Community 8 - "05 TopN and PDB Preparation"
Cohesion: 0.03
Nodes (114): _append_auxiliary_log(), _build_hd1_line(), check_prep_needed(), _check_residue_identity_guard(), cif_to_pdb_gemmi(), collect_best_cifs(), console_info(), console_separator() (+106 more)

### Community 9 - "06 Physics Validation"
Cohesion: 0.02
Nodes (221): _analysis_dir(), _avail_ram_gb(), _avg_dg(), _await_mmgbsa_jobserver(), _block_bootstrap_median_ci(), _boltzmann_mean_dg(), _build_msj(), _cancel_launched_jobs() (+213 more)

### Community 10 - "07 MD QMMM Defluorination"
Cohesion: 0.02
Nodes (156): _active_site_series(), _blockade_vec(), calculate_min_distance(), check_md_equilibration(), _collect_qsite_results(), compute_wm_csv_stats(), console_info(), console_qmm_ready() (+148 more)

### Community 11 - "Docs (README + deps)"
Cohesion: 0.01
Nodes (38): biopython 1.84 (sequence I/O), boltz 2.2.1 (Boltz-2 co-folding model), colabfold 1.6.1 (MSA generation), gemmi 0.6.5 (CIF to PDB conversion), hdbscan 0.8.43 + umap-learn 0.5.11 (embedding clustering), matplotlib 3.10.8 (figure rendering), MDAnalysis 2.9.0 (trajectory / geometry analysis), numpy 1.26.4 (numerics) (+30 more)

## Knowledge Gaps
- **57 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `load_and_prep_data (reads ranked CSV, renames cols)` (+52 more)
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Feasibility-weighted mechanistic tier ladder (Tier_1A..Tier_3/Decoy)` connect `Docs (README + deps)` to `00 00 run pipeline`, `00 01 Project Config`, `00 02 Project Utils`, `00 03 Environment`, `01 Merge`, `02 Production`, `03 Validation Figures`, `04 Dendrogram`, `05 TopN and PDB Preparation`, `06 Physics Validation`, `07 MD QMMM Defluorination`?**
  _High betweenness centrality (0.072) - this node is a cross-community bridge._
- **Why does `CFG` connect `00 01 Project Config` to `00 00 run pipeline`, `00 02 Project Utils`, `00 03 Environment`, `01 Merge`, `02 Production`, `03 Validation Figures`, `04 Dendrogram`, `05 TopN and PDB Preparation`, `06 Physics Validation`, `07 MD QMMM Defluorination`, `Docs (README + deps)`?**
  _High betweenness centrality (0.072) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `main()` (e.g. with `worker_task_wrapper()` and `analysis_worker_loop()`) actually correct?**
  _`main()` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 18 inferred relationships involving `_generate_comprehensive_figures_impl()` (e.g. with `_jf_volcano()` and `_jf_manhattan()`) actually correct?**
  _`_generate_comprehensive_figures_impl()` has 18 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _521 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `00 01 Project Config` be split into smaller, more focused modules?**
  _Cohesion score 0.09090909090909091 - nodes in this community are weakly interconnected._
- **Should `00 02 Project Utils` be split into smaller, more focused modules?**
  _Cohesion score 0.027771556550951846 - nodes in this community are weakly interconnected._