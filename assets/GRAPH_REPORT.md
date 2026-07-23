# Graph Report - .  (2026-07-23)

## Corpus Check
- 14 files · ~257,881 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1253 nodes · 2545 edges · 74 communities (70 shown, 4 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 166 edges (avg confidence: 0.72)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Ranking & Pareto Figures|Ranking & Pareto Figures]]
- [[_COMMUNITY_Step-06 Orchestration & Logging|Step-06 Orchestration & Logging]]
- [[_COMMUNITY_Tier Quality Figures|Tier Quality Figures]]
- [[_COMMUNITY_Merge & Alignment Cache|Merge & Alignment Cache]]
- [[_COMMUNITY_Figure Panel Helpers|Figure Panel Helpers]]
- [[_COMMUNITY_MD Heartbeat & Progress|MD Heartbeat & Progress]]
- [[_COMMUNITY_Figure Style & Console SSOT|Figure Style & Console SSOT]]
- [[_COMMUNITY_Validation Figure Suite (03)|Validation Figure Suite (03)]]
- [[_COMMUNITY_CFG Central Configuration|CFG Central Configuration]]
- [[_COMMUNITY_Report Manager & Data Loading|Report Manager & Data Loading]]
- [[_COMMUNITY_Step-07 QMMM Engine|Step-07 QM/MM Engine]]
- [[_COMMUNITY_Step-05 Prep & Interaction Diagrams|Step-05 Prep & Interaction Diagrams]]
- [[_COMMUNITY_MD Job Setup & Diagnostics|MD Job Setup & Diagnostics]]
- [[_COMMUNITY_Trajectory Geometry Analysis|Trajectory Geometry Analysis]]
- [[_COMMUNITY_Step-02 Boltz-2 Production|Step-02 Boltz-2 Production]]
- [[_COMMUNITY_MSA & Prediction Dependencies|MSA & Prediction Dependencies]]
- [[_COMMUNITY_Step-06 Physics Validation|Step-06 Physics Validation]]
- [[_COMMUNITY_Shared Utilities (00_02)|Shared Utilities (00_02)]]
- [[_COMMUNITY_SN2 Geometry Primitives|SN2 Geometry Primitives]]
- [[_COMMUNITY_Mechanistic Scoring Functions|Mechanistic Scoring Functions]]
- [[_COMMUNITY_Prep Phase & Pose Drift|Prep Phase & Pose Drift]]
- [[_COMMUNITY_WaterMap & QSite Inputs|WaterMap & QSite Inputs]]
- [[_COMMUNITY_Conda Environment Packages|Conda Environment Packages]]
- [[_COMMUNITY_Catalytic Triad & QMMM Literature|Catalytic Triad & QM/MM Literature]]
- [[_COMMUNITY_Feasibility Gates & SN2 Sterics|Feasibility Gates & SN2 Sterics]]
- [[_COMMUNITY_CIF Collection & Residue Guards|CIF Collection & Residue Guards]]
- [[_COMMUNITY_Figure Constants & Tooling|Figure Constants & Tooling]]
- [[_COMMUNITY_MM-GBSA Statistics|MM-GBSA Statistics]]
- [[_COMMUNITY_SID Event Analysis|SID Event Analysis]]
- [[_COMMUNITY_Step-01 Sequence Merge|Step-01 Sequence Merge]]
- [[_COMMUNITY_Step-04 Dendrogram|Step-04 Dendrogram]]
- [[_COMMUNITY_Tier Gates & Project Identity|Tier Gates & Project Identity]]
- [[_COMMUNITY_Pipeline Orchestrator (00_00)|Pipeline Orchestrator (00_00)]]
- [[_COMMUNITY_Pipeline Architecture Contracts|Pipeline Architecture Contracts]]
- [[_COMMUNITY_Dataset & Catalytic Figure Folders|Dataset & Catalytic Figure Folders]]
- [[_COMMUNITY_Control Labelling & Combined Figures|Control Labelling & Combined Figures]]
- [[_COMMUNITY_Extraction Progress & Sharding|Extraction Progress & Sharding]]
- [[_COMMUNITY_Per-Rank MM-GBSA Figures|Per-Rank MM-GBSA Figures]]
- [[_COMMUNITY_Geometry & PBC Math|Geometry & PBC Math]]
- [[_COMMUNITY_Auxiliary Group 39|Auxiliary Group 39]]
- [[_COMMUNITY_Auxiliary Group 40|Auxiliary Group 40]]
- [[_COMMUNITY_Auxiliary Group 41|Auxiliary Group 41]]
- [[_COMMUNITY_Auxiliary Group 42|Auxiliary Group 42]]
- [[_COMMUNITY_Auxiliary Group 43|Auxiliary Group 43]]
- [[_COMMUNITY_Auxiliary Group 44|Auxiliary Group 44]]
- [[_COMMUNITY_Auxiliary Group 45|Auxiliary Group 45]]
- [[_COMMUNITY_Auxiliary Group 46|Auxiliary Group 46]]
- [[_COMMUNITY_Auxiliary Group 47|Auxiliary Group 47]]
- [[_COMMUNITY_Auxiliary Group 48|Auxiliary Group 48]]
- [[_COMMUNITY_Auxiliary Group 49|Auxiliary Group 49]]
- [[_COMMUNITY_Auxiliary Group 50|Auxiliary Group 50]]
- [[_COMMUNITY_Auxiliary Group 51|Auxiliary Group 51]]
- [[_COMMUNITY_Auxiliary Group 52|Auxiliary Group 52]]
- [[_COMMUNITY_Auxiliary Group 53|Auxiliary Group 53]]
- [[_COMMUNITY_Auxiliary Group 54|Auxiliary Group 54]]
- [[_COMMUNITY_Auxiliary Group 55|Auxiliary Group 55]]
- [[_COMMUNITY_Auxiliary Group 56|Auxiliary Group 56]]
- [[_COMMUNITY_Auxiliary Group 57|Auxiliary Group 57]]
- [[_COMMUNITY_Auxiliary Group 58|Auxiliary Group 58]]
- [[_COMMUNITY_Auxiliary Group 59|Auxiliary Group 59]]
- [[_COMMUNITY_Auxiliary Group 60|Auxiliary Group 60]]
- [[_COMMUNITY_Auxiliary Group 61|Auxiliary Group 61]]
- [[_COMMUNITY_Auxiliary Group 62|Auxiliary Group 62]]
- [[_COMMUNITY_Auxiliary Group 63|Auxiliary Group 63]]
- [[_COMMUNITY_Auxiliary Group 64|Auxiliary Group 64]]
- [[_COMMUNITY_Auxiliary Group 65|Auxiliary Group 65]]
- [[_COMMUNITY_Auxiliary Group 66|Auxiliary Group 66]]
- [[_COMMUNITY_Auxiliary Group 67|Auxiliary Group 67]]
- [[_COMMUNITY_Auxiliary Group 68|Auxiliary Group 68]]
- [[_COMMUNITY_Auxiliary Group 69|Auxiliary Group 69]]
- [[_COMMUNITY_Auxiliary Group 70|Auxiliary Group 70]]
- [[_COMMUNITY_Auxiliary Group 71|Auxiliary Group 71]]
- [[_COMMUNITY_Auxiliary Group 72|Auxiliary Group 72]]
- [[_COMMUNITY_Auxiliary Group 73|Auxiliary Group 73]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 80 edges
2. `main()` - 34 edges
3. `main()` - 34 edges
4. `_generate_comprehensive_figures_impl()` - 32 edges
5. `_echo()` - 32 edges
6. `process_single_job()` - 30 edges
7. `prep_and_convert_phase()` - 25 edges
8. `console_info()` - 24 edges
9. `run_mmgbsa_phase()` - 22 edges
10. `process_single_job()` - 21 edges

## Surprising Connections (you probably didn't know these)
- `Step 07 - MD/NAC analysis and QSite QM/MM defluorination` --references--> `mdanalysis 2.9.0`  [INFERRED]
  README.md → PFAS.yml
- `calculate_sn2_metrics()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `check_catalytic_geometry()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `generate_scientific_ranking_csv()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `main()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Multi-axis tier gate: effective mechanistic score, raw angle, competence floor, constellation cap, BDE ceiling** — readme_mechanistic_score_effective, readme_raw_angle_gate, readme_competence_floor, readme_criterion_b, readme_elite_bde_ceiling, readme_tier_ladder [EXTRACTED 1.00]
- **Control calibration anchor: FA/DFA/TFA against 3R3U and DEHA4 fix the reference tier ordering** — readme_fluoroacetate_fa, readme_difluoroacetate_dfa, readme_trifluoroacetate_tfa, readme_3r3u, readme_deha4, readme_control_set_6 [EXTRACTED 1.00]
- **Phase 3 physics chain: ESP merge, WaterMap, System Builder, Desmond MD, SID, MM-GBSA, QSite QM/MM** — readme_esp_charges, readme_watermap, readme_system_builder, readme_desmond_md, readme_sid, readme_prime_mmgbsa, readme_qsite_qmmm [EXTRACTED 1.00]

## Communities (74 total, 4 thin omitted)

### Community 0 - "Ranking & Pareto Figures"
Cohesion: 0.07
Nodes (49): analyse_conflicts(), calculate_pareto_fronts(), _fig23_multitarget(), _fig24_sankey(), _fig25_pfas_size(), _fig_folder07_pfas(), _is_control_mask(), _md_ready_df() (+41 more)

### Community 1 - "Step-06 Orchestration & Logging"
Cohesion: 0.10
Nodes (37): _complex_label(), _echo(), _emit_timings(), _fail(), find_esp(), _fmt_dur(), _hydrate(), _log() (+29 more)

### Community 2 - "Tier Quality Figures"
Cohesion: 0.10
Nodes (34): _fig23b_toptier_breakdown(), _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic(), _fig_14b_tt_interactions(), _fig_folder04_ai_confidence(), _fig_folder06_ligand(), _fig_path(), _kruskal_by_tier() (+26 more)

### Community 3 - "Merge & Alignment Cache"
Cohesion: 0.10
Nodes (35): atomic_to_csv(), colabfold_a3m_path(), colabfold_meta_path(), console_info(), console_separator(), extract_3r3u_sequence(), generate_scientific_ranking_csv(), load_cached_alignments() (+27 more)

### Community 4 - "Figure Panel Helpers"
Cohesion: 0.07
Nodes (35): _ext_match_03_style(), Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present., Compact, publication-style p-value formatting., Bring an extended-analysis panel onto 03's typography before it is written., Persist a figure as a PNG at the pipeline's publication resolution, then free it, No-op: panel letters are not drawn.      No other figure in this set carries an, A sample size that fits the ~0.5 in a tier occupies on a two-panel figure. (+27 more)

### Community 5 - "MD Heartbeat & Progress"
Cohesion: 0.08
Nodes (17): _close_bar(), Heartbeat, MDHeartbeat, _progress_line(), Refresh the single shared \r progress line (terminal only, never the log)., Close the open \r progress line with a newline (phase change / completion)., Periodically report a long Schrödinger step's progress by tailing its log., Report progress as a full logged line (runs on the CPU worker, concurrent with t (+9 more)

### Community 6 - "Figure Style & Console SSOT"
Cohesion: 0.10
Nodes (34): apply_figure_style(), Remove all ANSI/VT100 escape sequences from a string., The pipeline's one typography and canvas definition, applied to matplotlib's rcP, _strip_ansi(), console_info(), console_title(), format_job_label(), generate_comparative_residue_engagement() (+26 more)

### Community 7 - "Validation Figure Suite (03)"
Cohesion: 0.12
Nodes (31): _auto_label_colour(), _jf_box(), _jf_circos(), _jf_fit_fs(), _jf_importance(), _jf_ligshort(), _jf_manhattan(), _jf_metric_matrix() (+23 more)

### Community 8 - "CFG Central Configuration"
Cohesion: 0.09
Nodes (22): CFG, ===============================================================================, Central Configuration Repository - FAcDs Pipeline.      Attribute prefix → secti, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, _avail_ram_gb(), _await_mmgbsa_jobserver(), _cleanup_mmgbsa_shards(), Remove the per-job _MMGBSA_Shards scratch once the merged CSV is settled.      T (+14 more)

### Community 9 - "Report Manager & Data Loading"
Cohesion: 0.11
Nodes (31): latest_by_mtime(), Simultaneous console + file logger shared by the pipeline steps.      The log pa, The newest existing file by modification time, or None.      Selection is by mti, ReportManager, _aux_dir(), console_info(), generate_comprehensive_figures(), generate_ramachandran_figures() (+23 more)

### Community 10 - "Step-07 QM/MM Engine"
Cohesion: 0.08
Nodes (30): calculate_min_distance(), console_qmm_ready(), console_separator(), _darken(), _engage_zones(), extract_hybrid_smart_system(), _fel_grid(), _fel_surface() (+22 more)

### Community 11 - "Step-05 Prep & Interaction Diagrams"
Cohesion: 0.11
Nodes (28): _append_auxiliary_log(), _draw_interaction_diagram(), _im_don_acc(), _im_parse_pdb(), _im_project(), _im_render_diagram(), _im_ring(), _im_separate_atoms() (+20 more)

### Community 12 - "MD Job Setup & Diagnostics"
Cohesion: 0.07
Nodes (30): _csv_nonempty(), _diagnose_mmgbsa_failure(), _diagnose_shard_failure(), discover_handover(), export_watermap_csv(), _free_swap_gb(), _load_module(), _md_msj() (+22 more)

### Community 13 - "Trajectory Geometry Analysis"
Cohesion: 0.08
Nodes (27): _blockade_vec(), compute_wm_csv_stats(), _eaf_at(), extract_8residue_indices(), identity_from_cms(), kabsch_transform(), load_eaf_scalar_series(), mapped_resnum() (+19 more)

### Community 14 - "Step-02 Boltz-2 Production"
Cohesion: 0.10
Nodes (24): check_job_status(), cpu_usage_summary(), extract_short_fasta_id(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), heal_smiles_in_files() (+16 more)

### Community 15 - "MSA & Prediction Dependencies"
Cohesion: 0.09
Nodes (26): boltz 2.2.1 (pip), colabfold 1.6.1 (pip), mmseqs2 / hmmer 3.4 / mafft 7.525 / trimal 1.5.1, umap-learn 0.5.11 / scikit-learn 1.6.1 / statsmodels 0.14.5 / scikit-posthocs 0.14.0, Bostock et al. (2011) IEEE TVCG - D3 data-driven documents, Henikoff & Henikoff (1992) PNAS - BLOSUM substitution matrices, McInnes et al. (2018) JOSS - UMAP, Mirdita et al. (2022) Nature Methods - ColabFold (+18 more)

### Community 16 - "Step-06 Physics Validation"
Cohesion: 0.11
Nodes (24): _build_msj(), _clean_job_env(), _defl_ligand_atoms(), _defl_mapped_residues(), _defl_min_image(), _defl_resid(), _lig_atom_indices(), _mmgbsa_csv() (+16 more)

### Community 17 - "Shared Utilities (00_02)"
Cohesion: 0.08
Nodes (23): atomic_write_csv(), auto_label_colour(), calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner() (+15 more)

### Community 18 - "SN2 Geometry Primitives"
Cohesion: 0.09
Nodes (24): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+16 more)

### Community 19 - "Mechanistic Scoring Functions"
Cohesion: 0.11
Nodes (18): Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge() (+10 more)

### Community 20 - "Prep Phase & Pose Drift"
Cohesion: 0.12
Nodes (21): console_separator(), _is_reference_control(), load_rank_map(), main(), plot_machinery_distribution(), plot_pose_drift(), prep_and_convert_phase(), _print_run_delta_table() (+13 more)

### Community 21 - "WaterMap & QSite Inputs"
Cohesion: 0.10
Nodes (21): check_md_equilibration(), find_eaf_file(), generate_qsite_inputs(), _load_module(), load_triad_mapping(), load_watermap_csv(), load_watermap_reference_ca(), load_watermap_sites() (+13 more)

### Community 22 - "Conda Environment Packages"
Cohesion: 0.11
Nodes (21): biopython 1.84, PFAS conda environment specification, gemmi 0.6.5, mdanalysis 2.9.0, plip 3.0.0, matplotlib 3.10.8 / seaborn 0.13.2 / pandas 2.3.3 / numpy 1.26.4 / scipy 1.13.1, pymol-open-source 3.1.0, python 3.10.19 (+13 more)

### Community 23 - "Catalytic Triad & QM/MM Literature"
Cohesion: 0.10
Nodes (21): Asp134 acid catalyst (Asp-His dyad), Catalytic triad Asp110 / Asp134 / His277, Becke (1993) J Chem Phys - exact exchange (B3LYP), Holmquist (2000) Curr Protein Pept Sci - alpha/beta-hydrolase fold enzymes, Lee, Yang & Parr (1988) Phys Rev B - LYP correlation functional, Murphy et al. (2000) J Comput Chem - QSite QM/MM implementation, Rosta, Klahn & Warshel (2006) J Phys Chem B - ab initio QM/MM free-energy profiles, Yue et al. (2021) Environ Sci Technol - FAcD QM/MM degradation of fluorocarboxylic acids (+13 more)

### Community 24 - "Feasibility Gates & SN2 Sterics"
Cohesion: 0.10
Nodes (21): sn2_backside_occlusion (SN2 dead-end indicator B; > 2.0 A blocked), Beta-withdrawal counter crossing a single ether oxygen, Bidentate carboxylate clamp Arg111 / Arg114, Bento & Bickelhaupt (2008) J Org Chem - backside SN2 sterics, Bondi (1964) J Phys Chem - van der Waals volumes and radii, O'Hagan (2008) Chem Soc Rev - understanding organofluorine chemistry / C-F strength, Substrate-feasibility floor TIER_COMP_MIN (0.50 / 0.40 at Tier_2A/2B), competence_score (gated continuous catalytic competence) (+13 more)

### Community 25 - "CIF Collection & Residue Guards"
Cohesion: 0.10
Nodes (20): _check_residue_identity_guard(), collect_best_cifs(), extract_chain_l_mol(), index_existing_files(), load_catalytic_anchor_map(), load_machinery_map(), load_md_selected_jobs(), load_metadata() (+12 more)

### Community 26 - "Figure Constants & Tooling"
Cohesion: 0.17
Nodes (13): console_info(), _fig_constants(), _fig_log(), load_reference_data(), Scans the 1_Input_Data folder to map names to Sequences and SMILES.     Returns:, Return a dict of figure-generation constants from CFG., Unified logging for figure generation - routes through console_info., Handles verification and automated installation of visual tools. (+5 more)

### Community 27 - "MM-GBSA Statistics"
Cohesion: 0.12
Nodes (20): _avg_dg(), _block_bootstrap_median_ci(), _boltzmann_mean_dg(), _cliffs_delta(), _delta_word(), _dg_failures(), _draw_time_cumulative(), _effective_n() (+12 more)

### Community 28 - "SID Event Analysis"
Cohesion: 0.13
Nodes (20): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Classify every desmond_md_job_R_* directory and return the subset     that still, Step 1: event_analysis.py - generates the SID-in.eaf descriptor., Step 2: analyze_simulation.py - generates the SID-out.eaf result vector.      Us (+12 more)

### Community 29 - "Step-01 Sequence Merge"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 30 - "Step-04 Dendrogram"
Cohesion: 0.18
Nodes (18): clean_id(), console_info(), console_separator(), generate_phylogenies(), generate_upgma_newick(), get_kmer_counts(), _load_module(), main() (+10 more)

### Community 31 - "Tier Gates & Project Identity"
Cohesion: 0.12
Nodes (18): PDB 3R3U crystal reference (R. palustris FAcD, 1.60 A), Shaban Ahmad (developer and project lead, University of Copenhagen), Tue K. Nielsen (scientific supervisor, University of Copenhagen), Size-fair backbone-clash veto (fraction >= 0.15 and count >= 3), Chan et al. (2011) JACS - mapping the reaction coordinates of enzymatic defluorination, Confidence demotion (Tier_1A with Boltz confidence < 0.85 drops to 1B), Coupled elite gate (MECH_ELITE_HI 0.90 OR 0.85 + constellation 0.74), Criterion A - active-site integrity (eight catalytic residues present) (+10 more)

### Community 32 - "Pipeline Orchestrator (00_00)"
Cohesion: 0.22
Nodes (13): _cancel_schrodinger_jobs(), _ensure_jobserver(), _print_step_menu(), _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_FAcDs.sh script, _sync_staging_log() (+5 more)

### Community 33 - "Pipeline Architecture Contracts"
Cohesion: 0.16
Nodes (17): Jacobson et al. (2004) Proteins - Schrodinger Prime refinement, Li et al. (2011) Proteins - VSGB 2.0 implicit solvent for MM-GBSA, Michaud-Agrawal et al. (2011) J Comput Chem - MDAnalysis, Pipelined GPU MD to CPU post-processing per rank, R_<Scientific_Rank> job-naming contract across Steps 05-07, LazyTrajectory chunk-streamed frame analysis (5000-frame chunks), systemd-oomd masking with crash-proof sentinel restore, 00_00_run_pipeline_FAcDs.sh (bash orchestrator) (+9 more)

### Community 34 - "Dataset & Catalytic Figure Folders"
Cohesion: 0.19
Nodes (16): _control_star_df(), _ctrl_star(), _diag10_model_agreement(), _fig_folder03_dataset(), _fig_folder05_catalytic(), generate_additional_figures(), _md_ready_stars_cat(), Folder 03_Dataset_and_Alignment_Overview - coverage, tiers, alignment grades. (+8 more)

### Community 35 - "Control Labelling & Combined Figures"
Cohesion: 0.13
Nodes (16): _analysis_dir(), _ctrl_label(), _ctrl_palette(), make_md_qc_figure(), _natural_rank(), plot_defluor_combined(), The ligand's short name (FA / DFA / TFA) from CFG - the same abbreviations every, 3R3U-FA' for a control rank, else the short PFAS name (fallback R_<rk>). (+8 more)

### Community 36 - "Extraction Progress & Sharding"
Cohesion: 0.14
Nodes (11): _eta_str(), _extract_with_progress(), Emit a CPU-worker progress line (extraction / SID / MM-GBSA) as a full logged li, Extract members from tgz into wd (dropping one leading path component), drawing, Put the production trajectory at the job-dir root as {job}_trj + {job}.ene and r, A ' · ETA <dur>' suffix from a linear extrapolation of the current rate - empty, One in-place progress line for a sharded MM-GBSA run.      Aggregates across the, (frames read across all shards, shards currently in Prime minimisation). (+3 more)

### Community 37 - "Per-Rank MM-GBSA Figures"
Cohesion: 0.14
Nodes (16): _failure_windows(), _lookup_controls(), _lookup_ligands(), _lookup_tiers(), _mmgbsa_dg_series(), _plot_rank_mmgbsa(), Extract the per-frame ΔG_bind series, tolerant of column-name variants., Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier     co (+8 more)

### Community 38 - "Geometry & PBC Math"
Cohesion: 0.17
Nodes (15): _box_diag(), calculate_improper_dihedral(), calculate_min_distance(), distance(), _ensure_box_3x3(), get_mic_vector(), mic_dists_2d(), ndarray (+7 more)

### Community 39 - "Auxiliary Group 39"
Cohesion: 0.16
Nodes (14): clean_spines(), Academic-style axes: remove top/right spines, thin the remaining borders.      C, _active_site_series(), _draw_trajectory_figures(), generate_active_site_dynamics(), generate_individual_dashboard(), generate_mmgbsa_trace(), 2-panel per-job dashboard: SN2 scatter and dual-trace anchoring time series. (+6 more)

### Community 40 - "Auxiliary Group 40"
Cohesion: 0.19
Nodes (14): compute_ligand_properties(), _draw_rama_background(), Any, Path, _rama_classify(), _rama_stats(), Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, Classify a phi/psi pair as Favored, Allowed, or Outlier. (+6 more)

### Community 41 - "Auxiliary Group 41"
Cohesion: 0.16
Nodes (5): _ConsoleRuleFilter, install_console_rule_filter(), A stdout wrapper that collapses consecutive separator rules.      The logs grow, Collapse stacked separator rules for the rest of this process's output., Record a line in the log file WITHOUT printing it to the console.          For o

### Community 42 - "Auxiliary Group 42"
Cohesion: 0.18
Nodes (14): _generate_comprehensive_figures_impl(), _mm_load_variance_df(), _mm_nac_lines(), _mm_variance_by_model(), _mm_variance_conf90(), _mm_variance_cut(), _mm_variance_iptm90(), Per-model variance CSV (one row per model per complex) joined to the pipeline ti (+6 more)

### Community 43 - "Auxiliary Group 43"
Cohesion: 0.15
Nodes (14): _collect_qsite_results(), _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), parse_qsite_barrier(), parse_qsite_profile(), plot_qsite_reaction_profile(), Per-scan-point energy series, returned in KCAL/MOL, from a Jaguar/QSite     rela (+6 more)

### Community 44 - "Auxiliary Group 44"
Cohesion: 0.20
Nodes (12): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi, Classify a ligand atom as 'anion'- or 'cation'-capable for salt-bridge     detec (+4 more)

### Community 45 - "Auxiliary Group 45"
Cohesion: 0.18
Nodes (8): _cancel_launched_jobs(), _install_job_cleanup(), _kill_launched_procs(), OomdGuard, Cancel every still-registered job on the job server. Idempotent (guarded), safe, Terminate every still-registered LOCAL subprocess (SID, MM-GBSA drivers) and its, Cancel the run's job-server jobs AND kill its local subprocesses on normal exit, Optionally mask systemd-oomd to prevent Out-Of-Memory kills during long SID runs

### Community 46 - "Auxiliary Group 46"
Cohesion: 0.17
Nodes (12): Farajollahi et al. (2024) ACS Omega - D4B dehalogenase defluorination, Khusnutdinova et al. (2023) FEBS J - difluoroacetate is a genuine FAcD substrate, Wackett (2022) Microb Biotechnol - FAcD substrate scope and TFA recalcitrance, Six control cases (3R3U and DEHA4 x FA/DFA/TFA), DEHA4 (Delftia acidovorans D4B) control sequence, Difluoroacetate (DFA) - positive control, bent ~108 deg, Tier_2B, Ether / next-generation PFAS (GenX, ADONA, C6O4), Fluoroacetate (FA) - positive control, reference Tier_2A (+4 more)

### Community 47 - "Auxiliary Group 47"
Cohesion: 0.18
Nodes (11): graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit(), Constructs a three-dimensional RDKit molecule from a SMILES string.     This is, Extracts raw XYZ coordination matrices from RDKit molecule data blocks., Determines the optimal rotation/translation matrix required to superimpose spati, Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not b, FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perce (+3 more)

### Community 48 - "Auxiliary Group 48"
Cohesion: 0.18
Nodes (11): Binding-probability logit centring and gain calibration, Cation_Capped frame flag (CATION_CAP_DIST 3.0 A), chem_verified flag bars unresolved scissile centres from Tier_1A, Lu et al. (2021) JCTC - OPLS4 force field, Ion/salt exclusion within 5 A of the ligand, Missing is not zero doctrine, Rationale - a Na+ on the Asp-Odelta screens the nucleophile and corrupts NAC, Rationale - TIP3P/SPC chosen because OPLS4, Prime and WaterMap are validated against them (+3 more)

### Community 49 - "Auxiliary Group 49"
Cohesion: 0.20
Nodes (10): Distance (Å) -> 0-1 proximity; NaN (sentinel / not engaged) -> 0., Occupancy adequacy tent: too empty and overflow both score low; a substrate-size, Per-complex sub-scores + geometry_readiness, feasibility (=02 competence_score), Pocket + 8-residue reaction-readiness per ligand, ordered by molecular size (4 p, _rr_family(), _rr_occ(), _rr_prox(), _rr_scores() (+2 more)

### Community 50 - "Auxiliary Group 50"
Cohesion: 0.20
Nodes (10): _build_hd1_line(), enforce_catalytic_protonation(), measure_sn2_geometry(), preparation_step(), Measure the SN2 geometry of a structure - CIF or PDB - with the pipeline's own j, Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev, Task 2: Protein Preparation.     Logic: Checks consistency between Raw and Prep (+2 more)

### Community 51 - "Auxiliary Group 51"
Cohesion: 0.20
Nodes (10): check_prep_needed(), cif_to_pdb_gemmi(), generate_raw_step(), get_source_tag(), Reads the REMARK 999 SOURCE_CIF tag from a PDB file., Injects metadata into the PDB Header.     Crucially, adds 'SOURCE_CIF' to track, Converts mmCIF to PDB using Gemmi.     - Moves Ligands (Non-Standard Residues) t, Task 1: Generate Raw PDB from CIF.     CIF is passed directly from 2_Best_Comple (+2 more)

### Community 52 - "Auxiliary Group 52"
Cohesion: 0.22
Nodes (9): console_info(), console_separator(), console_title(), Logger, Canonical file logger for all pipeline scripts.      Creates a file-only handler, Bold section header - printed and optionally written to log file., Two-space-indented info line - printed and optionally written to log file., Horizontal rule - printed and optionally written to log file.      Parameters (+1 more)

### Community 53 - "Auxiliary Group 53"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 54 - "Auxiliary Group 54"
Cohesion: 0.22
Nodes (9): Bowers et al. (2006) SC'06 - Desmond scalable MD algorithms, Olsson et al. (2011) JCTC - PROPKA3, Sastry et al. (2013) J Comput Aided Mol Des - protein and ligand preparation, Desmond MD production (NPT, 300 K, 1000 ns, 100000 frames), Positional restraint (ligand 5.0, backbone 2.0 kcal/mol/A^2), Pose drift CIF to prepared (6.3 deg mean, max 18.5 deg; +0.30 A), Schrodinger PrepWizard preparation (PropKa pH 8.0, 0.15 A restrained min), Rationale - prepared-pose drift is recorded and flagged, never gated on (+1 more)

### Community 55 - "Auxiliary Group 55"
Cohesion: 0.22
Nodes (9): Bruice (2002) Acc Chem Res - efficiency of enzymatic catalysis, Hur & Bruice (2003) PNAS - near attack conformation approach, Lightstone & Bruice (1996) JACS - NAC ground-state conformations, Is_Defluorinating verdict (QSITE_F_CHARGE_CLEAVED <= -0.5 e), Near Attack Conformation (NAC) criteria, NAC dwell time and occupancy over the trajectory, Nucleophile-C distance (attacking Asp-Odelta to C-alpha), Rationale - one oxygen must satisfy both halves of the NAC (no chimeric nucleophile) (+1 more)

### Community 56 - "Auxiliary Group 56"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 57 - "Auxiliary Group 57"
Cohesion: 0.25
Nodes (8): Abel et al. (2008) JACS - active-site solvent thermodynamics (WaterMap), Explicit-solvent droplet MM region (15 A, frozen surface, no implicit solvation), Fold_RMSD_A flag (MD_FOLD_RMSD_MAX 3.0 A withholds WaterMap term), Kabsch C-alpha superposition of frames onto the WaterMap reference, Rationale - the C-alpha centroid moves 30-45 A over 1 us, so sites must be carried across, Rationale - no implicit solvation keyword, the droplet would double-count, Rationale - ligand retained so the site map describes the same holo system as MD, WaterMap GCMC hydration-site mapping (holo, ligand retained)

### Community 58 - "Auxiliary Group 58"
Cohesion: 0.29
Nodes (7): Boltz-2 diffusion co-folding model, Boltz-2 confidence metrics (ipTM, pLDDT, cross-PAE, pTM), Passaro et al. (2025) bioRxiv - Boltz-2, Wohlwend et al. (2024) bioRxiv - Boltz-1, model_degrader_consensus (reproducibility tie-breaker), Rationale - confidence metrics assess plausibility, not catalytic competence, Representative-pose selection across diffusion samples

### Community 59 - "Auxiliary Group 59"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 60 - "Auxiliary Group 60"
Cohesion: 0.33
Nodes (6): analysis_worker_loop(), append_rows_to_csv(), flatten_job_result(), Flattens a process_single_job result dict into a single-level row for CSV writin, Appends a list of result rows to the master CSV, deduplicates, sorts by job_name, Background worker that analyses freshly GPU-predicted jobs and appends rows to t

### Community 61 - "Auxiliary Group 61"
Cohesion: 0.33
Nodes (6): fetch_msa_direct(), Generates a stable hash for FASTA sequences to verify consistency during resumed, Validates the A3M MSA file to detect corruption. Returns True if the file struct, Fetches the MSA directly from the ColabFold REST API. This circumvents the Boltz, sequence_hash(), validate_a3m_file()

### Community 62 - "Auxiliary Group 62"
Cohesion: 0.33
Nodes (6): _esp_run(), generate_esp_charges(), plot_esp_alpha_carbon(), Run a snippet under $SCHRODINGER/run - its interpreter, not ours., QM (Jaguar ESP) partial charges for every prepared MD ligand. Returns the summar, The charge on the carbon the nucleophile actually attacks - the number OPLS4 can

### Community 63 - "Auxiliary Group 63"
Cohesion: 0.33
Nodes (6): Asp110 nucleophile (backside attack on C-alpha), Jansen, van Beers & Mayer (2026) Angew Chem - engineering FAcD on non-natural organofluorides, Jitsumori et al. (2009) J Bacteriol - FAcD from Burkholderia sp. FA1, Verschueren et al. (1993) Nature - haloalkane dehalogenase mechanism, Fluoroacetate dehalogenase (FAcD), SN2 Walden-inversion defluorination mechanism

### Community 65 - "Auxiliary Group 65"
Cohesion: 0.50
Nodes (4): calculate_angle(), calculate_burgi_dunitz(), Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit, Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for

### Community 66 - "Auxiliary Group 66"
Cohesion: 0.50
Nodes (4): analyse_candidate_structure(), load_structure_safe(), Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilit, Superimposes the candidate structure onto the DeHa4 control model.     Computes

### Community 67 - "Auxiliary Group 67"
Cohesion: 0.50
Nodes (4): analyse_model_task(), Global worker routine functionally decoupled for robust serialisation capabiliti, Selects the representative pose across the Boltz diffusion samples tier-first: t, select_best_degrader_model()

### Community 68 - "Auxiliary Group 68"
Cohesion: 0.50
Nodes (4): compute_cross_interface_pae(), load_extra_boltz_metrics(), Calculates the Predicted Aligned Error (PAE), focusing specifically on the prote, Parses background statistical validation tensors directly from Boltz NPZ file du

### Community 69 - "Auxiliary Group 69"
Cohesion: 0.50
Nodes (4): init_gpu_reservations(), query_gpu_memory(), A resilient system query that safely falls back to CPU memory context if nvidia-, Maps unallocated VRAM memory segments natively upon system initialisation.

### Community 70 - "Auxiliary Group 70"
Cohesion: 0.50
Nodes (4): map_active_site_residues(), Append a newly computed alignment to the on-disk cache file, under an exclusive, Performs a global sequence alignment between the canonical reference sequence, update_alignment_stats()

## Knowledge Gaps
- **89 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `Phase 1 - high-throughput screening (steps 01-04)` (+84 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `CFG Central Configuration` to `Ranking & Pareto Figures`, `Step-06 Orchestration & Logging`, `Tier Quality Figures`, `Merge & Alignment Cache`, `Figure Panel Helpers`, `Figure Style & Console SSOT`, `Validation Figure Suite (03)`, `Report Manager & Data Loading`, `Step-07 QM/MM Engine`, `MD Job Setup & Diagnostics`, `Trajectory Geometry Analysis`, `Step-06 Physics Validation`, `Mechanistic Scoring Functions`, `Prep Phase & Pose Drift`, `WaterMap & QSite Inputs`, `CIF Collection & Residue Guards`, `MM-GBSA Statistics`, `SID Event Analysis`, `Step-04 Dendrogram`, `Control Labelling & Combined Figures`, `Per-Rank MM-GBSA Figures`, `Auxiliary Group 39`, `Auxiliary Group 42`, `Auxiliary Group 43`, `Auxiliary Group 62`, `Auxiliary Group 70`?**
  _High betweenness centrality (0.467) - this node is a cross-community bridge._
- **Why does `_generate_comprehensive_figures_impl()` connect `Auxiliary Group 42` to `Ranking & Pareto Figures`, `Dataset & Catalytic Figure Folders`, `Tier Quality Figures`, `Figure Panel Helpers`, `Validation Figure Suite (03)`, `CFG Central Configuration`, `Report Manager & Data Loading`, `Auxiliary Group 49`?**
  _High betweenness centrality (0.074) - this node is a cross-community bridge._
- **Why does `main()` connect `Merge & Alignment Cache` to `Auxiliary Group 66`, `Auxiliary Group 69`, `CFG Central Configuration`, `Step-02 Boltz-2 Production`, `Shared Utilities (00_02)`, `Auxiliary Group 60`, `Auxiliary Group 61`?**
  _High betweenness centrality (0.073) - this node is a cross-community bridge._
- **Are the 65 inferred relationships involving `CFG` (e.g. with `calculate_sn2_metrics()` and `check_catalytic_geometry()`) actually correct?**
  _`CFG` has 65 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `main()` (e.g. with `CFG` and `safe_name()`) actually correct?**
  _`main()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _538 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Ranking & Pareto Figures` be split into smaller, more focused modules?**
  _Cohesion score 0.07397959183673469 - nodes in this community are weakly interconnected._