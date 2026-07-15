# Graph Report - .  (2026-07-15)

## Corpus Check
- 10 files · ~99,999 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 826 nodes · 1701 edges · 51 communities (50 shown, 1 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 106 edges (avg confidence: 0.73)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_05 TopN and PDB Preparation|05 TopN and PDB Preparation]]
- [[_COMMUNITY_03 Validation Figures|03 Validation Figures]]
- [[_COMMUNITY_04 Dendrogram|04 Dendrogram]]
- [[_COMMUNITY_03 Validation Figures (2)|03 Validation Figures (2)]]
- [[_COMMUNITY_02 Production|02 Production]]
- [[_COMMUNITY_03 Validation Figures (3)|03 Validation Figures (3)]]
- [[_COMMUNITY_03 Validation Figures (4)|03 Validation Figures (4)]]
- [[_COMMUNITY_07 MD QMMM Defluorination|07 MD QMMM Defluorination]]
- [[_COMMUNITY_00_02 Project Utils|00_02 Project Utils]]
- [[_COMMUNITY_02 Production (2)|02 Production (2)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA|06 SID Prime-MMGBSA]]
- [[_COMMUNITY_00_01 Project Config|00_01 Project Config]]
- [[_COMMUNITY_02 Production (3)|02 Production (3)]]
- [[_COMMUNITY_03 Validation Figures (5)|03 Validation Figures (5)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (2)|05 TopN and PDB Preparation (2)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (3)|05 TopN and PDB Preparation (3)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (4)|05 TopN and PDB Preparation (4)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (2)|06 SID Prime-MMGBSA (2)]]
- [[_COMMUNITY_01 Merge|01 Merge]]
- [[_COMMUNITY_07 MD QMMM Defluorination (2)|07 MD QMMM Defluorination (2)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (3)|07 MD QMMM Defluorination (3)]]
- [[_COMMUNITY_02 Production (4)|02 Production (4)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (3)|06 SID Prime-MMGBSA (3)]]
- [[_COMMUNITY_03 Validation Figures (6)|03 Validation Figures (6)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (4)|07 MD QMMM Defluorination (4)]]
- [[_COMMUNITY_02 Production (5)|02 Production (5)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (4)|06 SID Prime-MMGBSA (4)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (5)|06 SID Prime-MMGBSA (5)]]
- [[_COMMUNITY_00_02 Project Utils (2)|00_02 Project Utils (2)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (5)|07 MD QMMM Defluorination (5)]]
- [[_COMMUNITY_02 Production (6)|02 Production (6)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (6)|06 SID Prime-MMGBSA (6)]]
- [[_COMMUNITY_00_02 Project Utils (3)|00_02 Project Utils (3)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (6)|07 MD QMMM Defluorination (6)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (7)|07 MD QMMM Defluorination (7)]]
- [[_COMMUNITY_00_02 Project Utils (4)|00_02 Project Utils (4)]]
- [[_COMMUNITY_00_03 Environment|00_03 Environment]]
- [[_COMMUNITY_05 TopN and PDB Preparation (5)|05 TopN and PDB Preparation (5)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (7)|06 SID Prime-MMGBSA (7)]]
- [[_COMMUNITY_06 SID Prime-MMGBSA (8)|06 SID Prime-MMGBSA (8)]]
- [[_COMMUNITY_00_02 Project Utils (5)|00_02 Project Utils (5)]]
- [[_COMMUNITY_00_02 Project Utils (6)|00_02 Project Utils (6)]]
- [[_COMMUNITY_03 Validation Figures (7)|03 Validation Figures (7)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (6)|05 TopN and PDB Preparation (6)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (7)|05 TopN and PDB Preparation (7)]]
- [[_COMMUNITY_07 MD QMMM Defluorination (8)|07 MD QMMM Defluorination (8)]]
- [[_COMMUNITY_00_02 Project Utils (7)|00_02 Project Utils (7)]]
- [[_COMMUNITY_02 Production (7)|02 Production (7)]]
- [[_COMMUNITY_02 Production (8)|02 Production (8)]]
- [[_COMMUNITY_03 Validation Figures (8)|03 Validation Figures (8)]]
- [[_COMMUNITY_05 TopN and PDB Preparation (8)|05 TopN and PDB Preparation (8)]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 65 edges
2. `main()` - 34 edges
3. `process_single_job()` - 30 edges
4. `_generate_comprehensive_figures_impl()` - 27 edges
5. `prep_and_convert_phase()` - 24 edges
6. `console_info()` - 22 edges
7. `process_single_job()` - 21 edges
8. `_echo()` - 19 edges
9. `run_mmgbsa_phase()` - 19 edges
10. `check_catalytic_geometry()` - 18 edges

## Surprising Connections (you probably didn't know these)
- `calculate_sn2_metrics()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `check_catalytic_geometry()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `generate_scientific_ranking_csv()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `main()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `map_active_site_residues()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_FAcDs.py → 00_01_Project_Config_FAcDs.py

## Import Cycles
- None detected.

## Communities (51 total, 1 thin omitted)

### Community 0 - "05 TopN and PDB Preparation"
Cohesion: 0.08
Nodes (39): check_prep_needed(), _check_residue_identity_guard(), collect_best_cifs(), console_info(), generate_esp_charges(), get_source_tag(), index_existing_files(), load_catalytic_anchor_map() (+31 more)

### Community 1 - "03 Validation Figures"
Cohesion: 0.10
Nodes (36): _kruskal_by_tier(), _median_ci95(), Register one test into the family. extra carries effect size, group ns, medians, Kruskal–Wallis across tiers with epsilon-squared effect size. Returns annotation, Paired Wilcoxon signed-rank test of ipTM vs pTM across complexes., Silhouette of tier labels in UMAP space with a label-permutation p-value.     Su, Bootstrap 95% confidence interval of the MEDIAN (percentile method).      Return, Nucleophile–substrate distance, resolving the several historical names.      The (+28 more)

### Community 2 - "04 Dendrogram"
Cohesion: 0.09
Nodes (25): _ConsoleRuleFilter, install_console_rule_filter(), Collapse stacked separator rules for the rest of this process's output., Simultaneous console + file logger shared by the pipeline steps.      The log pa, Record a line in the log file WITHOUT printing it to the console.          For o, A stdout wrapper that collapses consecutive separator rules.      The logs grow, ReportManager, clean_id() (+17 more)

### Community 3 - "03 Validation Figures (2)"
Cohesion: 0.09
Nodes (31): Nucleophile–substrate distance, resolving the several historical names.      The, Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present., Coerce a column to numeric, returning an empty series when it is absent., Compact, publication-style p-value formatting., No-op: panel letters are not drawn.      No other figure in this set carries an, A sample size that fits the ~0.5 in a tier occupies on a two-panel figure., The panel's test result, on its own line ABOVE the axes.      Inside the axes it (+23 more)

### Community 4 - "02 Production"
Cohesion: 0.11
Nodes (29): analysis_worker_loop(), append_rows_to_csv(), atomic_to_csv(), console_info(), console_separator(), flatten_job_result(), generate_scientific_ranking_csv(), load_cached_alignments() (+21 more)

### Community 5 - "03 Validation Figures (3)"
Cohesion: 0.12
Nodes (27): analyse_conflicts(), _aux_dir(), console_info(), generate_comprehensive_figures(), generate_ramachandran_figures(), load_and_prep_data(), _load_module(), main() (+19 more)

### Community 6 - "03 Validation Figures (4)"
Cohesion: 0.09
Nodes (26): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+18 more)

### Community 7 - "07 MD QMMM Defluorination"
Cohesion: 0.09
Nodes (25): _blockade_vec(), compute_wm_csv_stats(), _eaf_at(), identity_from_cms(), kabsch_transform(), load_eaf_scalar_series(), mapped_resnum(), _nac_dwell_stats() (+17 more)

### Community 8 - "00_02 Project Utils"
Cohesion: 0.08
Nodes (23): atomic_write_csv(), auto_label_colour(), calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), _jobserver_addresses(), print_script_banner() (+15 more)

### Community 9 - "02 Production (2)"
Cohesion: 0.10
Nodes (24): calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge(), graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit() (+16 more)

### Community 10 - "06 SID Prime-MMGBSA"
Cohesion: 0.11
Nodes (22): _block_bootstrap_median_ci(), _cliffs_delta(), _csv_nonempty(), _delta_word(), _diagnose_shard_failure(), _draw_time_cumulative(), _effective_n(), _failure_windows() (+14 more)

### Community 11 - "00_01 Project Config"
Cohesion: 0.10
Nodes (9): CFG, ===============================================================================, Single source of truth for the holistic mechanistic score (0–1) — the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, Central Configuration Repository — FAcDs Pipeline.      Attribute prefix → secti, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, _auto_label_colour() (+1 more)

### Community 12 - "02 Production (3)"
Cohesion: 0.11
Nodes (22): analyse_candidate_structure(), colabfold_a3m_path(), colabfold_meta_path(), extract_3r3u_sequence(), heal_smiles_in_files(), load_structure_safe(), map_active_site_residues(), perform_control_calibration() (+14 more)

### Community 13 - "03 Validation Figures (5)"
Cohesion: 0.13
Nodes (21): calculate_pareto_fronts(), _diag10_model_agreement(), _fig24_sankey(), _fig25_pfas_size(), generate_additional_figures(), generate_extended_figures(), _md_ready_df(), _md_ready_stars_cat() (+13 more)

### Community 14 - "05 TopN and PDB Preparation (2)"
Cohesion: 0.13
Nodes (19): console_separator(), _esp_run(), extract_chain_l_mol(), load_reference_data(), main(), _plip_angles_ok(), plot_machinery_distribution(), Run a snippet under $SCHRODINGER/run — its interpreter, not ours. (+11 more)

### Community 15 - "05 TopN and PDB Preparation (3)"
Cohesion: 0.12
Nodes (20): _append_auxiliary_log(), _draw_interaction_diagram(), _im_don_acc(), _im_parse_pdb(), _im_project(), _im_render_diagram(), _im_ring(), _im_separate_atoms() (+12 more)

### Community 16 - "05 TopN and PDB Preparation (4)"
Cohesion: 0.16
Nodes (13): _fig_constants(), _fig_log(), load_metadata(), Return a dict of figure-generation constants from CFG., Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir., Unified logging for figure generation — routes through console_info., Handles verification and automated installation of visual tools., Check if PLIP is installed as a Python module. (+5 more)

### Community 17 - "06 SID Prime-MMGBSA (2)"
Cohesion: 0.12
Nodes (20): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Count frames written into a SID-out.eaf (the analysed result vector).      Parse, Ground-truth frame count of a Desmond trajectory via the Schrödinger     traj AP, A SID-out.eaf is complete when its ``Result=[…]`` vector holds one full     toke (+12 more)

### Community 18 - "01 Merge"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 19 - "07 MD QMMM Defluorination (2)"
Cohesion: 0.13
Nodes (19): _darken(), _engage_zones(), generate_reactive_pose_figures(), _load_reactive_pose_data(), _master_container(), _mmgbsa_components(), plot_machinery_engagement(), plot_mmgbsa_decomposition() (+11 more)

### Community 20 - "07 MD QMMM Defluorination (3)"
Cohesion: 0.14
Nodes (18): apply_figure_style(), clean_spines(), The pipeline's one typography and canvas definition, applied to matplotlib's rcP, Academic-style axes: remove top/right spines, thin the remaining borders.      C, format_job_label(), generate_comparative_residue_engagement(), generate_defluorination_landscape(), generate_global_comparative_dashboard() (+10 more)

### Community 21 - "02 Production (4)"
Cohesion: 0.13
Nodes (16): cpu_usage_summary(), extract_short_fasta_id(), fetch_msa_direct(), init_gpu_reservations(), _load_module(), query_gpu_memory(), Load a Python file as a module using its filesystem path, regardless of filename, A resilient system query that safely falls back to CPU memory context if nvidia- (+8 more)

### Community 22 - "06 SID Prime-MMGBSA (3)"
Cohesion: 0.15
Nodes (13): print_elapsed(), Print total pipeline-script elapsed time bookended by heavy separators.      Par, cleanup_schrodinger_dirs(), _echo(), main(), OomdGuard, _open_step_log(), Open 6_Physics_Validation/00_SID_MMGBSA.log for this run (fresh each run). (+5 more)

### Community 23 - "03 Validation Figures (6)"
Cohesion: 0.27
Nodes (17): _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic(), _fig_14b_tt_interactions(), _fig_18b_tt_landscape(), _generate_comprehensive_figures_impl(), ndarray, Render the full publication figure suite (folders 02–07).      Builds the datase, Draw a boxed statistics annotation inside an axis (no effect if text is empty). (+9 more)

### Community 24 - "07 MD QMMM Defluorination (4)"
Cohesion: 0.16
Nodes (16): calculate_min_distance(), console_qmm_ready(), console_separator(), console_title(), extract_hybrid_smart_system(), load_triad_mapping(), main(), _print_labeled() (+8 more)

### Community 25 - "02 Production (5)"
Cohesion: 0.12
Nodes (16): check_job_status(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), mirror_best_cif(), process_single_job(), Helper tool designed to retrieve alignment data employing secure string-safe key (+8 more)

### Community 26 - "06 SID Prime-MMGBSA (4)"
Cohesion: 0.17
Nodes (7): Heartbeat, Periodically report a long Schrödinger step's progress by tailing its log., Rewrite the single progress line in place with a carriage return., Close the open in-place \r line with a newline (phase change / exit)., Largest last-match across all configured patterns (phase-tolerant)., Best-effort (done, total) for the Prime phase from a ``*-prime*.log`` in, Highest 'Structure N (of M)' seen in the main log or any ``*-prime*.log``

### Community 27 - "06 SID Prime-MMGBSA (5)"
Cohesion: 0.15
Nodes (14): _avail_ram_gb(), _diagnose_mmgbsa_failure(), _mmgbsa_csv(), Score the whole trajectory as concurrent contiguous frame shards.      Each shar, Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:     ret, Locate a thermal_mmgbsa results CSV in `job_dir` (name varies by version).     R, Best-effort human-readable cause when MM-GBSA exits non-zero, read from the, Free RAM the kernel expects to hand out without swapping (MemAvailable). (+6 more)

### Community 28 - "00_02 Project Utils (2)"
Cohesion: 0.17
Nodes (15): _box_diag(), calculate_improper_dihedral(), calculate_min_distance(), distance(), _ensure_box_3x3(), get_mic_vector(), mic_dists_2d(), ndarray (+7 more)

### Community 29 - "07 MD QMMM Defluorination (5)"
Cohesion: 0.14
Nodes (15): find_eaf_file(), _load_module(), load_watermap_csv(), load_watermap_reference_ca(), Path, _qsite_scan_failure_reason(), The Cα coordinates of the structure the WaterMap sites were computed IN, keyed b, Locate protein-ligand EAF file with priority:       P1: *-out*pl*.eaf  (Rank_1: (+7 more)

### Community 30 - "02 Production (6)"
Cohesion: 0.16
Nodes (14): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), load_atoms_from_structure(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi (+6 more)

### Community 31 - "06 SID Prime-MMGBSA (6)"
Cohesion: 0.20
Nodes (14): _load_module(), _lookup_ligands(), _lookup_tiers(), _mmgbsa_dg_series(), _natural_rank(), Path, Extract the per-frame ΔG_bind series, tolerant of column-name variants., Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier     co (+6 more)

### Community 32 - "00_02 Project Utils (3)"
Cohesion: 0.24
Nodes (11): _draw_rama_background(), Any, _rama_classify(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favoured, allowed, and outlier regions., Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes., Save a side-by-side comparison Ramachandran PNG. (+3 more)

### Community 33 - "07 MD QMMM Defluorination (6)"
Cohesion: 0.18
Nodes (11): Remove all ANSI/VT100 escape sequences from a string., _strip_ansi(), check_md_equilibration(), console_info(), extract_8residue_indices(), generate_qsite_inputs(), load_watermap_sites(), Returns {role: [sidechain_heavy_atom_indices]} for all Dream Team     catalytic (+3 more)

### Community 34 - "07 MD QMMM Defluorination (7)"
Cohesion: 0.20
Nodes (10): _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), parse_qsite_barrier(), parse_qsite_profile(), Per-scan-point energy series, returned in KCAL/MOL, from a Jaguar/QSite     rela, The constrained scan value actually used at each point, read from the Jaguar out, The departing fluorine's Mulliken charge at each scan point.      Jaguar prints (+2 more)

### Community 35 - "00_02 Project Utils (4)"
Cohesion: 0.22
Nodes (9): console_info(), console_separator(), console_title(), Logger, Canonical file logger for all pipeline scripts.      Creates a file-only handler, Bold section header — printed and optionally written to log file., Two-space-indented info line — printed and optionally written to log file., Horizontal rule — printed and optionally written to log file.      Parameters (+1 more)

### Community 36 - "00_03 Environment"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verifies installed packages in the current environment and prints a checklist., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 37 - "05 TopN and PDB Preparation (5)"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 38 - "06 SID Prime-MMGBSA (7)"
Cohesion: 0.25
Nodes (8): _avg_dg(), _boltzmann_mean_dg(), _dg_failures(), plot_mmgbsa_individual(), Log-sum-exp ("Boltzmann") ensemble mean binding free energy over frames., Headline per-job ΔG_bind estimator selected by CFG.MMGBSA_AVERAGING.      "mean", Frames whose ΔG_bind is a failed Prime minimisation rather than physics.      A, Per-job MM-GBSA: ΔG_bind against simulation time, plus its distribution.      Tw

### Community 39 - "06 SID Prime-MMGBSA (8)"
Cohesion: 0.29
Nodes (3): One in-place progress line for a sharded MM-GBSA run.      Aggregates across the, (frames read across all shards, shards currently in Prime minimisation)., ShardHeartbeat

### Community 40 - "00_02 Project Utils (5)"
Cohesion: 0.33
Nodes (7): compute_ligand_properties(), ensure_jobserver_on_working_disk(), _fs_device(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, st_dev of the nearest existing ancestor of `path` (identifies its mounted     fi, Force ALL Schrödinger job scratch onto the working disk so 06/07 run     entirel

### Community 41 - "00_02 Project Utils (6)"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 42 - "03 Validation Figures (7)"
Cohesion: 0.33
Nodes (6): _ext_match_03_style(), Bring an extended-analysis panel onto 03's typography before it is written., Persist a figure as a PNG at the pipeline's publication resolution, then free it, Create Tier_1A enzyme × ligand heatmap for lab validation targets., _xn_figure_06a(), _xn__save()

### Community 43 - "05 TopN and PDB Preparation (6)"
Cohesion: 0.33
Nodes (6): cif_to_pdb_gemmi(), generate_raw_step(), Injects metadata into the PDB Header.     Crucially, adds 'SOURCE_CIF' to track, Converts mmCIF to PDB using Gemmi.     - Moves Ligands (Non-Standard Residues) t, Task 1: Generate Raw PDB from CIF.     CIF is passed directly from 2_Best_Comple, update_pdb_header()

### Community 44 - "05 TopN and PDB Preparation (7)"
Cohesion: 0.40
Nodes (5): _parse_plip_xml(), _plip_cfg_cutoffs(), _plip_coo(), Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]., Parse a PLIP XML report into a contacts list compatible with _im_project /     _

### Community 46 - "00_02 Project Utils (7)"
Cohesion: 0.50
Nodes (4): calculate_angle(), calculate_burgi_dunitz(), Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit, Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for

### Community 47 - "02 Production (7)"
Cohesion: 0.50
Nodes (4): analyse_model_task(), Global worker routine functionally decoupled for robust serialisation capabiliti, Selects the representative pose across the Boltz diffusion samples tier-first: t, select_best_degrader_model()

### Community 48 - "02 Production (8)"
Cohesion: 0.50
Nodes (4): compute_cross_interface_pae(), load_extra_boltz_metrics(), Calculates the Predicted Aligned Error (PAE), focusing specifically on the prote, Parses background statistical validation tensors directly from Boltz NPZ file du

### Community 49 - "03 Validation Figures (8)"
Cohesion: 0.50
Nodes (4): _fig23_multitarget(), _fig23b_toptier_breakdown(), Figure 24: Top 25 multi-target proteins (stacked bar)., Figure 24b: for every protein that reaches the TOP tier (CFG.TIER_TOP) on at

### Community 50 - "05 TopN and PDB Preparation (8)"
Cohesion: 0.50
Nodes (4): _build_hd1_line(), enforce_catalytic_protonation(), Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev

## Knowledge Gaps
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `00_01 Project Config` to `05 TopN and PDB Preparation`, `04 Dendrogram`, `03 Validation Figures (2)`, `02 Production`, `03 Validation Figures (3)`, `07 MD QMMM Defluorination`, `02 Production (2)`, `06 SID Prime-MMGBSA`, `02 Production (3)`, `03 Validation Figures (5)`, `05 TopN and PDB Preparation (2)`, `06 SID Prime-MMGBSA (2)`, `07 MD QMMM Defluorination (2)`, `07 MD QMMM Defluorination (3)`, `06 SID Prime-MMGBSA (3)`, `03 Validation Figures (6)`, `07 MD QMMM Defluorination (4)`, `06 SID Prime-MMGBSA (5)`, `07 MD QMMM Defluorination (5)`, `06 SID Prime-MMGBSA (6)`, `07 MD QMMM Defluorination (6)`, `07 MD QMMM Defluorination (7)`, `06 SID Prime-MMGBSA (7)`, `03 Validation Figures (7)`, `03 Validation Figures (8)`?**
  _High betweenness centrality (0.645) - this node is a cross-community bridge._
- **Why does `main()` connect `02 Production` to `00_02 Project Utils`, `00_01 Project Config`, `02 Production (3)`, `02 Production (4)`, `02 Production (5)`?**
  _High betweenness centrality (0.121) - this node is a cross-community bridge._
- **Why does `process_single_job()` connect `07 MD QMMM Defluorination` to `07 MD QMMM Defluorination (6)`, `07 MD QMMM Defluorination (7)`, `00_02 Project Utils`, `00_01 Project Config`, `07 MD QMMM Defluorination (8)`, `07 MD QMMM Defluorination (2)`, `07 MD QMMM Defluorination (3)`, `07 MD QMMM Defluorination (4)`, `00_02 Project Utils (2)`, `07 MD QMMM Defluorination (5)`?**
  _High betweenness centrality (0.105) - this node is a cross-community bridge._
- **Are the 50 inferred relationships involving `CFG` (e.g. with `calculate_sn2_metrics()` and `check_catalytic_geometry()`) actually correct?**
  _`CFG` has 50 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `main()` (e.g. with `CFG` and `safe_name()`) actually correct?**
  _`main()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **What connects `===============================================================================`, `Central Configuration Repository — FAcDs Pipeline.      Attribute prefix → secti`, `Single source of truth for the holistic mechanistic score (0–1) — the raw geomet` to the rest of the system?**
  _354 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `05 TopN and PDB Preparation` be split into smaller, more focused modules?**
  _Cohesion score 0.08097165991902834 - nodes in this community are weakly interconnected._