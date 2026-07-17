# Graph Report - .  (2026-07-16)

## Corpus Check
- 13 files · ~60,000 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1007 nodes · 2011 edges · 62 communities (60 shown, 2 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 110 edges (avg confidence: 0.73)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 28|Community 28]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 32|Community 32]]
- [[_COMMUNITY_Community 33|Community 33]]
- [[_COMMUNITY_Community 34|Community 34]]
- [[_COMMUNITY_Community 35|Community 35]]
- [[_COMMUNITY_Community 36|Community 36]]
- [[_COMMUNITY_Community 37|Community 37]]
- [[_COMMUNITY_Community 38|Community 38]]
- [[_COMMUNITY_Community 39|Community 39]]
- [[_COMMUNITY_Community 40|Community 40]]
- [[_COMMUNITY_Community 41|Community 41]]
- [[_COMMUNITY_Community 42|Community 42]]
- [[_COMMUNITY_Community 43|Community 43]]
- [[_COMMUNITY_Community 44|Community 44]]
- [[_COMMUNITY_Community 45|Community 45]]
- [[_COMMUNITY_Community 46|Community 46]]
- [[_COMMUNITY_Community 47|Community 47]]
- [[_COMMUNITY_Community 48|Community 48]]
- [[_COMMUNITY_Community 49|Community 49]]
- [[_COMMUNITY_Community 50|Community 50]]
- [[_COMMUNITY_Community 51|Community 51]]
- [[_COMMUNITY_Community 52|Community 52]]
- [[_COMMUNITY_Community 53|Community 53]]
- [[_COMMUNITY_Community 54|Community 54]]
- [[_COMMUNITY_Community 55|Community 55]]
- [[_COMMUNITY_Community 56|Community 56]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 58|Community 58]]
- [[_COMMUNITY_Community 59|Community 59]]
- [[_COMMUNITY_Community 60|Community 60]]
- [[_COMMUNITY_Community 61|Community 61]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 63 edges
2. `main()` - 34 edges
3. `process_single_job()` - 30 edges
4. `main()` - 30 edges
5. `_echo()` - 25 edges
6. `prep_and_convert_phase()` - 24 edges
7. `console_info()` - 23 edges
8. `process_single_job()` - 21 edges
9. `_generate_comprehensive_figures_impl()` - 21 edges
10. `run_mmgbsa_phase()` - 20 edges

## Surprising Connections (you probably didn't know these)
- `apply_clean_spines()` --calls--> `clean_spines()`  [INFERRED]
  01_Merge_FAcDs.py → 00_02_Project_Utils_FAcDs.py
- `load_md_selected_jobs()` --indirect_call--> `CFG`  [INFERRED]
  05_TopN_and_PDB_Preparation_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `plot_pose_drift()` --indirect_call--> `CFG`  [INFERRED]
  05_TopN_and_PDB_Preparation_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `plot_machinery_distribution()` --indirect_call--> `CFG`  [INFERRED]
  05_TopN_and_PDB_Preparation_FAcDs.py → 00_01_Project_Config_FAcDs.py
- `plot_esp_alpha_carbon()` --indirect_call--> `CFG`  [INFERRED]
  05_TopN_and_PDB_Preparation_FAcDs.py → 00_01_Project_Config_FAcDs.py

## Import Cycles
- None detected.

## Communities (62 total, 2 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.08
Nodes (39): check_prep_needed(), _check_residue_identity_guard(), collect_best_cifs(), console_info(), generate_esp_charges(), get_source_tag(), index_existing_files(), load_catalytic_anchor_map() (+31 more)

### Community 1 - "Community 1"
Cohesion: 0.08
Nodes (36): _blockade_vec(), compute_wm_csv_stats(), console_qmm_ready(), _eaf_at(), find_eaf_file(), identity_from_cms(), kabsch_transform(), load_eaf_scalar_series() (+28 more)

### Community 2 - "Community 2"
Cohesion: 0.06
Nodes (37): angle_effective (Sidak multiplicity-corrected), Arg111/Arg114 Carboxylate Clamp, Backside SN2 Occlusion (SN2 dead-end B), Elite Bond-Strength Ceiling (TIER_ELITE_BDE_MAX 128), Bidentate Carboxylate Clamp Gate, Graded Chemistry Penalty (C-F BDE + occlusion), Bento & Bickelhaupt 2008 (backside SN2 sterics), O'Hagan 2008 (C-F bond strength) (+29 more)

### Community 3 - "Community 3"
Cohesion: 0.09
Nodes (25): _ConsoleRuleFilter, install_console_rule_filter(), Collapse stacked separator rules for the rest of this process's output., Simultaneous console + file logger shared by the pipeline steps.      The log pa, Record a line in the log file WITHOUT printing it to the console.          For o, A stdout wrapper that collapses consecutive separator rules.      The logs grow, ReportManager, clean_id() (+17 more)

### Community 4 - "Community 4"
Cohesion: 0.08
Nodes (35): analysis_worker_loop(), append_rows_to_csv(), atomic_to_csv(), console_info(), console_separator(), cpu_usage_summary(), extract_short_fasta_id(), flatten_job_result() (+27 more)

### Community 5 - "Community 5"
Cohesion: 0.08
Nodes (33): Nucleophile–substrate distance, resolving the several historical names.      The, Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present., Coerce a column to numeric, returning an empty series when it is absent., Compact, publication-style p-value formatting., No-op: panel letters are not drawn.      No other figure in this set carries an, A sample size that fits the ~0.5 in a tier occupies on a two-panel figure., The panel's test result, on its own line ABOVE the axes.      Inside the axes it (+25 more)

### Community 6 - "Community 6"
Cohesion: 0.08
Nodes (32): _auto_label_colour(), _fig23_multitarget(), _fig23b_toptier_breakdown(), _fig_folder03_dataset(), _fig_folder04_ai_confidence(), _fig_folder05_catalytic(), _fig_folder06_ligand(), _fig_path() (+24 more)

### Community 7 - "Community 7"
Cohesion: 0.12
Nodes (31): calculate_pareto_fronts(), perform_advanced_ranking(), DataFrame, Identifies non-dominated sorting fronts (Pareto Frontiers)., Multi-objective ranking of complexes.      Scales the feature matrix (robust, wi, Nucleophile–substrate distance, resolving the several historical names.      The, Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present. (+23 more)

### Community 8 - "Community 8"
Cohesion: 0.11
Nodes (29): _diag10_model_agreement(), _ext_match_03_style(), _fig24_sankey(), _fig25_pfas_size(), _fig_folder07_pfas(), generate_additional_figures(), generate_comprehensive_figures(), _generate_comprehensive_figures_impl() (+21 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (24): _cancel_launched_jobs(), _complex_label(), _echo(), _fail(), find_esp(), _hydrate(), _install_job_cleanup(), main() (+16 more)

### Community 10 - "Community 10"
Cohesion: 0.14
Nodes (24): _clean_job_env(), export_watermap_csv(), _log(), _md_msj(), _md_production(), _ok(), _phase_build(), _phase_md() (+16 more)

### Community 11 - "Community 11"
Cohesion: 0.08
Nodes (23): atomic_write_csv(), auto_label_colour(), calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner() (+15 more)

### Community 12 - "Community 12"
Cohesion: 0.09
Nodes (24): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+16 more)

### Community 13 - "Community 13"
Cohesion: 0.13
Nodes (20): check_job_status(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), heal_smiles_in_files(), _load_module(), mirror_best_cif() (+12 more)

### Community 14 - "Community 14"
Cohesion: 0.10
Nodes (22): analyse_model_task(), colabfold_a3m_path(), colabfold_meta_path(), compute_cross_interface_pae(), extract_3r3u_sequence(), load_extra_boltz_metrics(), map_active_site_residues(), Path (+14 more)

### Community 15 - "Community 15"
Cohesion: 0.12
Nodes (22): analyse_conflicts(), _aux_dir(), console_info(), generate_ramachandran_figures(), load_and_prep_data(), main(), _make_reporter(), The full statistical battery, registered into the same BH-corrected family as th (+14 more)

### Community 16 - "Community 16"
Cohesion: 0.16
Nodes (21): apply_figure_style(), clean_spines(), The pipeline's one typography and canvas definition, applied to matplotlib's rcP, Academic-style axes: remove top/right spines, thin the remaining borders.      C, console_separator(), format_job_label(), generate_comparative_residue_engagement(), generate_defluorination_landscape() (+13 more)

### Community 17 - "Community 17"
Cohesion: 0.13
Nodes (19): console_separator(), _esp_run(), extract_chain_l_mol(), load_reference_data(), main(), _plip_angles_ok(), plot_machinery_distribution(), Run a snippet under $SCHRODINGER/run — its interpreter, not ours. (+11 more)

### Community 18 - "Community 18"
Cohesion: 0.12
Nodes (20): _append_auxiliary_log(), _draw_interaction_diagram(), _im_don_acc(), _im_parse_pdb(), _im_project(), _im_render_diagram(), _im_ring(), _im_separate_atoms() (+12 more)

### Community 19 - "Community 19"
Cohesion: 0.16
Nodes (13): _fig_constants(), _fig_log(), load_metadata(), Return a dict of figure-generation constants from CFG., Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir., Unified logging for figure generation — routes through console_info., Handles verification and automated installation of visual tools., Check if PLIP is installed as a Python module. (+5 more)

### Community 20 - "Community 20"
Cohesion: 0.11
Nodes (20): _avg_dg(), _boltzmann_mean_dg(), _failure_windows(), _lookup_ligands(), _lookup_tiers(), _mmgbsa_csv(), _mmgbsa_dg_series(), _natural_rank() (+12 more)

### Community 21 - "Community 21"
Cohesion: 0.11
Nodes (20): _build_msj(), _csv_nonempty(), _diagnose_shard_failure(), discover_handover(), _lig_atom_indices(), _load_module(), merge_esp(), newest_ranked_csv() (+12 more)

### Community 22 - "Community 22"
Cohesion: 0.12
Nodes (20): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Add the `Frame` column to an MM-GBSA CSV that was written before stamping existe, Count frames written into a SID-out.eaf (the analysed result vector).      Parse, Ground-truth frame count of a Desmond trajectory via the Schrödinger     traj AP (+12 more)

### Community 23 - "Community 23"
Cohesion: 0.13
Nodes (16): Single source of truth for the holistic mechanistic score (0–1) — the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge() (+8 more)

### Community 24 - "Community 24"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 25 - "Community 25"
Cohesion: 0.12
Nodes (6): CFG, ===============================================================================, Central Configuration Repository — FAcDs Pipeline.      Attribute prefix → secti, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, generate_qsite_inputs(), Write a valid QSite/Jaguar QM/MM relaxed-scan .in for the SN2     dehalogenation

### Community 26 - "Community 26"
Cohesion: 0.12
Nodes (16): Remove all ANSI/VT100 escape sequences from a string., _strip_ansi(), calculate_min_distance(), check_md_equilibration(), console_info(), console_title(), extract_8residue_indices(), extract_hybrid_smart_system() (+8 more)

### Community 27 - "Community 27"
Cohesion: 0.15
Nodes (16): _avail_ram_gb(), _await_mmgbsa_jobserver(), _diagnose_mmgbsa_failure(), Record which trajectory frame each Prime row actually came from.      thermal_mm, Score the whole trajectory as concurrent contiguous frame shards.      Each shar, Block until no active job-server job whose name contains `job_prefix` remains (m, Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:     ret, Wall-clock length of a Desmond trajectory in ns, read from the frames' own     t (+8 more)

### Community 28 - "Community 28"
Cohesion: 0.17
Nodes (7): Heartbeat, Periodically report a long Schrödinger step's progress by tailing its log., Rewrite the single progress line in place with a carriage return., Close the open in-place \r line with a newline (phase change / exit)., Largest last-match across all configured patterns (phase-tolerant)., Best-effort (done, total) for the Prime phase from a ``*-prime*.log`` in, Highest 'Structure N (of M)' seen in the main log or any ``*-prime*.log``

### Community 29 - "Community 29"
Cohesion: 0.14
Nodes (16): _darken(), _engage_zones(), generate_reactive_pose_figures(), _load_reactive_pose_data(), _master_container(), plot_machinery_engagement(), plot_mmgbsa_decomposition(), The criteria that define engagement — every one of them a CFG constant. (+8 more)

### Community 30 - "Community 30"
Cohesion: 0.14
Nodes (16): Abel et al. 2008 (WaterMap thermodynamics), Bowers et al. 2006 (Desmond MD engine), Li et al. 2011 (MM-GBSA VSGB 2.0), Lu et al. 2021 / Roos et al. 2019 (OPLS4), Defluor_Propensity + Is_Defluorinating Verdict, Desmond MD (NPT, OPLS4, TIP3P), HID Base Tautomer Enforcement, Prime MM-GBSA Binding Free Energy (+8 more)

### Community 31 - "Community 31"
Cohesion: 0.17
Nodes (15): _box_diag(), calculate_improper_dihedral(), calculate_min_distance(), distance(), _ensure_box_3x3(), get_mic_vector(), mic_dists_2d(), ndarray (+7 more)

### Community 32 - "Community 32"
Cohesion: 0.41
Nodes (13): _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic(), _fig_14b_tt_interactions(), _fig_18b_tt_landscape(), ndarray, _tt_add_legend(), _tt_draw_scatter(), _tt_draw_stars() (+5 more)

### Community 33 - "Community 33"
Cohesion: 0.14
Nodes (15): _load_module(), load_triad_mapping(), plot_qsite_reaction_profile(), Path, _qsite_scan_failure_reason(), Loads catalytic triad residue numbers and static metrics from master CSV., Write an uncompressed QSite .mae trimmed to a solvation droplet: the full     pr, Launch the QSite/Jaguar executable on a generated .in, writing all output     in (+7 more)

### Community 34 - "Community 34"
Cohesion: 0.14
Nodes (15): Lightstone & Bruice 1996 / Bruice 2002 (NAC criteria), Mirdita et al. 2022 (ColabFold), ColabFold MSA Server, ESP QM Ligand Charges (Jaguar, default ON), Near Attack Conformation (NAC) Criteria, Phase 1 — High-Throughput Screening, Phase 2 — Select, Prepare & Filter, Screened-Pose-Is-Not-Simulated-Pose Drift (+7 more)

### Community 35 - "Community 35"
Cohesion: 0.16
Nodes (14): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), load_atoms_from_structure(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi (+6 more)

### Community 36 - "Community 36"
Cohesion: 0.18
Nodes (14): _block_bootstrap_median_ci(), _cliffs_delta(), _delta_word(), _dg_failures(), _draw_time_cumulative(), _effective_n(), plot_mmgbsa_combined(), Frames whose ΔG_bind is a failed Prime minimisation rather than physics.      A (+6 more)

### Community 37 - "Community 37"
Cohesion: 0.24
Nodes (11): _draw_rama_background(), Any, _rama_classify(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favoured, allowed, and outlier regions., Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes., Save a side-by-side comparison Ramachandran PNG. (+3 more)

### Community 38 - "Community 38"
Cohesion: 0.18
Nodes (11): graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit(), Constructs a three-dimensional RDKit molecule from a SMILES string.     This is, Extracts raw XYZ coordination matrices from RDKit molecule data blocks., Determines the optimal rotation/translation matrix required to superimpose spati, Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not b, FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perce (+3 more)

### Community 39 - "Community 39"
Cohesion: 0.20
Nodes (11): Beta-Fluorine/Ether Withdrawal Penalty, Khusnutdinova et al. 2023 (DFA is genuine substrate), ADONA (Ether-PFCA), C6O4 (Cyclic perfluoroether), Difluoroacetate (DFA, positive control), GenX (Ether-PFCA), PFOA (Perfluorooctanoic acid), PFOS (Perfluorooctanesulfonic acid) (+3 more)

### Community 40 - "Community 40"
Cohesion: 0.20
Nodes (10): _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), parse_qsite_barrier(), parse_qsite_profile(), Per-scan-point energy series, returned in KCAL/MOL, from a Jaguar/QSite     rela, The constrained scan value actually used at each point, read from the Jaguar out, The departing fluorine's Mulliken charge at each scan point.      Jaguar prints (+2 more)

### Community 41 - "Community 41"
Cohesion: 0.50
Nodes (7): _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_FAcDs.sh script, _sync_staging_log(), _tee(), WARN_EXIT_CODE

### Community 42 - "Community 42"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 43 - "Community 43"
Cohesion: 0.22
Nodes (9): Alpha-Carbon Reactive Centre, Asp110 Nucleophile, Boltz-2 Diffusion Co-Folding Model, Carbon-Fluorine Bond (~544 kJ/mol BDE), Passaro et al. 2025 (Boltz-2), FAcDs PFAS-27 Defluorination Pipeline, Mechanistic Filter (beyond structural confidence), PFAS Environmental Problem (+1 more)

### Community 44 - "Community 44"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 45 - "Community 45"
Cohesion: 0.29
Nodes (3): One in-place progress line for a sharded MM-GBSA run.      Aggregates across the, (frames read across all shards, shards currently in Prime minimisation)., ShardHeartbeat

### Community 46 - "Community 46"
Cohesion: 0.29
Nodes (8): Asp134 Acid Catalyst (Asp-His dyad), Catalytic Triad (Asp110/Asp134/His277), Chan et al. 2011 JACS (3R3U FAcD structure), PDB 3R3U Crystal Reference, Fluoroacetate Dehalogenase (FAcD), His277 Base Catalyst, Monomeric Active-Site Scope Limitation, Residue Numbering +2 Offset (predicted vs crystal)

### Community 47 - "Community 47"
Cohesion: 0.29
Nodes (7): console_info(), console_separator(), console_title(), Logger, Bold section header — printed and optionally written to log file., Two-space-indented info line — printed and optionally written to log file., Horizontal rule — printed and optionally written to log file.      Parameters

### Community 48 - "Community 48"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 49 - "Community 49"
Cohesion: 0.33
Nodes (6): fetch_msa_direct(), Generates a stable hash for FASTA sequences to verify consistency during resumed, Validates the A3M MSA file to detect corruption. Returns True if the file struct, Fetches the MSA directly from the ColabFold REST API. This circumvents the Boltz, sequence_hash(), validate_a3m_file()

### Community 50 - "Community 50"
Cohesion: 0.33
Nodes (6): cif_to_pdb_gemmi(), generate_raw_step(), Injects metadata into the PDB Header.     Crucially, adds 'SOURCE_CIF' to track, Converts mmCIF to PDB using Gemmi.     - Moves Ligands (Non-Standard Residues) t, Task 1: Generate Raw PDB from CIF.     CIF is passed directly from 2_Best_Comple, update_pdb_header()

### Community 51 - "Community 51"
Cohesion: 0.33
Nodes (6): Becke 1993 / Lee-Yang-Parr 1988 (B3LYP), Murphy et al. 2000 (QSite implementation), Yue et al. 2021 (FAcD SN2 QM/MM energetics), Electronic (dE) not Free-Energetic (dG) Barrier Limitation, QSite QM/MM Relaxed Scan (B3LYP/6-31G**), Rationale: B3LYP required (frozen-orbital cuts reject M06-2X/B3LYP-D3)

### Community 52 - "Community 52"
Cohesion: 0.40
Nodes (5): compute_ligand_properties(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, Canonical file logger for all pipeline scripts.      Creates a file-only handler, setup_logging()

### Community 53 - "Community 53"
Cohesion: 0.40
Nodes (5): _parse_plip_xml(), _plip_cfg_cutoffs(), _plip_coo(), Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]., Parse a PLIP XML report into a contacts list compatible with _im_project /     _

### Community 55 - "Community 55"
Cohesion: 0.40
Nodes (5): 00_01 Project Config (CFG single source of truth), Code Architecture Graph (826 nodes, CFG hub), 00_03 Environment (conda/pip pinning), Missing-Is-Not-Zero Defence (SENTINEL_UNDEFINED), Rationale: missing measurement must not become optimal value

### Community 56 - "Community 56"
Cohesion: 0.50
Nodes (4): calculate_angle(), calculate_burgi_dunitz(), Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit, Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for

### Community 57 - "Community 57"
Cohesion: 0.50
Nodes (4): analyse_candidate_structure(), load_structure_safe(), Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilit, Superimposes the candidate structure onto the DeHa4 control model.     Computes

### Community 58 - "Community 58"
Cohesion: 0.50
Nodes (4): init_gpu_reservations(), query_gpu_memory(), A resilient system query that safely falls back to CPU memory context if nvidia-, Maps unallocated VRAM memory segments natively upon system initialisation.

### Community 59 - "Community 59"
Cohesion: 0.50
Nodes (4): _build_hd1_line(), enforce_catalytic_protonation(), Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev

### Community 60 - "Community 60"
Cohesion: 0.50
Nodes (4): 6 Control Cases (TFA/FA/DFA + 3R3U + DEHA4), Farajollahi et al. 2024 (DEHA4 validation), DEHA4 (Delftia acidovorans D4B), Fluoroacetate (FA, positive control)

## Knowledge Gaps
- **58 isolated node(s):** `WARN_EXIT_CODE`, `git_push_FAcDs.sh script`, `Trifluoroacetate (TFA, alpha-CF3 substrate)`, `PFOA (Perfluorooctanoic acid)`, `PFOS (Perfluorooctanesulfonic acid)` (+53 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `Community 25` to `Community 0`, `Community 1`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 8`, `Community 14`, `Community 15`, `Community 16`, `Community 17`, `Community 20`, `Community 22`, `Community 23`, `Community 27`, `Community 29`, `Community 33`, `Community 36`, `Community 40`?**
  _High betweenness centrality (0.513) - this node is a cross-community bridge._
- **Why does `process_single_job()` connect `Community 1` to `Community 33`, `Community 40`, `Community 11`, `Community 16`, `Community 54`, `Community 25`, `Community 26`, `Community 31`?**
  _High betweenness centrality (0.096) - this node is a cross-community bridge._
- **Why does `main()` connect `Community 4` to `Community 11`, `Community 13`, `Community 14`, `Community 49`, `Community 25`, `Community 58`, `Community 57`?**
  _High betweenness centrality (0.096) - this node is a cross-community bridge._
- **Are the 48 inferred relationships involving `CFG` (e.g. with `calculate_sn2_metrics()` and `check_catalytic_geometry()`) actually correct?**
  _`CFG` has 48 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `main()` (e.g. with `CFG` and `safe_name()`) actually correct?**
  _`main()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Initialises a dual-handler logger (File + Console).     Log is saved beside the`, `Normalisation:     1. Convert to String     2. Upper Case     3. Remove Gaps (-)`, `Quality Control (QC):     Checks for internal stops and excessive ambiguous resi` to the rest of the system?**
  _436 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.08097165991902834 - nodes in this community are weakly interconnected._