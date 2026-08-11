# Graph Report - .  (2026-08-11)

## Corpus Check
- Large corpus: 21 files · ~514,705 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 1236 nodes · 2662 edges · 67 communities (64 shown, 3 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 191 edges (avg confidence: 0.73)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_07 QMMM console & geometry|07 QM/MM console & geometry]]
- [[_COMMUNITY_Literature citations|Literature citations]]
- [[_COMMUNITY_03 Comparative figures|03 Comparative figures]]
- [[_COMMUNITY_05 Structure preparation|05 Structure preparation]]
- [[_COMMUNITY_07 Figure styling (SSOT)|07 Figure styling (SSOT)]]
- [[_COMMUNITY_03 Pareto ranking|03 Pareto ranking]]
- [[_COMMUNITY_06 Analysis labelling|06 Analysis labelling]]
- [[_COMMUNITY_02 ColabFold MSA IO|02 ColabFold MSA I/O]]
- [[_COMMUNITY_06 Progress heartbeat|06 Progress heartbeat]]
- [[_COMMUNITY_07 Reactive-pose figures|07 Reactive-pose figures]]
- [[_COMMUNITY_03 Top-tier figures|03 Top-tier figures]]
- [[_COMMUNITY_07 QSite job control|07 QSite job control]]
- [[_COMMUNITY_07 QSite scan parsing|07 QSite scan parsing]]
- [[_COMMUNITY_03 Column resolvers|03 Column resolvers]]
- [[_COMMUNITY_03 Logging & file discovery|03 Logging & file discovery]]
- [[_COMMUNITY_03 Journal figures|03 Journal figures]]
- [[_COMMUNITY_05 ESP charge prep|05 ESP charge prep]]
- [[_COMMUNITY_CFG tier definitions|CFG tier definitions]]
- [[_COMMUNITY_06 MD job build|06 MD job build]]
- [[_COMMUNITY_02 Job status & format|02 Job status & format]]
- [[_COMMUNITY_05 Figure constants|05 Figure constants]]
- [[_COMMUNITY_06 Defluorination geometry|06 Defluorination geometry]]
- [[_COMMUNITY_03 SN2 geometry|03 SN2 geometry]]
- [[_COMMUNITY_02 Mechanistic scoring|02 Mechanistic scoring]]
- [[_COMMUNITY_06 Desmond job management|06 Desmond job management]]
- [[_COMMUNITY_06 MM-GBSA jobserver|06 MM-GBSA jobserver]]
- [[_COMMUNITY_06 Statistics & CI|06 Statistics & CI]]
- [[_COMMUNITY_01 Sequence mergeQC|01 Sequence merge/QC]]
- [[_COMMUNITY_04 Phylogenydendrogram|04 Phylogeny/dendrogram]]
- [[_COMMUNITY_Utils nucleophile geometry|Utils: nucleophile geometry]]
- [[_COMMUNITY_06 WaterMap & MD phases|06 WaterMap & MD phases]]
- [[_COMMUNITY_06 Progress ETA|06 Progress ETA]]
- [[_COMMUNITY_00_00 Pipeline runner|00_00 Pipeline runner]]
- [[_COMMUNITY_03 AI-confidence figures|03 AI-confidence figures]]
- [[_COMMUNITY_06 MM-GBSA plots|06 MM-GBSA plots]]
- [[_COMMUNITY_Utils PBC distance math|Utils: PBC distance math]]
- [[_COMMUNITY_07 Label colour utils|07 Label colour utils]]
- [[_COMMUNITY_06 Job cleanup guards|06 Job cleanup guards]]
- [[_COMMUNITY_Utils console logging|Utils: console logging]]
- [[_COMMUNITY_Utils Ramachandran|Utils: Ramachandran]]
- [[_COMMUNITY_02 Pi-interaction analysis|02 Pi-interaction analysis]]
- [[_COMMUNITY_07 Defluorination verdict|07 Defluorination verdict]]
- [[_COMMUNITY_Architecture diagram|Architecture diagram]]
- [[_COMMUNITY_02 RDKit structure mapping|02 RDKit structure mapping]]
- [[_COMMUNITY_05 Interaction diagrams|05 Interaction diagrams]]
- [[_COMMUNITY_03 Per-model variance|03 Per-model variance]]
- [[_COMMUNITY_Utils file logger|Utils: file logger]]
- [[_COMMUNITY_00_03 Environment setup|00_03 Environment setup]]
- [[_COMMUNITY_Utils ligand properties|Utils: ligand properties]]
- [[_COMMUNITY_05 Contact geometry|05 Contact geometry]]
- [[_COMMUNITY_Utils dihedral angles|Utils: dihedral angles]]
- [[_COMMUNITY_02 Analysis worker IO|02 Analysis worker I/O]]
- [[_COMMUNITY_02 MSA fetch & validate|02 MSA fetch & validate]]
- [[_COMMUNITY_05 PDB donoracceptor|05 PDB donor/acceptor]]
- [[_COMMUNITY_05 PLIP parsing|05 PLIP parsing]]
- [[_COMMUNITY_07 Lazy trajectory|07 Lazy trajectory]]
- [[_COMMUNITY_Workflow phases (logo)|Workflow phases (logo)]]
- [[_COMMUNITY_Utils Burgi-Dunitz angle|Utils: Burgi-Dunitz angle]]
- [[_COMMUNITY_02 Structure loading|02 Structure loading]]
- [[_COMMUNITY_02 Model analysis task|02 Model analysis task]]
- [[_COMMUNITY_02 PAEBoltz metrics|02 PAE/Boltz metrics]]
- [[_COMMUNITY_02 GPU reservation|02 GPU reservation]]
- [[_COMMUNITY_02 Active-site mapping|02 Active-site mapping]]
- [[_COMMUNITY_05 His protonation|05 His protonation]]
- [[_COMMUNITY_DeFluorX logobrand|DeFluorX logo/brand]]
- [[_COMMUNITY_README overview|README overview]]
- [[_COMMUNITY_git push helper|git push helper]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 91 edges
2. `main()` - 36 edges
3. `process_single_job()` - 35 edges
4. `main()` - 34 edges
5. `_generate_comprehensive_figures_impl()` - 32 edges
6. `console_info()` - 32 edges
7. `_echo()` - 31 edges
8. `prep_and_convert_phase()` - 25 edges
9. `deflx_fig_name()` - 22 edges
10. `run_mmgbsa_phase()` - 22 edges

## Surprising Connections (you probably didn't know these)
- `generate_plots()` --indirect_call--> `CFG`  [INFERRED]
  01_Merge_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `calculate_sn2_metrics()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `check_catalytic_geometry()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `generate_scientific_ranking_csv()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `main()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **DeFluorX four-phase data flow (01 to 07)** — step_01_merge, step_02_production, step_03_validation_figures, step_05_topn_prep, step_06_physics, step_07_qmmm [EXTRACTED 0.90]
- **Foundation (CFG/utils/env) imported by every step** — step_0001_config, step_0002_utils, step_0003_env, step_02_production, step_06_physics [EXTRACTED 0.85]
- **NAC/SN2 mechanistic gating machinery** — readme_nac, readme_sn2_walden, readme_catalytic_triad, readme_fluoride_cradle, readme_mech_score [EXTRACTED 0.85]
- **Foundation modules imported by every step** — arch_svg_00_01_config, arch_svg_00_02_utils, arch_svg_00_03_environment [EXTRACTED 0.90]
- **Main screening-to-defluorination pipeline chain** — arch_svg_02_production, arch_svg_03_validation_figures, arch_svg_05_topn_prep, arch_svg_06_physics_validation, arch_svg_07_qmmm_defluorination [EXTRACTED 0.85]

## Communities (67 total, 3 thin omitted)

### Community 0 - "07 QM/MM console & geometry"
Cohesion: 0.05
Nodes (52): Remove all ANSI/VT100 escape sequences from a string., _strip_ansi(), _blockade_vec(), calculate_min_distance(), check_md_equilibration(), compute_wm_csv_stats(), console_info(), _eaf_at() (+44 more)

### Community 1 - "Literature citations"
Cohesion: 0.06
Nodes (50): Auffinger et al. 2004 PNAS (halogen bonds), Passaro et al. 2025 Boltz-2, Burgi-Dunitz 1974 Tetrahedron (attack trajectory), Chan et al. 2011 JACS (3R3U, enzymatic defluorination), Donald 2011 (salt-bridge geometry), Farajollahi et al. 2024 ACS Omega (DeHa4/D4B), Kurihara & Esaki 2008 (haloalkanoate dehalogenase), Michaud-Agrawal et al. 2011 MDAnalysis (+42 more)

### Community 2 - "03 Comparative figures"
Cohesion: 0.10
Nodes (43): _control_star_df(), _ctrl_star(), _diag10_model_agreement(), _fig23_multitarget(), _fig24_sankey(), _fig25_pfas_size(), _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic() (+35 more)

### Community 3 - "05 Structure preparation"
Cohesion: 0.08
Nodes (42): check_prep_needed(), _check_residue_identity_guard(), cif_to_pdb_gemmi(), collect_best_cifs(), console_info(), generate_esp_charges(), generate_raw_step(), get_source_tag() (+34 more)

### Community 4 - "07 Figure styling (SSOT)"
Cohesion: 0.07
Nodes (40): apply_figure_style(), clean_spines(), deflx_fig_name(), Return a figure path with its suffix swapped to the ACTIVE SSOT format (CFG.VIS_, The pipeline's one typography and canvas definition, applied to matplotlib's rcP, Academic-style axes: remove top/right spines, thin the remaining borders.      C, _active_site_series(), _draw_trajectory_figures() (+32 more)

### Community 5 - "03 Pareto ranking"
Cohesion: 0.10
Nodes (37): calculate_pareto_fronts(), _is_control_mask(), perform_advanced_ranking(), DataFrame, Nucleophile–substrate distance, resolving the several historical names.      The, Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present., Coerce a column to numeric, returning an empty series when it is absent. (+29 more)

### Community 6 - "06 Analysis labelling"
Cohesion: 0.09
Nodes (37): _analysis_dir(), _complex_label(), _ctrl_label(), _ctrl_palette(), _echo(), _emit_timings(), _fail(), find_esp() (+29 more)

### Community 7 - "02 ColabFold MSA I/O"
Cohesion: 0.10
Nodes (35): atomic_to_csv(), colabfold_a3m_path(), colabfold_meta_path(), console_info(), console_separator(), extract_3r3u_sequence(), generate_scientific_ranking_csv(), load_cached_alignments() (+27 more)

### Community 8 - "06 Progress heartbeat"
Cohesion: 0.08
Nodes (17): _close_bar(), Heartbeat, MDHeartbeat, _progress_line(), (current relaxation stage, production stage number) for the pre-production phase, Refresh the single shared \r progress line (terminal only, never the log)., Close the open \r progress line with a newline (phase change / completion)., Periodically report a long Schrödinger step's progress by tailing its log. (+9 more)

### Community 9 - "07 Reactive-pose figures"
Cohesion: 0.09
Nodes (32): _draw_reactive_pose_for_rank(), generate_reactive_pose_figures(), _load_module(), _load_ranked_df(), _load_reactive_pose_data(), _mae_av(), _mae_parse_cts(), _mech_ligand_geom() (+24 more)

### Community 10 - "03 Top-tier figures"
Cohesion: 0.09
Nodes (29): _auto_label_colour(), _ext_match_03_style(), _fig23b_toptier_breakdown(), _lig_short(), _load_module(), ndarray, Place a boxed statistics annotation on an axis., Auto-contrast label colour for hex_bg (delegates to _auto_label_colour). (+21 more)

### Community 11 - "07 QSite job control"
Cohesion: 0.10
Nodes (29): console_qmm_ready(), console_separator(), console_title(), _kill_all_qsite_jobs(), _kill_qsite_job(), load_triad_mapping(), main(), _mask_oomd_at_start() (+21 more)

### Community 12 - "07 QSite scan parsing"
Cohesion: 0.10
Nodes (30): _extract_cf_distance_series(), _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), _extract_scan_points(), _parse_qsite_pes_raw(), parse_qsite_profile(), _qsite_frame_summary() (+22 more)

### Community 13 - "03 Column resolvers"
Cohesion: 0.09
Nodes (29): Return the first candidate column that exists in ``df`` (else ``None``)., Resolve a headline-pillar name to whichever concrete column is present., Compact, publication-style p-value formatting., No-op: panel letters are not drawn.      No other figure in this set carries an, A sample size that fits the ~0.5 in a tier occupies on a two-panel figure., The panel's test result, on its own line ABOVE the axes.      Inside the axes it, A per-group random subsample for strip overlays (see _xn_STRIP_MAX_PER_GROUP)., Kruskal-Wallis across the ordered groups, formatted for an on-panel annotation. (+21 more)

### Community 14 - "03 Logging & file discovery"
Cohesion: 0.10
Nodes (28): latest_by_mtime(), Simultaneous console + file logger shared by the pipeline steps.      The log pa, The newest existing file by modification time, or None.      Selection is by mti, ReportManager, analyse_conflicts(), _aux_dir(), console_info(), generate_comprehensive_figures() (+20 more)

### Community 15 - "03 Journal figures"
Cohesion: 0.12
Nodes (28): _generate_comprehensive_figures_impl(), _jf_box(), _jf_circos(), _jf_fit_fs(), _jf_importance(), _jf_ligshort(), _jf_manhattan(), _jf_metric_matrix() (+20 more)

### Community 16 - "05 ESP charge prep"
Cohesion: 0.10
Nodes (26): console_separator(), _esp_run(), extract_chain_l_mol(), _is_reference_control(), _lig_short_token(), _load_module(), load_reference_data(), main() (+18 more)

### Community 17 - "CFG tier definitions"
Cohesion: 0.08
Nodes (19): CFG, ===============================================================================, Central Configuration Repository - DeFluorX Pipeline.      Attribute prefix → se, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, _aggregate_mmgbsa_stats(), plot_qsite_ensemble_profiles(), plot_qsite_reaction_profile(), _qsite_ensemble_barrier() (+11 more)

### Community 18 - "06 MD job build"
Cohesion: 0.08
Nodes (27): _build_msj(), _csv_nonempty(), _diagnose_shard_failure(), discover_handover(), _lig_atom_indices(), _load_module(), merge_esp(), _mmgbsa_complete() (+19 more)

### Community 19 - "02 Job status & format"
Cohesion: 0.10
Nodes (24): check_job_status(), cpu_usage_summary(), extract_short_fasta_id(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), heal_smiles_in_files() (+16 more)

### Community 20 - "05 Figure constants"
Cohesion: 0.12
Nodes (19): _append_auxiliary_log(), _fig_constants(), _fig_log(), load_metadata(), Return a dict of figure-generation constants from CFG., Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir., Unified logging for figure generation - routes through console_info., Reads a specific tool's output log and appends it to the main report. (+11 more)

### Community 21 - "06 Defluorination geometry"
Cohesion: 0.11
Nodes (25): _defl_ligand_atoms(), _defl_mapped_residues(), _defl_min_image(), _defl_resid(), _md_msj(), _md_production(), _mmgbsa_csv(), newest_ranked_csv() (+17 more)

### Community 22 - "03 SN2 geometry"
Cohesion: 0.09
Nodes (24): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+16 more)

### Community 23 - "02 Mechanistic scoring"
Cohesion: 0.11
Nodes (18): Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge() (+10 more)

### Community 24 - "06 Desmond job management"
Cohesion: 0.11
Nodes (22): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Classify every desmond_md_job_R_* directory and return the subset     that still, Step 1: event_analysis.py - generates the SID-in.eaf descriptor., Step 2: analyze_simulation.py - generates the SID-out.eaf result vector.      Us (+14 more)

### Community 25 - "06 MM-GBSA jobserver"
Cohesion: 0.11
Nodes (20): _avail_ram_gb(), _await_mmgbsa_jobserver(), _cleanup_mmgbsa_shards(), _diagnose_mmgbsa_failure(), _free_swap_gb(), Remove the per-job _MMGBSA_Shards scratch once the merged CSV is settled.      T, Best-effort human-readable cause when MM-GBSA exits non-zero, read from the, Free RAM the kernel expects to hand out without swapping (MemAvailable). (+12 more)

### Community 26 - "06 Statistics & CI"
Cohesion: 0.12
Nodes (20): _avg_dg(), _block_bootstrap_median_ci(), _boltzmann_mean_dg(), _cliffs_delta(), _delta_word(), _dg_failures(), _draw_time_cumulative(), _effective_n() (+12 more)

### Community 27 - "01 Sequence merge/QC"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 28 - "04 Phylogeny/dendrogram"
Cohesion: 0.18
Nodes (18): clean_id(), console_info(), console_separator(), generate_phylogenies(), generate_upgma_newick(), get_kmer_counts(), _load_module(), main() (+10 more)

### Community 29 - "Utils: nucleophile geometry"
Cohesion: 0.11
Nodes (17): calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner(), ===============================================================================, Geometry-based Oδ fallback nucleophile search for the Smart-Lock algorithm. (+9 more)

### Community 30 - "06 WaterMap & MD phases"
Cohesion: 0.16
Nodes (18): _clean_job_env(), export_watermap_csv(), _log(), _ok(), _phase_build(), _phase_md(), _phase_watermap(), A fresh login-shell environment for launching WaterMap, detached from the $SCHRO (+10 more)

### Community 31 - "06 Progress ETA"
Cohesion: 0.12
Nodes (13): _eta_str(), _extract_with_progress(), _fmt_dur(), Refresh the shared single \r progress line from a CPU worker (extraction / SID /, Extract members from tgz into wd (dropping one leading path component), drawing, Put the production trajectory at the job-dir root as {job}_trj + {job}.ene and r, Seconds → compact H/M/S (2h07m03s · 8m12s · 41s)., A ' · ETA <dur>' suffix from a linear extrapolation of the current rate - empty (+5 more)

### Community 32 - "00_00 Pipeline runner"
Cohesion: 0.22
Nodes (13): _cancel_schrodinger_jobs(), _ensure_jobserver(), _print_step_menu(), _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_DeFluorX.sh script, _sync_staging_log() (+5 more)

### Community 33 - "03 AI-confidence figures"
Cohesion: 0.12
Nodes (16): _fig_folder04_ai_confidence(), _kruskal_by_tier(), _median_ci95(), Folder 04_AI_Confidence_Quality - Boltz-2 confidence + pTM/ipTM., Draw a boxed statistics annotation inside an axis (no effect if text is empty)., Register one test into the family. extra carries effect size, group ns, medians, The full statistical battery, registered into the same BH-corrected family as th, Kruskal–Wallis across tiers with epsilon-squared effect size. Returns annotation (+8 more)

### Community 34 - "06 MM-GBSA plots"
Cohesion: 0.14
Nodes (16): _failure_windows(), _lookup_controls(), _lookup_ligands(), _lookup_tiers(), _mmgbsa_dg_series(), _plot_rank_mmgbsa(), Extract the per-frame ΔG_bind series, tolerant of column-name variants., Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier     co (+8 more)

### Community 35 - "Utils: PBC distance math"
Cohesion: 0.17
Nodes (15): _box_diag(), calculate_improper_dihedral(), calculate_min_distance(), distance(), _ensure_box_3x3(), get_mic_vector(), mic_dists_2d(), ndarray (+7 more)

### Community 36 - "07 Label colour utils"
Cohesion: 0.18
Nodes (14): auto_label_colour(), The contrast colour for a label written ON a filled mark, from the fill's lumina, _darken(), _engage_zones(), _master_container(), plot_machinery_engagement(), plot_mmgbsa_decomposition(), The criteria that define engagement - every one of them a CFG constant. (+6 more)

### Community 37 - "06 Job cleanup guards"
Cohesion: 0.15
Nodes (10): _cancel_launched_jobs(), _install_job_cleanup(), _kill_desmond_jobs(), _kill_launched_procs(), OomdGuard, Cancel every still-registered job on the job server. Idempotent (guarded), safe, Force-terminate any Desmond MD process chain still on the GPU for a job THIS run, Terminate every still-registered LOCAL subprocess (SID, MM-GBSA drivers) and its (+2 more)

### Community 38 - "Utils: console logging"
Cohesion: 0.18
Nodes (5): _ConsoleRuleFilter, install_console_rule_filter(), A stdout wrapper that collapses consecutive separator rules.      The logs grow, Collapse stacked separator rules for the rest of this process's output., Record a line in the log file WITHOUT printing it to the console.          For o

### Community 39 - "Utils: Ramachandran"
Cohesion: 0.22
Nodes (13): _draw_rama_background(), Any, _rama_classify(), _rama_palette(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favoured, allowed, and outlier regions., Ramachandran colours from CFG.VIS_RAMA (the single source of truth); the built-i (+5 more)

### Community 40 - "02 Pi-interaction analysis"
Cohesion: 0.20
Nodes (12): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi, Classify a ligand atom as 'anion'- or 'cation'-capable for salt-bridge     detec (+4 more)

### Community 41 - "07 Defluorination verdict"
Cohesion: 0.17
Nodes (12): _collect_qsite_results(), defluor_propensity(), defluor_verdict(), parse_qsite_barrier(), _qsite_qm_water_count(), Re-derive the full QM-region provenance from a job's QSite .in for the scan CSV', Scalar barrier/reaction-energy/fluoride subset of parse_qsite_profile, for     t, QM waters in a frame's QSite input = full-QM molids beyond the ligand (molid 2) (+4 more)

### Community 42 - "Architecture diagram"
Cohesion: 0.21
Nodes (12): 00_00 Run Pipeline Orchestrator (bash), 00_01 Project Config (SSOT), 00_02 Project Utils, 00_03 Environment, 01 Merge, 02 Production (Boltz-2 Screening), 03 Validation Figures, 04 Dendrogram (+4 more)

### Community 43 - "02 RDKit structure mapping"
Cohesion: 0.18
Nodes (11): graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit(), Constructs a three-dimensional RDKit molecule from a SMILES string.     This is, Extracts raw XYZ coordination matrices from RDKit molecule data blocks., Determines the optimal rotation/translation matrix required to superimpose spati, Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not b, FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perce (+3 more)

### Community 44 - "05 Interaction diagrams"
Cohesion: 0.18
Nodes (11): _draw_interaction_diagram(), _im_project(), _im_render_diagram(), _im_resnum(), _im_role_resnums(), _im_separate_atoms(), Shared 2D interaction renderer for InteractionMap (mode='distance') and PLIP (mo, {job_name: {residue_number: role_key}} from the ranked CSV Mapped_* columns, so (+3 more)

### Community 45 - "03 Per-model variance"
Cohesion: 0.22
Nodes (10): _mm_load_variance_df(), _mm_nac_lines(), _mm_variance_by_model(), _mm_variance_cut(), _mm_variance_iptm90(), Per-model variance CSV (one row per model per complex) joined to the pipeline ti, Strict-NAC reference lines shared by the multi-model panels (CFG-sourced)., Fig 13 - per-Boltz-model SN2 distance/angle spread by tier (are the 5 models con (+2 more)

### Community 46 - "Utils: file logger"
Cohesion: 0.22
Nodes (9): console_info(), console_separator(), console_title(), Logger, Canonical file logger for all pipeline scripts.      Creates a file-only handler, Bold section header - printed and optionally written to log file., Two-space-indented info line - printed and optionally written to log file., Horizontal rule - printed and optionally written to log file.      Parameters (+1 more)

### Community 47 - "00_03 Environment setup"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 48 - "Utils: ligand properties"
Cohesion: 0.25
Nodes (7): atomic_write_csv(), compute_ligand_properties(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, Write a JSON file so a reader never sees a half-written one.      The write goes, Write a DataFrame to CSV so a reader never sees a half-written file.      Same g, write_json_atomic()

### Community 49 - "05 Contact geometry"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 50 - "Utils: dihedral angles"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 51 - "02 Analysis worker I/O"
Cohesion: 0.33
Nodes (6): analysis_worker_loop(), append_rows_to_csv(), flatten_job_result(), Flattens a process_single_job result dict into a single-level row for CSV writin, Appends a list of result rows to the master CSV, deduplicates, sorts by job_name, Background worker that analyses freshly GPU-predicted jobs and appends rows to t

### Community 52 - "02 MSA fetch & validate"
Cohesion: 0.33
Nodes (6): fetch_msa_direct(), Generates a stable hash for FASTA sequences to verify consistency during resumed, Validates the A3M MSA file to detect corruption. Returns True if the file struct, Fetches the MSA directly from the ColabFold REST API. This circumvents the Boltz, sequence_hash(), validate_a3m_file()

### Community 53 - "05 PDB donor/acceptor"
Cohesion: 0.33
Nodes (6): _im_don_acc(), _im_parse_pdb(), _im_ring(), Split heavy N/O atoms into donors (with bonded H) and acceptors.      Returns (d, Aromatic ring (centroid, normal) for an aromatic residue, else None., Return (lig_atoms, pro_residues, lig_hb) from a prepared (protonated) PDB.

### Community 54 - "05 PLIP parsing"
Cohesion: 0.40
Nodes (5): _parse_plip_xml(), _plip_cfg_cutoffs(), _plip_coo(), Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]., Parse a PLIP XML report into a contacts list compatible with _im_project /     _

### Community 56 - "Workflow phases (logo)"
Cohesion: 0.40
Nodes (5): F- release (C-F bond cleavage payoff), Phase 1 - Predict (Boltz-2 MSA + 5-model diffusion), Phase 2 - Screen & Tier (competence-floor tiers, SN2 geometry, Pareto fronts), Phase 3 - Simulate (Desmond MD 1 us, SID, Prime MM-GBSA, WaterMap, NAC dwell), Phase 4 - Defluorinate (B3LYP QM/MM relaxed scan, dE barrier, Mulliken F- charge)

### Community 57 - "Utils: Burgi-Dunitz angle"
Cohesion: 0.50
Nodes (4): calculate_angle(), calculate_burgi_dunitz(), Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for, Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit

### Community 58 - "02 Structure loading"
Cohesion: 0.50
Nodes (4): analyse_candidate_structure(), load_structure_safe(), Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilit, Superimposes the candidate structure onto the DeHa4 control model.     Computes

### Community 59 - "02 Model analysis task"
Cohesion: 0.50
Nodes (4): analyse_model_task(), Global worker routine functionally decoupled for robust serialisation capabiliti, Selects the representative pose across the Boltz diffusion samples tier-first: t, select_best_degrader_model()

### Community 60 - "02 PAE/Boltz metrics"
Cohesion: 0.50
Nodes (4): compute_cross_interface_pae(), load_extra_boltz_metrics(), Calculates the Predicted Aligned Error (PAE), focusing specifically on the prote, Parses background statistical validation tensors directly from Boltz NPZ file du

### Community 61 - "02 GPU reservation"
Cohesion: 0.50
Nodes (4): init_gpu_reservations(), query_gpu_memory(), A resilient system query that safely falls back to CPU memory context if nvidia-, Maps unallocated VRAM memory segments natively upon system initialisation.

### Community 62 - "02 Active-site mapping"
Cohesion: 0.50
Nodes (4): map_active_site_residues(), Append a newly computed alignment to the on-disk cache file, under an exclusive, Performs a global sequence alignment between the canonical reference sequence, update_alignment_stats()

### Community 63 - "05 His protonation"
Cohesion: 0.50
Nodes (4): _build_hd1_line(), enforce_catalytic_protonation(), Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev

### Community 64 - "DeFluorX logo/brand"
Cohesion: 0.50
Nodes (4): DeFluorX, 7-step automated pipeline, 44,543 lines, Boltz-2 -> Desmond MD -> B3LYP QM/MM, CFG single source of truth, 2,150 FAcD variants x 27 PFAS = 58,050 predicted complexes, Discovering FAcDs for PFAS Defluorination (fully automated end-to-end enzyme-discovery pipeline)

## Knowledge Gaps
- **32 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `git_push_DeFluorX.sh script` (+27 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `CFG tier definitions` to `07 QM/MM console & geometry`, `03 Comparative figures`, `05 Structure preparation`, `07 Figure styling (SSOT)`, `06 Analysis labelling`, `02 ColabFold MSA I/O`, `07 Reactive-pose figures`, `03 Top-tier figures`, `07 QSite job control`, `07 QSite scan parsing`, `03 Column resolvers`, `03 Logging & file discovery`, `03 Journal figures`, `05 ESP charge prep`, `06 MD job build`, `06 Defluorination geometry`, `02 Mechanistic scoring`, `06 Desmond job management`, `06 MM-GBSA jobserver`, `06 Statistics & CI`, `01 Sequence merge/QC`, `04 Phylogeny/dendrogram`, `06 WaterMap & MD phases`, `03 AI-confidence figures`, `06 MM-GBSA plots`, `07 Label colour utils`, `02 Active-site mapping`?**
  _High betweenness centrality (0.561) - this node is a cross-community bridge._
- **Why does `main()` connect `02 ColabFold MSA I/O` to `CFG tier definitions`, `02 Job status & format`, `02 Analysis worker I/O`, `02 GPU reservation`, `02 MSA fetch & validate`, `02 Structure loading`, `Utils: nucleophile geometry`?**
  _High betweenness centrality (0.107) - this node is a cross-community bridge._
- **Why does `_generate_comprehensive_figures_impl()` connect `03 Journal figures` to `03 AI-confidence figures`, `03 Comparative figures`, `03 Pareto ranking`, `03 Top-tier figures`, `03 Per-model variance`, `03 Logging & file discovery`, `03 Column resolvers`, `CFG tier definitions`?**
  _High betweenness centrality (0.077) - this node is a cross-community bridge._
- **Are the 76 inferred relationships involving `CFG` (e.g. with `generate_plots()` and `calculate_sn2_metrics()`) actually correct?**
  _`CFG` has 76 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `main()` (e.g. with `print_elapsed()` and `print_script_banner()`) actually correct?**
  _`main()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _528 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `07 QM/MM console & geometry` be split into smaller, more focused modules?**
  _Cohesion score 0.049019607843137254 - nodes in this community are weakly interconnected._