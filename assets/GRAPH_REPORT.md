# Graph Report - .  (2026-08-07)

## Corpus Check
- 15 files · ~476,019 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1388 nodes · 2833 edges · 80 communities (67 shown, 13 thin omitted)
- Extraction: 91% EXTRACTED · 8% INFERRED · 0% AMBIGUOUS · INFERRED: 240 edges (avg confidence: 0.74)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Pareto ranking & conflicts|Pareto ranking & conflicts]]
- [[_COMMUNITY_MD geometry & NAC helpers|MD geometry & NAC helpers]]
- [[_COMMUNITY_Pipeline data artifacts|Pipeline data artifacts]]
- [[_COMMUNITY_CFG Boltz & QSite config|CFG Boltz & QSite config]]
- [[_COMMUNITY_03 figure rendering helpers|03 figure rendering helpers]]
- [[_COMMUNITY_Top-N PDB preparation|Top-N PDB preparation]]
- [[_COMMUNITY_06 orchestration & logging|06 orchestration & logging]]
- [[_COMMUNITY_CFG tier definitions|CFG tier definitions]]
- [[_COMMUNITY_Physics run artifacts|Physics run artifacts]]
- [[_COMMUNITY_Console log filtering|Console log filtering]]
- [[_COMMUNITY_Boltz production & ranking|Boltz production & ranking]]
- [[_COMMUNITY_06 defluorination geometry|06 defluorination geometry]]
- [[_COMMUNITY_06 progress heartbeat|06 progress heartbeat]]
- [[_COMMUNITY_03 top-tier figures|03 top-tier figures]]
- [[_COMMUNITY_Catalytic mechanism concepts|Catalytic mechanism concepts]]
- [[_COMMUNITY_07 figure style SSOT|07 figure style SSOT]]
- [[_COMMUNITY_03 figure data loading|03 figure data loading]]
- [[_COMMUNITY_06 MDWaterMap job build|06 MD/WaterMap job build]]
- [[_COMMUNITY_07 QSite job control|07 QSite job control]]
- [[_COMMUNITY_07 reactive-pose figures|07 reactive-pose figures]]
- [[_COMMUNITY_07 QSite CFG & data|07 QSite CFG & data]]
- [[_COMMUNITY_QSite scan parsing|QSite scan parsing]]
- [[_COMMUNITY_05 ESP & chain extraction|05 ESP & chain extraction]]
- [[_COMMUNITY_02 job status & justification|02 job status & justification]]
- [[_COMMUNITY_05 metadata & fig logging|05 metadata & fig logging]]
- [[_COMMUNITY_SN2 geometry extraction|SN2 geometry extraction]]
- [[_COMMUNITY_Defluorination verdict & ensemble|Defluorination verdict & ensemble]]
- [[_COMMUNITY_Mechanistic scoring engine|Mechanistic scoring engine]]
- [[_COMMUNITY_06 MM-GBSA statistics|06 MM-GBSA statistics]]
- [[_COMMUNITY_06 MD job scheduling|06 MD job scheduling]]
- [[_COMMUNITY_FASTA merge & QC|FASTA merge & QC]]
- [[_COMMUNITY_Shared utils & banners|Shared utils & banners]]
- [[_COMMUNITY_03 stats annotation helpers|03 stats annotation helpers]]
- [[_COMMUNITY_Pipeline shell driver|Pipeline shell driver]]
- [[_COMMUNITY_06 trajectory extraction|06 trajectory extraction]]
- [[_COMMUNITY_06 MM-GBSA dG series|06 MM-GBSA dG series]]
- [[_COMMUNITY_07 machinery & MMGBSA figures|07 machinery & MMGBSA figures]]
- [[_COMMUNITY_03 additional figures|03 additional figures]]
- [[_COMMUNITY_Ramachandran analysis|Ramachandran analysis]]
- [[_COMMUNITY_Interaction classification|Interaction classification]]
- [[_COMMUNITY_mmCIF to RDKit mapping|mmCIF to RDKit mapping]]
- [[_COMMUNITY_05 interaction diagrams|05 interaction diagrams]]
- [[_COMMUNITY_Geometry primitives|Geometry primitives]]
- [[_COMMUNITY_Reaction-readiness scoring|Reaction-readiness scoring]]
- [[_COMMUNITY_06 MD trajectory QC|06 MD trajectory QC]]
- [[_COMMUNITY_Atomic file IO & logging|Atomic file IO & logging]]
- [[_COMMUNITY_PBC distance utilities|PBC distance utilities]]
- [[_COMMUNITY_Console output helpers|Console output helpers]]
- [[_COMMUNITY_Conda environment mgmt|Conda environment mgmt]]
- [[_COMMUNITY_Maestro interaction detection|Maestro interaction detection]]
- [[_COMMUNITY_Dihedral & Ramachandran angles|Dihedral & Ramachandran angles]]
- [[_COMMUNITY_02 analysis worker loop|02 analysis worker loop]]
- [[_COMMUNITY_ColabFold MSA fetch|ColabFold MSA fetch]]
- [[_COMMUNITY_03 extended-analysis panels|03 extended-analysis panels]]
- [[_COMMUNITY_05 PDB atom parsing|05 PDB atom parsing]]
- [[_COMMUNITY_OOM guard|OOM guard]]
- [[_COMMUNITY_Pipeline design & controls|Pipeline design & controls]]
- [[_COMMUNITY_PLIP report parsing|PLIP report parsing]]
- [[_COMMUNITY_Free-energy landscapes|Free-energy landscapes]]
- [[_COMMUNITY_Structure superposition|Structure superposition]]
- [[_COMMUNITY_Best-pose selection|Best-pose selection]]
- [[_COMMUNITY_Boltz confidence metrics|Boltz confidence metrics]]
- [[_COMMUNITY_GPU memory management|GPU memory management]]
- [[_COMMUNITY_Active-site residue mapping|Active-site residue mapping]]
- [[_COMMUNITY_Catalytic protonation|Catalytic protonation]]
- [[_COMMUNITY_ConsoleColours drift guard|ConsoleColours drift guard]]
- [[_COMMUNITY_Boltz-2 co-folding (README)|Boltz-2 co-folding (README)]]
- [[_COMMUNITY_CFG SSOT & figures (README)|CFG SSOT & figures (README)]]
- [[_COMMUNITY_Alignment grade bins|Alignment grade bins]]
- [[_COMMUNITY_Merge QC plotting|Merge QC plotting]]
- [[_COMMUNITY_Sequence QC thresholds|Sequence QC thresholds]]
- [[_COMMUNITY_Tier standardisation|Tier standardisation]]
- [[_COMMUNITY_Ramachandran palette|Ramachandran palette]]
- [[_COMMUNITY_Git push helper|Git push helper]]
- [[_COMMUNITY_3R3U control calibration|3R3U control calibration]]
- [[_COMMUNITY_FAcD SN2 mechanism (README)|FAcD SN2 mechanism (README)]]
- [[_COMMUNITY_Alpha-carbon gating (README)|Alpha-carbon gating (README)]]
- [[_COMMUNITY_Step 01 Merge (README)|Step 01 Merge (README)]]
- [[_COMMUNITY_Step 04 Dendrogram (README)|Step 04 Dendrogram (README)]]
- [[_COMMUNITY_Console rule filter|Console rule filter]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 90 edges
2. `main()` - 36 edges
3. `main()` - 34 edges
4. `process_single_job()` - 34 edges
5. `_generate_comprehensive_figures_impl()` - 32 edges
6. `_echo()` - 31 edges
7. `console_info()` - 31 edges
8. `prep_and_convert_phase()` - 25 edges
9. `deflx_fig_name()` - 22 edges
10. `run_mmgbsa_phase()` - 22 edges

## Surprising Connections (you probably didn't know these)
- `requirements.txt (pip freeze)` --conceptually_related_to--> `02_Production_DeFluorX.py`  [INFERRED]
  requirements.txt → 00_00_run_pipeline_DeFluorX.sh
- `main (Boltz-2 orchestrator)` --conceptually_related_to--> `PFAS conda environment`  [INFERRED]
  02_Production_DeFluorX.py → PFAS.yml
- `process_and_write` --conceptually_related_to--> `biopython==1.84 (SeqIO / PairwiseAligner)`  [INFERRED]
  01_Merge_DeFluorX.py → PFAS.yml
- `boltz predict subprocess` --references--> `boltz==2.2.1 (dependency)`  [INFERRED]
  02_Production_DeFluorX.py → PFAS.yml
- `map_active_site_residues` --conceptually_related_to--> `biopython==1.84 (SeqIO / PairwiseAligner)`  [INFERRED]
  02_Production_DeFluorX.py → PFAS.yml

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **CFG single source of truth consumed pipeline-wide** — cfg_dataclass, step02_production, step03_validation_figures, step06_physics, step07_md_qmmm, utils_apply_figure_style [INFERRED 0.85]
- **Figure typography/format SSOT (CFG + utils router)** — cfg_VIS_FONT, cfg_VIS_FIGURE_FORMAT, utils_apply_figure_style, utils_deflx_fig_name [INFERRED 0.85]
- **Schrodinger job-server + OOMD lifecycle for steps 05-07** — run_pipeline_ensure_jobserver, run_pipeline_cancel_schrodinger, run_pipeline_oomd_mask, step05_topn_prep, step06_physics, step07_md_qmmm [INFERRED 0.80]
- **Boltz-2 co-fold requires PFAS conda runtime** — prod_run_boltz_predict, env_boltz, env_PFAS [INFERRED 0.80]
- **shared Biopython sequence handling across merge and production** — merge_is_valid_protein, prod_map_active_site_residues, env_biopython [INFERRED 0.60]
- **Step 06 physics per-rank pipeline (WaterMap->Build->MD->MM-GBSA->Defluorination)** — phys06_run_watermap, phys06_run_build, phys06_run_md, phys06_run_mmgbsa, phys06_run_defluorination [INFERRED 0.85]
- **CFG-SSOT no-title SVG figures (palette+VIS_FONT_* from CFG)** — prep05_plot_pose_drift, phys06_make_md_qc_figure, phys06_plot_mmgbsa_combined, phys06_plot_watermap_combined [INFERRED 0.70]
- **Geometric C-F cleavage verdict flow** — qmmm_qsite_scan_geometries, qmmm_extract_cf_distance_series, qmmm_qsite_landmarks, qmmm_qsite_val_at, qmmm_qsite_is_cleaved_geom, qmmm_defluor_verdict, cfg_QSITE_CF_CLEAVED_A [EXTRACTED 0.90]
- **SVG no-title figure builders** — qmmm_plot_qsite_reaction_profile, qmmm_plot_qsite_ensemble_profiles, qmmm_plot_qsite_profiles_all_jobs, qmmm_plot_mmgbsa_decomposition, qmmm_plot_machinery_engagement, qmmm_plot_md_qsite_timeline, qmmm_generate_global_comparative_dashboard, qmmm_generate_viability_bar_chart, qmmm_generate_comparative_residue_engagement, qmmm_generate_defluorination_landscape [EXTRACTED 0.90]
- **CFG label/colour SSOT consumers** — qmmm_format_job_label, qmmm_format_job_label_short, qmmm_CTRL_COLOUR, cfg_VIS_LIGAND_SHORT, cfg_DEFLUOR_JOB_PALETTE, cfg_VIS_ACCENT_control [INFERRED 0.85]
- **DeFluorX 8-step pipeline flow (00_01 Config -> 07 QM/MM)** — readme_pipeline_00_config, readme_pipeline_01_merge, readme_pipeline_02_production, readme_pipeline_03_figures, readme_pipeline_04_dendrogram, readme_pipeline_05_topn_prep, readme_pipeline_06_physics, readme_pipeline_07_qmmm [EXTRACTED 1.00]
- **Tier gating axes (mechanistic score + angle + competence + BDE ceiling + criteria A/B)** — readme_tier_system, readme_mechanistic_score_effective, readme_sn2_angle_gate, readme_competence_score, readme_cf_bde_ceiling, readme_criterion_a, readme_criterion_b [EXTRACTED 0.90]
- **Eight mapped catalytic residues (triad + clamp + fluoride cradle)** — readme_catalytic_triad, readme_asp110_nucleophile, readme_carboxylate_clamp, readme_fluoride_cradle [EXTRACTED 1.00]

## Communities (80 total, 13 thin omitted)

### Community 0 - "Pareto ranking & conflicts"
Cohesion: 0.05
Nodes (62): analyse_conflicts(), calculate_pareto_fronts(), _fig24_sankey(), _fig25_pfas_size(), _fig_folder07_pfas(), _is_control_mask(), _md_ready_df(), perform_advanced_ranking() (+54 more)

### Community 1 - "MD geometry & NAC helpers"
Cohesion: 0.05
Nodes (52): _blockade_vec(), calculate_min_distance(), check_md_equilibration(), compute_wm_csv_stats(), console_info(), _eaf_at(), extract_8residue_indices(), extract_hybrid_smart_system() (+44 more)

### Community 2 - "Pipeline data artifacts"
Cohesion: 0.05
Nodes (47): Active_Site_Alignments stats (FILE_ALIGNMENT_STATS), 1_Boltz2_Production input FASTA panel, 5_Boltz2_DeFluorX_Master_*.csv, C_INP_Merged_for_Boltz-2.fasta, 03_Figure_Enriched_Dataset.csv (FILE_VALIDATED_MASTER), CFG.COL_TIER, CFG.DENDRO_DISTANCE_METRIC, CFG.DENDRO_KMER_SIZE (+39 more)

### Community 3 - "CFG Boltz & QSite config"
Cohesion: 0.05
Nodes (46): CFG BOLTZ_* prediction/MSA params, CFG.BOLTZ_EXECUTABLE / BOLTZ_MODEL_VERSION, CFG control ligand panel + expected verdicts, CFG QSITE_* QM/MM extraction settings, CFG REF_ACTIVE_SITE_MAP (3R3U residues), CFG SMART_LOCK_* nucleophile-search bias, CFG interaction thresholds (HB/salt/pi/halogen), CFG THRESHOLD_* interaction geometry (Maestro) (+38 more)

### Community 4 - "03 figure rendering helpers"
Cohesion: 0.10
Nodes (41): _auto_label_colour(), _generate_comprehensive_figures_impl(), _jf_box(), _jf_circos(), _jf_fit_fs(), _jf_importance(), _jf_ligshort(), _jf_manhattan() (+33 more)

### Community 5 - "Top-N PDB preparation"
Cohesion: 0.08
Nodes (42): check_prep_needed(), _check_residue_identity_guard(), cif_to_pdb_gemmi(), collect_best_cifs(), console_info(), generate_esp_charges(), generate_raw_step(), get_source_tag() (+34 more)

### Community 6 - "06 orchestration & logging"
Cohesion: 0.10
Nodes (38): _analysis_dir(), _complex_label(), _echo(), _emit_timings(), _fail(), _fmt_dur(), _log(), main() (+30 more)

### Community 7 - "CFG tier definitions"
Cohesion: 0.07
Nodes (26): CFG, ===============================================================================, Central Configuration Repository - DeFluorX Pipeline.      Attribute prefix → se, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, _avail_ram_gb(), _await_mmgbsa_jobserver(), _cleanup_mmgbsa_shards(), _diagnose_mmgbsa_failure() (+18 more)

### Community 8 - "Physics run artifacts"
Cohesion: 0.08
Nodes (36): 2_Best_Complexes_CIFs/*.cif, R_N_<stem>_ESP_Complex.mae, <stem>_ESP.mae ligand ESP charges, 05_MD_Simulations desmond -out.cms + trj + mmgbsa csv, 06_Analysis/00_MMGBSA_Summary.csv, R{N}_<stem>.pdb prepared handover, 6_Boltz2_DeFluorX_Ranked_*.csv, 03_WaterMaps/watermap_R_N *_wm.maegz + csv (+28 more)

### Community 9 - "Console log filtering"
Cohesion: 0.09
Nodes (25): _ConsoleRuleFilter, install_console_rule_filter(), A stdout wrapper that collapses consecutive separator rules.      The logs grow, Collapse stacked separator rules for the rest of this process's output., Simultaneous console + file logger shared by the pipeline steps.      The log pa, Record a line in the log file WITHOUT printing it to the console.          For o, ReportManager, clean_id() (+17 more)

### Community 10 - "Boltz production & ranking"
Cohesion: 0.10
Nodes (35): atomic_to_csv(), colabfold_a3m_path(), colabfold_meta_path(), console_info(), console_separator(), extract_3r3u_sequence(), generate_scientific_ranking_csv(), load_cached_alignments() (+27 more)

### Community 11 - "06 defluorination geometry"
Cohesion: 0.08
Nodes (34): _cancel_launched_jobs(), _clean_job_env(), _defl_ligand_atoms(), _defl_mapped_residues(), _defl_min_image(), _defl_resid(), find_esp(), _hydrate() (+26 more)

### Community 12 - "06 progress heartbeat"
Cohesion: 0.08
Nodes (17): _close_bar(), Heartbeat, MDHeartbeat, _progress_line(), (current relaxation stage, production stage number) for the pre-production phase, Refresh the single shared \r progress line (terminal only, never the log)., Close the open \r progress line with a newline (phase change / completion)., Periodically report a long Schrödinger step's progress by tailing its log. (+9 more)

### Community 13 - "03 top-tier figures"
Cohesion: 0.10
Nodes (34): _fig23_multitarget(), _fig23b_toptier_breakdown(), _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic(), _fig_14b_tt_interactions(), _fig_folder04_ai_confidence(), _fig_folder06_ligand(), _fig_path() (+26 more)

### Community 14 - "Catalytic mechanism concepts"
Cohesion: 0.06
Nodes (34): Asp110 nucleophile (backside attack on Ca), Bidentate carboxylate clamp (Arg111/Arg114), Catalytic triad Asp110 (nuc) / Asp134 (acid) / His280 (base), Elite C-F BDE ceiling (TIER_ELITE_BDE_MAX 128 kcal/mol), competence_score substrate-feasibility floor, Criterion A - active-site integrity (8 residues), Criterion B - catalytic constellation RMSD, Crystal reference PDB 3R3U (R. palustris FAcD, 1.60 A) (+26 more)

### Community 15 - "07 figure style SSOT"
Cohesion: 0.09
Nodes (33): apply_figure_style(), clean_spines(), deflx_fig_name(), Return a figure path with its suffix swapped to the ACTIVE SSOT format (CFG.VIS_, The pipeline's one typography and canvas definition, applied to matplotlib's rcP, Academic-style axes: remove top/right spines, thin the remaining borders.      C, _active_site_series(), format_job_label_short() (+25 more)

### Community 16 - "03 figure data loading"
Cohesion: 0.09
Nodes (33): latest_by_mtime(), The newest existing file by modification time, or None.      Selection is by mti, _aux_dir(), console_info(), generate_comprehensive_figures(), generate_ramachandran_figures(), load_and_prep_data(), _load_module() (+25 more)

### Community 17 - "06 MD/WaterMap job build"
Cohesion: 0.07
Nodes (33): _build_msj(), _csv_nonempty(), _diagnose_shard_failure(), discover_handover(), export_watermap_csv(), _lig_atom_indices(), _load_module(), _md_msj() (+25 more)

### Community 18 - "07 QSite job control"
Cohesion: 0.09
Nodes (31): console_qmm_ready(), console_separator(), console_title(), format_job_label(), _kill_all_qsite_jobs(), _kill_qsite_job(), load_triad_mapping(), main() (+23 more)

### Community 19 - "07 reactive-pose figures"
Cohesion: 0.09
Nodes (32): _draw_reactive_pose_for_rank(), _draw_trajectory_figures(), generate_reactive_pose_figures(), _load_module(), _load_ranked_df(), _load_reactive_pose_data(), _mae_av(), _mae_parse_cts() (+24 more)

### Community 20 - "07 QSite CFG & data"
Cohesion: 0.12
Nodes (32): CFG.DEFLUOR_JOB_PALETTE, CFG.QSITE_CF_CLEAVED_A, CFG.VIS_ACCENT['control'], CFG.VIS_LIGAND_SHORT, 01_MD_Stats.json, 02_NAC_Data.csv, 11_QSite_Scan_Data.csv, _CTRL_COLOUR (+24 more)

### Community 21 - "QSite scan parsing"
Cohesion: 0.10
Nodes (30): _extract_cf_distance_series(), _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), _extract_scan_points(), _parse_qsite_pes_raw(), parse_qsite_profile(), _qsite_frame_summary() (+22 more)

### Community 22 - "05 ESP & chain extraction"
Cohesion: 0.10
Nodes (26): console_separator(), _esp_run(), extract_chain_l_mol(), _is_reference_control(), _lig_short_token(), _load_module(), load_reference_data(), main() (+18 more)

### Community 23 - "02 job status & justification"
Cohesion: 0.10
Nodes (24): check_job_status(), cpu_usage_summary(), extract_short_fasta_id(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), heal_smiles_in_files() (+16 more)

### Community 24 - "05 metadata & fig logging"
Cohesion: 0.12
Nodes (19): _append_auxiliary_log(), _fig_constants(), _fig_log(), load_metadata(), Return a dict of figure-generation constants from CFG., Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir., Unified logging for figure generation - routes through console_info., Reads a specific tool's output log and appends it to the main report. (+11 more)

### Community 25 - "SN2 geometry extraction"
Cohesion: 0.09
Nodes (24): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+16 more)

### Community 26 - "Defluorination verdict & ensemble"
Cohesion: 0.10
Nodes (23): _collect_qsite_results(), defluor_propensity(), defluor_verdict(), parse_qsite_barrier(), plot_qsite_ensemble_profiles(), plot_qsite_reaction_profile(), _qsite_ensemble_barrier(), _qsite_ensemble_derxn() (+15 more)

### Community 27 - "Mechanistic scoring engine"
Cohesion: 0.11
Nodes (18): Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge() (+10 more)

### Community 28 - "06 MM-GBSA statistics"
Cohesion: 0.10
Nodes (22): _avg_dg(), _block_bootstrap_median_ci(), _boltzmann_mean_dg(), _cliffs_delta(), _ctrl_palette(), _delta_word(), _dg_failures(), _draw_time_cumulative() (+14 more)

### Community 29 - "06 MD job scheduling"
Cohesion: 0.13
Nodes (20): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Classify every desmond_md_job_R_* directory and return the subset     that still, Step 1: event_analysis.py - generates the SID-in.eaf descriptor., Step 2: analyze_simulation.py - generates the SID-out.eaf result vector.      Us (+12 more)

### Community 30 - "FASTA merge & QC"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 31 - "Shared utils & banners"
Cohesion: 0.11
Nodes (17): calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner(), ===============================================================================, Geometry-based Oδ fallback nucleophile search for the Smart-Lock algorithm. (+9 more)

### Community 32 - "03 stats annotation helpers"
Cohesion: 0.16
Nodes (18): Return the canonical tier ordering restricted to tiers actually present., Compact, publication-style p-value formatting., Place a boxed statistics annotation on an axis., One combined box: the legend entries, then the statistics lines as blank-handle, Persist a figure at publication resolution and free its memory., Binding_Affinity_Score violin + Affinity/Pocket-ratio/Density mean ± 95% CI line, Vertical dotted separators between the n categorical groups on a per-tier x-axis, _tier_seps() (+10 more)

### Community 33 - "Pipeline shell driver"
Cohesion: 0.22
Nodes (13): _cancel_schrodinger_jobs(), _ensure_jobserver(), _print_step_menu(), _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_DeFluorX.sh script, _sync_staging_log() (+5 more)

### Community 34 - "06 trajectory extraction"
Cohesion: 0.14
Nodes (11): _eta_str(), _extract_with_progress(), Refresh the shared single \r progress line from a CPU worker (extraction / SID /, Extract members from tgz into wd (dropping one leading path component), drawing, Put the production trajectory at the job-dir root as {job}_trj + {job}.ene and r, A ' · ETA <dur>' suffix from a linear extrapolation of the current rate - empty, One in-place progress line for a sharded MM-GBSA run.      Aggregates across the, (frames read across all shards, shards currently in Prime minimisation). (+3 more)

### Community 35 - "06 MM-GBSA dG series"
Cohesion: 0.14
Nodes (16): _failure_windows(), _lookup_controls(), _lookup_ligands(), _lookup_tiers(), _mmgbsa_dg_series(), _plot_rank_mmgbsa(), Extract the per-frame ΔG_bind series, tolerant of column-name variants., Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier     co (+8 more)

### Community 36 - "07 machinery & MMGBSA figures"
Cohesion: 0.18
Nodes (14): auto_label_colour(), The contrast colour for a label written ON a filled mark, from the fill's lumina, _darken(), _engage_zones(), _master_container(), plot_machinery_engagement(), plot_mmgbsa_decomposition(), The criteria that define engagement - every one of them a CFG constant. (+6 more)

### Community 37 - "03 additional figures"
Cohesion: 0.21
Nodes (14): _control_star_df(), _ctrl_star(), _diag10_model_agreement(), _fig_folder03_dataset(), _fig_folder05_catalytic(), generate_additional_figures(), _md_ready_stars_cat(), Folder 03_Dataset_and_Alignment_Overview - coverage, tiers, alignment grades. (+6 more)

### Community 38 - "Ramachandran analysis"
Cohesion: 0.22
Nodes (13): _draw_rama_background(), Any, _rama_classify(), _rama_palette(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favoured, allowed, and outlier regions., Ramachandran colours from CFG.VIS_RAMA (the single source of truth); the built-i (+5 more)

### Community 39 - "Interaction classification"
Cohesion: 0.20
Nodes (12): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi, Classify a ligand atom as 'anion'- or 'cation'-capable for salt-bridge     detec (+4 more)

### Community 40 - "mmCIF to RDKit mapping"
Cohesion: 0.18
Nodes (11): graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit(), Constructs a three-dimensional RDKit molecule from a SMILES string.     This is, Extracts raw XYZ coordination matrices from RDKit molecule data blocks., Determines the optimal rotation/translation matrix required to superimpose spati, Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not b, FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perce (+3 more)

### Community 41 - "05 interaction diagrams"
Cohesion: 0.18
Nodes (11): _draw_interaction_diagram(), _im_project(), _im_render_diagram(), _im_resnum(), _im_role_resnums(), _im_separate_atoms(), Shared 2D interaction renderer for InteractionMap (mode='distance') and PLIP (mo, {job_name: {residue_number: role_key}} from the ranked CSV Mapped_* columns, so (+3 more)

### Community 42 - "Geometry primitives"
Cohesion: 0.20
Nodes (10): calculate_angle(), calculate_burgi_dunitz(), calculate_improper_dihedral(), distance(), get_mic_vector(), Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for, MIC displacement vector from pos2 → pos1, supporting orthorhombic and     tricli, PBC-corrected or Euclidean distance between two Cartesian coordinates (Å). (+2 more)

### Community 43 - "Reaction-readiness scoring"
Cohesion: 0.20
Nodes (10): Distance (Å) -> 0-1 proximity; NaN (sentinel / not engaged) -> 0., Occupancy adequacy tent: too empty and overflow both score low; a substrate-size, Per-complex sub-scores + geometry_readiness, feasibility (=02 competence_score), Pocket + 8-residue reaction-readiness per ligand, ordered by molecular size (4 p, _rr_family(), _rr_occ(), _rr_prox(), _rr_scores() (+2 more)

### Community 44 - "06 MD trajectory QC"
Cohesion: 0.20
Nodes (10): _ctrl_label(), make_md_qc_figure(), _natural_rank(), The ligand's short name (FA / DFA / TFA) from CFG - the same abbreviations every, 3R3U-FA' for a control rank, else the short PFAS name (fallback R_<rk>)., Read one finished MD job's SID .eaf + Desmond .ene into the trajectory-QC arrays, Draw the MD trajectory-QC figure - Cα-RMSD, ligand RMSD, temperature, Cα-RMSF -, Sort key for desmond_md_job_R_N directories (numeric, -V style).      Matches th (+2 more)

### Community 45 - "Atomic file IO & logging"
Cohesion: 0.22
Nodes (9): atomic_write_csv(), compute_ligand_properties(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, Canonical file logger for all pipeline scripts.      Creates a file-only handler, Write a JSON file so a reader never sees a half-written one.      The write goes, Write a DataFrame to CSV so a reader never sees a half-written file.      Same g, setup_logging() (+1 more)

### Community 46 - "PBC distance utilities"
Cohesion: 0.31
Nodes (9): _box_diag(), calculate_min_distance(), _ensure_box_3x3(), mic_dists_2d(), ndarray, Minimum PBC-corrected distance between two sets of Cartesian positions., Coerce a periodic box to a (3,3) float64 matrix.     Handles: (3,3) arrays, flat, Orthorhombic fallback box lengths as a safe (3,) diagonal.      Used when the fu (+1 more)

### Community 47 - "Console output helpers"
Cohesion: 0.25
Nodes (9): console_info(), console_separator(), console_title(), Logger, Remove all ANSI/VT100 escape sequences from a string., Bold section header - printed and optionally written to log file., Two-space-indented info line - printed and optionally written to log file., Horizontal rule - printed and optionally written to log file.      Parameters (+1 more)

### Community 48 - "Conda environment mgmt"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 49 - "Maestro interaction detection"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 50 - "Dihedral & Ramachandran angles"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 51 - "02 analysis worker loop"
Cohesion: 0.33
Nodes (6): analysis_worker_loop(), append_rows_to_csv(), flatten_job_result(), Flattens a process_single_job result dict into a single-level row for CSV writin, Appends a list of result rows to the master CSV, deduplicates, sorts by job_name, Background worker that analyses freshly GPU-predicted jobs and appends rows to t

### Community 52 - "ColabFold MSA fetch"
Cohesion: 0.33
Nodes (6): fetch_msa_direct(), Generates a stable hash for FASTA sequences to verify consistency during resumed, Validates the A3M MSA file to detect corruption. Returns True if the file struct, Fetches the MSA directly from the ColabFold REST API. This circumvents the Boltz, sequence_hash(), validate_a3m_file()

### Community 53 - "03 extended-analysis panels"
Cohesion: 0.33
Nodes (6): _ext_match_03_style(), Create Tier_1A enzyme × ligand heatmap for lab validation targets., Bring an extended-analysis panel onto 03's typography before it is written., Persist a figure at the pipeline's publication resolution, then free its memory., _xn_figure_06a(), _xn__save()

### Community 54 - "05 PDB atom parsing"
Cohesion: 0.33
Nodes (6): _im_don_acc(), _im_parse_pdb(), _im_ring(), Split heavy N/O atoms into donors (with bonded H) and acceptors.      Returns (d, Aromatic ring (centroid, normal) for an aromatic residue, else None., Return (lig_atoms, pro_residues, lig_hb) from a prepared (protonated) PDB.

### Community 56 - "Pipeline design & controls"
Cohesion: 0.33
Nodes (6): Control ligands FA/DFA (positive) and TFA (decoy), DeFluorX pipeline (FAcD/PFAS defluorination discovery), Four-phase pipeline design (Foundation/Ingest/Screening/Dynamics-QM), 27 PFAS ligand screening panel, PFAS persistence problem (stable C-F bond, ~544 kJ/mol BDE), WHY FA/DFA positive controls + TFA decoy (threshold sanity check before trusting screen)

### Community 57 - "PLIP report parsing"
Cohesion: 0.40
Nodes (5): _parse_plip_xml(), _plip_cfg_cutoffs(), _plip_coo(), Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]., Parse a PLIP XML report into a contacts list compatible with _im_project /     _

### Community 58 - "Free-energy landscapes"
Cohesion: 0.40
Nodes (5): _fel_grid(), _fel_surface(), generate_free_energy_landscapes(), 2D histogram → free energy G = −kT ln P (kcal/mol), minimum shifted to zero., Two free-energy surfaces in one figure, on one shared discrete energy scale.

### Community 59 - "Structure superposition"
Cohesion: 0.50
Nodes (4): analyse_candidate_structure(), load_structure_safe(), Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilit, Superimposes the candidate structure onto the DeHa4 control model.     Computes

### Community 60 - "Best-pose selection"
Cohesion: 0.50
Nodes (4): analyse_model_task(), Global worker routine functionally decoupled for robust serialisation capabiliti, Selects the representative pose across the Boltz diffusion samples tier-first: t, select_best_degrader_model()

### Community 61 - "Boltz confidence metrics"
Cohesion: 0.50
Nodes (4): compute_cross_interface_pae(), load_extra_boltz_metrics(), Calculates the Predicted Aligned Error (PAE), focusing specifically on the prote, Parses background statistical validation tensors directly from Boltz NPZ file du

### Community 62 - "GPU memory management"
Cohesion: 0.50
Nodes (4): init_gpu_reservations(), query_gpu_memory(), A resilient system query that safely falls back to CPU memory context if nvidia-, Maps unallocated VRAM memory segments natively upon system initialisation.

### Community 63 - "Active-site residue mapping"
Cohesion: 0.50
Nodes (4): map_active_site_residues(), Append a newly computed alignment to the on-disk cache file, under an exclusive, Performs a global sequence alignment between the canonical reference sequence, update_alignment_stats()

### Community 64 - "Catalytic protonation"
Cohesion: 0.50
Nodes (4): _build_hd1_line(), enforce_catalytic_protonation(), Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev

### Community 65 - "ConsoleColours drift guard"
Cohesion: 0.67
Nodes (3): 00_03 local ConsoleColours copy, 00_03 ConsoleColours drift assertion, 00_02 ConsoleColours (canonical)

### Community 66 - "Boltz-2 co-folding (README)"
Cohesion: 0.67
Nodes (3): Boltz-2 diffusion co-folding structure prediction (step 02), ColabFold MSA evolutionary context, Step 02 Production (Boltz-2 co-fold + scoring engine)

### Community 67 - "CFG SSOT & figures (README)"
Cohesion: 0.67
Nodes (3): Step 00_01 Central Configuration (CFG SSOT), Step 03 Validation Figures, CFG.VIS_FIGURE_FORMAT single-switch figure output (SVG SSOT)

## Ambiguous Edges - Review These
- `run_prepwizard (Schrodinger PrepWizard)` → `OPLS_2005 PrepWizard vs OPLS4 downstream mismatch`  [AMBIGUOUS]
  05_TopN_and_PDB_Preparation_DeFluorX.py · relation: conceptually_related_to

## Knowledge Gaps
- **98 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `git_push_DeFluorX.sh script` (+93 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **13 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `run_prepwizard (Schrodinger PrepWizard)` and `OPLS_2005 PrepWizard vs OPLS4 downstream mismatch`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `CFG` connect `CFG tier definitions` to `Pareto ranking & conflicts`, `MD geometry & NAC helpers`, `03 figure rendering helpers`, `Top-N PDB preparation`, `06 orchestration & logging`, `Console log filtering`, `Boltz production & ranking`, `06 defluorination geometry`, `03 top-tier figures`, `07 figure style SSOT`, `03 figure data loading`, `06 MD/WaterMap job build`, `07 QSite job control`, `07 reactive-pose figures`, `QSite scan parsing`, `05 ESP & chain extraction`, `Defluorination verdict & ensemble`, `Mechanistic scoring engine`, `06 MM-GBSA statistics`, `06 MD job scheduling`, `FASTA merge & QC`, `06 MM-GBSA dG series`, `07 machinery & MMGBSA figures`, `06 MD trajectory QC`, `03 extended-analysis panels`, `Free-energy landscapes`, `Active-site residue mapping`?**
  _High betweenness centrality (0.420) - this node is a cross-community bridge._
- **Why does `main()` connect `Boltz production & ranking` to `CFG tier definitions`, `02 analysis worker loop`, `ColabFold MSA fetch`, `02 job status & justification`, `Structure superposition`, `GPU memory management`, `Shared utils & banners`?**
  _High betweenness centrality (0.069) - this node is a cross-community bridge._
- **Why does `_generate_comprehensive_figures_impl()` connect `03 figure rendering helpers` to `Pareto ranking & conflicts`, `03 stats annotation helpers`, `03 additional figures`, `CFG tier definitions`, `Console log filtering`, `Reaction-readiness scoring`, `03 top-tier figures`, `03 figure data loading`, `03 extended-analysis panels`?**
  _High betweenness centrality (0.062) - this node is a cross-community bridge._
- **Are the 75 inferred relationships involving `CFG` (e.g. with `generate_plots()` and `calculate_sn2_metrics()`) actually correct?**
  _`CFG` has 75 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `main()` (e.g. with `print_elapsed()` and `print_script_banner()`) actually correct?**
  _`main()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _592 weakly-connected nodes found - possible documentation gaps or missing edges._