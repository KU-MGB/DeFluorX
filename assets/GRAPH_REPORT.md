# Graph Report - /tmp/claude-1001/-mnt-wdpassport-02-DeFluorX-02-Data/ca98d275-4986-4c47-ab4c-ddf45e6999c4/scratchpad/dfx_src  (2026-10-19)

## Corpus Check
- 12 files · ~269,193 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1250 nodes · 2733 edges · 62 communities
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 177 edges (avg confidence: 0.73)
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
1. `CFG` - 90 edges
2. `main()` - 36 edges
3. `process_single_job()` - 35 edges
4. `main()` - 34 edges
5. `_generate_comprehensive_figures_impl()` - 32 edges
6. `_echo()` - 31 edges
7. `console_info()` - 31 edges
8. `prep_and_convert_phase()` - 25 edges
9. `run_mmgbsa_phase()` - 22 edges
10. `deflx_fig_name()` - 21 edges

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

## Communities (62 total, 0 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.05
Nodes (52): _aggregate_mmgbsa_stats(), _blockade_vec(), calculate_min_distance(), check_md_equilibration(), compute_wm_csv_stats(), console_info(), _eaf_at(), extract_8residue_indices() (+44 more)

### Community 1 - "Community 1"
Cohesion: 0.07
Nodes (52): _auto_label_colour(), _fig23b_toptier_breakdown(), _generate_comprehensive_figures_impl(), _jf_box(), _jf_circos(), _jf_fit_fs(), _jf_importance(), _jf_ligshort() (+44 more)

### Community 2 - "Community 2"
Cohesion: 0.08
Nodes (48): _control_star_df(), _ctrl_star(), _diag10_model_agreement(), _fig25_pfas_size(), _fig_05b_tt_ai_quality(), _fig_13b_tt_mechanistic(), _fig_14b_tt_interactions(), _fig_folder03_dataset() (+40 more)

### Community 3 - "Community 3"
Cohesion: 0.09
Nodes (44): _analysis_dir(), _complex_label(), _ctrl_label(), _ctrl_palette(), _echo(), _emit_timings(), _fail(), _log() (+36 more)

### Community 4 - "Community 4"
Cohesion: 0.08
Nodes (41): analyse_conflicts(), calculate_pareto_fronts(), _fig23_multitarget(), _fig24_sankey(), _fig_folder07_pfas(), _is_control_mask(), perform_advanced_ranking(), DataFrame (+33 more)

### Community 5 - "Community 5"
Cohesion: 0.08
Nodes (40): apply_figure_style(), clean_spines(), deflx_fig_name(), Return a figure path with its suffix swapped to the ACTIVE SSOT format (CFG.VIS_, The pipeline's one typography and canvas definition, applied to matplotlib's rcP, Academic-style axes: remove top/right spines, thin the remaining borders.      C, _active_site_series(), _draw_trajectory_figures() (+32 more)

### Community 6 - "Community 6"
Cohesion: 0.07
Nodes (28): CFG, ===============================================================================, Central Configuration Repository - DeFluorX Pipeline.      Attribute prefix → se, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, _avg_dg(), _cleanup_mmgbsa_shards(), _dg_failures(), _mmgbsa_complete() (+20 more)

### Community 7 - "Community 7"
Cohesion: 0.09
Nodes (25): _ConsoleRuleFilter, install_console_rule_filter(), A stdout wrapper that collapses consecutive separator rules.      The logs grow, Collapse stacked separator rules for the rest of this process's output., Simultaneous console + file logger shared by the pipeline steps.      The log pa, Record a line in the log file WITHOUT printing it to the console.          For o, ReportManager, clean_id() (+17 more)

### Community 8 - "Community 8"
Cohesion: 0.10
Nodes (35): atomic_to_csv(), colabfold_a3m_path(), colabfold_meta_path(), console_info(), console_separator(), extract_3r3u_sequence(), generate_scientific_ranking_csv(), load_cached_alignments() (+27 more)

### Community 9 - "Community 9"
Cohesion: 0.08
Nodes (34): _check_residue_identity_guard(), collect_best_cifs(), console_separator(), extract_chain_l_mol(), index_existing_files(), load_catalytic_anchor_map(), load_machinery_map(), load_md_selected_jobs() (+26 more)

### Community 10 - "Community 10"
Cohesion: 0.09
Nodes (33): console_qmm_ready(), console_separator(), console_title(), format_job_label(), _kill_all_qsite_jobs(), _kill_qsite_job(), _latest_schrodinger(), load_triad_mapping() (+25 more)

### Community 11 - "Community 11"
Cohesion: 0.07
Nodes (33): _build_msj(), _csv_nonempty(), _diagnose_mmgbsa_failure(), _diagnose_shard_failure(), discover_handover(), export_watermap_csv(), _lig_atom_indices(), _load_module() (+25 more)

### Community 12 - "Community 12"
Cohesion: 0.08
Nodes (15): _close_bar(), Heartbeat, MDHeartbeat, True once multisim has entered the Production stage (and not yet completed it)., (current relaxation stage, production stage number) for the pre-production phase, Close the open \r progress line with a newline (phase change / completion)., Periodically report a long Schrödinger step's progress by tailing its log., Report progress as a full logged line (runs on the CPU worker, concurrent with t (+7 more)

### Community 13 - "Community 13"
Cohesion: 0.08
Nodes (32): _ext_match_03_style(), Create Tier_1A enzyme × ligand heatmap for lab validation targets., Compact, publication-style p-value formatting., Bring an extended-analysis panel onto 03's typography before it is written., Persist a figure at the pipeline's publication resolution, then free its memory., No-op: panel letters are not drawn.      No other figure in this set carries an, A sample size that fits the ~0.5 in a tier occupies on a two-panel figure., The panel's test result, on its own line ABOVE the axes.      Inside the axes it (+24 more)

### Community 14 - "Community 14"
Cohesion: 0.09
Nodes (32): _extract_cf_distance_series(), _extract_fluoride_charge_series(), _extract_scan_coordinates(), _extract_scan_energies(), _extract_scan_points(), parse_qsite_barrier(), _parse_qsite_pes_raw(), parse_qsite_profile() (+24 more)

### Community 15 - "Community 15"
Cohesion: 0.10
Nodes (30): latest_by_mtime(), The newest existing file by modification time, or None.      Selection is by mti, _aux_dir(), console_info(), generate_comprehensive_figures(), generate_ramachandran_figures(), load_and_prep_data(), _load_module() (+22 more)

### Community 16 - "Community 16"
Cohesion: 0.09
Nodes (28): _clean_job_env(), _defl_ligand_atoms(), _defl_mapped_residues(), _defl_min_image(), _defl_resid(), find_esp(), _hydrate(), _latest_schrodinger() (+20 more)

### Community 17 - "Community 17"
Cohesion: 0.11
Nodes (25): _draw_reactive_pose_for_rank(), generate_reactive_pose_figures(), _load_module(), _load_ranked_df(), _load_reactive_pose_data(), _mae_av(), _mae_parse_cts(), _mech_smiles_for() (+17 more)

### Community 18 - "Community 18"
Cohesion: 0.09
Nodes (24): An atom plus the minimum metadata needed to locate it: residue, sequence id, cha, The best SN2 geometry found in one model., The p1-p2-p3 angle in degrees.      For the SN2: p1 = the attacking Oδ of the ca, Anything that is not a standard amino acid or water is treated as the ligand., Split a Boltz-2 CIF into its ligand atoms and its protein atoms.      A parse fa, Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).      The, The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to, Rank the two carboxylate oxygens: short distance AND an angle near 180°.      Th (+16 more)

### Community 19 - "Community 19"
Cohesion: 0.10
Nodes (24): _append_auxiliary_log(), _draw_interaction_diagram(), _im_don_acc(), _im_parse_pdb(), _im_project(), _im_render_diagram(), _im_resnum(), _im_ring() (+16 more)

### Community 20 - "Community 20"
Cohesion: 0.10
Nodes (23): _collect_qsite_results(), defluor_propensity(), defluor_verdict(), plot_qsite_ensemble_profiles(), plot_qsite_reaction_profile(), _qsite_ensemble_barrier(), _qsite_ensemble_derxn(), _qsite_point_role() (+15 more)

### Community 21 - "Community 21"
Cohesion: 0.15
Nodes (22): _acc_of(), _apply_supp_layout(), build_cascade(), _c_scan_inset(), _job_class(), _key(), load_qm(), load_scan() (+14 more)

### Community 22 - "Community 22"
Cohesion: 0.11
Nodes (18): Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, calculate_sn2_metrics(), check_catalytic_geometry(), compute_pocket_fit(), _derive_burgi_dunitz(), _derive_flippin_lodge() (+10 more)

### Community 23 - "Community 23"
Cohesion: 0.13
Nodes (18): _is_reference_control(), _latest_schrodinger(), _lig_short_token(), _parse_plip_xml(), _plip_angles_ok(), _plip_cfg_cutoffs(), _plip_coo(), plot_machinery_distribution() (+10 more)

### Community 24 - "Community 24"
Cohesion: 0.17
Nodes (13): console_info(), _fig_constants(), _fig_log(), load_reference_data(), Scans the 1_Input_Data folder to map names to Sequences and SMILES.     Returns:, Return a dict of figure-generation constants from CFG., Unified logging for figure generation - routes through console_info., Handles verification and automated installation of visual tools. (+5 more)

### Community 25 - "Community 25"
Cohesion: 0.11
Nodes (15): _eta_str(), _extract_with_progress(), _fmt_dur(), _progress_line(), Refresh the single shared \r progress line (terminal only, never the log)., Refresh the shared single \r progress line from a CPU worker (extraction / SID /, Extract members from tgz into wd (dropping one leading path component), drawing, Put the production trajectory at the job-dir root as {job}_trj + {job}.ene and r (+7 more)

### Community 26 - "Community 26"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 27 - "Community 27"
Cohesion: 0.11
Nodes (17): calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner(), ===============================================================================, Geometry-based Oδ fallback nucleophile search for the Smart-Lock algorithm. (+9 more)

### Community 28 - "Community 28"
Cohesion: 0.13
Nodes (16): cpu_usage_summary(), extract_short_fasta_id(), fetch_msa_direct(), init_gpu_reservations(), _load_module(), query_gpu_memory(), Load a Python file as a module using its filesystem path, regardless of filename, A resilient system query that safely falls back to CPU memory context if nvidia- (+8 more)

### Community 29 - "Community 29"
Cohesion: 0.11
Nodes (18): check_job_status(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), heal_smiles_in_files(), mirror_best_cif(), process_single_job() (+10 more)

### Community 30 - "Community 30"
Cohesion: 0.13
Nodes (18): is_eaf_complete(), out_eaf_frames(), _proc_alive(), process_jobs(), _rank_of(), Classify every desmond_md_job_R_* directory and return the subset     that still, Step 1: event_analysis.py - generates the SID-in.eaf descriptor., Step 2: analyze_simulation.py - generates the SID-out.eaf result vector.      Us (+10 more)

### Community 31 - "Community 31"
Cohesion: 0.22
Nodes (13): _cancel_schrodinger_jobs(), _ensure_jobserver(), _print_step_menu(), _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_DeFluorX.sh script, _sync_staging_log() (+5 more)

### Community 32 - "Community 32"
Cohesion: 0.16
Nodes (17): _add_svg_picture(), _decode_data_uri(), _latest_run(), montage(), _pptx_convert(), _ptmp(), Path, _rasterise() (+9 more)

### Community 33 - "Community 33"
Cohesion: 0.12
Nodes (17): _boost_fonts(), _montage_svg(), Intrinsic (viewBox) width/height of an SVG root, so a nested panel keeps its tru, Prefix every id and every reference to it so panels never collide on shared matp, Panel aspect ratio (w/h): from the SVG viewBox for vector panels, from pixel siz, Enlarge every <text>/<tspan> font-size by factor f in place (labels are editable, One placed panel: a nested live <svg> for a vector source, an <image> for a rast, Vector panel letter: a bold near-black glyph with a thin white halo for legibili (+9 more)

### Community 34 - "Community 34"
Cohesion: 0.17
Nodes (16): Return the canonical tier ordering restricted to tiers actually present., Compact, publication-style p-value formatting., Place a boxed statistics annotation on an axis., One combined box: the legend entries, then the statistics lines as blank-handle, Persist a figure at publication resolution and free its memory., Binding_Affinity_Score violin + Affinity/Pocket-ratio/Density mean ± 95% CI line, _xo__annotate(), _xo__fig_02A_binding_affinity_metrics() (+8 more)

### Community 35 - "Community 35"
Cohesion: 0.13
Nodes (16): check_prep_needed(), cif_to_pdb_gemmi(), generate_raw_step(), get_source_tag(), measure_sn2_geometry(), preparation_step(), Measure the SN2 geometry of a structure - CIF or PDB - with the pipeline's own j, Task 2: Protein Preparation.     Logic: Checks consistency between Raw and Prep (+8 more)

### Community 36 - "Community 36"
Cohesion: 0.17
Nodes (15): _box_diag(), calculate_improper_dihedral(), calculate_min_distance(), distance(), _ensure_box_3x3(), get_mic_vector(), mic_dists_2d(), ndarray (+7 more)

### Community 37 - "Community 37"
Cohesion: 0.18
Nodes (14): auto_label_colour(), The contrast colour for a label written ON a filled mark, from the fill's lumina, _darken(), _engage_zones(), plot_machinery_engagement(), plot_mmgbsa_decomposition(), _primary_container(), The criteria that define engagement - every one of them a CFG constant. (+6 more)

### Community 38 - "Community 38"
Cohesion: 0.15
Nodes (14): _avail_ram_gb(), _await_mmgbsa_jobserver(), _free_swap_gb(), Free RAM the kernel expects to hand out without swapping (MemAvailable)., Free swap (SwapFree). Only added to the MM-GBSA budget when MMGBSA_RAM_SWAP_FRAC, Split [0, total) into contiguous frame shards and decide how many run at once., Record which trajectory frame each Prime row actually came from.      thermal_mm, Score the whole trajectory as concurrent contiguous frame shards.      Each shar (+6 more)

### Community 39 - "Community 39"
Cohesion: 0.14
Nodes (14): _boltzmann_mean_dg(), _failure_windows(), _lookup_controls(), _lookup_ligands(), _lookup_tiers(), _natural_rank(), Log-sum-exp ("Boltzmann") ensemble mean binding free energy over frames., Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier     co (+6 more)

### Community 40 - "Community 40"
Cohesion: 0.15
Nodes (10): _cancel_launched_jobs(), _install_job_cleanup(), _kill_desmond_jobs(), _kill_launched_procs(), OomdGuard, Cancel every still-registered job on the job server. Idempotent (guarded), safe, Force-terminate any Desmond MD process chain still on the GPU for a job THIS run, Terminate every still-registered LOCAL subprocess (SID, MM-GBSA drivers) and its (+2 more)

### Community 41 - "Community 41"
Cohesion: 0.26
Nodes (14): atom(), bond(), _c_band(), _c_clamp(), _c_disc(), _c_dyad(), _c_earrow(), _c_hbond() (+6 more)

### Community 42 - "Community 42"
Cohesion: 0.22
Nodes (13): _draw_rama_background(), Any, _rama_classify(), _rama_palette(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favoured, allowed, and outlier regions., Ramachandran colours from CFG.VIS_RAMA (the single source of truth); the built-i (+5 more)

### Community 43 - "Community 43"
Cohesion: 0.20
Nodes (12): analyse_pi_interactions(), _canonical_resname(), classify_pair(), generate_detailed_interactions(), get_plane_normal(), _ligand_ionisable(), Calculates the optimal best-fit plane normal vector for Pi-stacking analysis usi, Classify a ligand atom as 'anion'- or 'cation'-capable for salt-bridge     detec (+4 more)

### Community 44 - "Community 44"
Cohesion: 0.20
Nodes (12): _block_bootstrap_median_ci(), _cliffs_delta(), _delta_word(), _draw_time_cumulative(), _effective_n(), plot_mmgbsa_combined(), Confidence interval on the median by MOVING-BLOCK bootstrap.      Frames are aut, Cliff's delta: P(a > b) − P(a < b). A non-parametric effect size that needs no (+4 more)

### Community 45 - "Community 45"
Cohesion: 0.18
Nodes (11): graph_map_mmcif_to_rdkit(), kabsch_transform(), map_mmcif_to_rdkit(), Constructs a three-dimensional RDKit molecule from a SMILES string.     This is, Extracts raw XYZ coordination matrices from RDKit molecule data blocks., Determines the optimal rotation/translation matrix required to superimpose spati, Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not b, FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perce (+3 more)

### Community 46 - "Community 46"
Cohesion: 0.22
Nodes (9): atomic_write_csv(), compute_ligand_properties(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, Canonical file logger for all pipeline scripts.      Creates a file-only handler, Write a JSON file so a reader never sees a half-written one.      The write goes, Write a DataFrame to CSV so a reader never sees a half-written file.      Same g, setup_logging() (+1 more)

### Community 47 - "Community 47"
Cohesion: 0.25
Nodes (9): console_info(), console_separator(), console_title(), Logger, Remove all ANSI/VT100 escape sequences from a string., Bold section header - printed and optionally written to log file., Two-space-indented info line - printed and optionally written to log file., Horizontal rule - printed and optionally written to log file.      Parameters (+1 more)

### Community 48 - "Community 48"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 49 - "Community 49"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 50 - "Community 50"
Cohesion: 0.43
Nodes (7): _av(), bonded(), _c_legend_arrow(), _cascade_build(), geom_attack(), _parse(), retained_substituents()

### Community 51 - "Community 51"
Cohesion: 0.29
Nodes (7): _font(), _montage_png(), _pack_rows(), _partition_rows(), Split the ORDERED panels into r contiguous rows so the row aspect-sums are as ba, Lay the panels out as justified rows that FILL the page (like a packed collage),, Justified-row (masonry) raster montage (legacy FIG_FORMAT='png' fallback): read

### Community 52 - "Community 52"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 53 - "Community 53"
Cohesion: 0.33
Nodes (6): analysis_worker_loop(), append_rows_to_csv(), flatten_job_result(), Flattens a process_single_job result dict into a single-level row for CSV writin, Appends a list of result rows to the primary CSV, deduplicates, sorts by job_nam, Background worker that analyses freshly GPU-predicted jobs and appends rows to t

### Community 54 - "Community 54"
Cohesion: 0.33
Nodes (6): _esp_run(), generate_esp_charges(), plot_esp_alpha_carbon(), Run a snippet under $SCHRODINGER/run - its interpreter, not ours., QM (Jaguar ESP) partial charges for every prepared MD ligand. Returns the summar, The charge on the carbon the nucleophile actually attacks - the number OPLS4 can

### Community 55 - "Community 55"
Cohesion: 0.50
Nodes (4): calculate_angle(), calculate_burgi_dunitz(), Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for, Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit

### Community 56 - "Community 56"
Cohesion: 0.50
Nodes (4): analyse_candidate_structure(), load_structure_safe(), Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilit, Superimposes the candidate structure onto the DeHa4 control model.     Computes

### Community 57 - "Community 57"
Cohesion: 0.50
Nodes (4): analyse_model_task(), Global worker routine functionally decoupled for robust serialisation capabiliti, Selects the representative pose across the Boltz diffusion samples tier-first: t, select_best_degrader_model()

### Community 58 - "Community 58"
Cohesion: 0.50
Nodes (4): compute_cross_interface_pae(), load_extra_boltz_metrics(), Calculates the Predicted Aligned Error (PAE), focusing specifically on the prote, Parses background statistical validation tensors directly from Boltz NPZ file du

### Community 59 - "Community 59"
Cohesion: 0.50
Nodes (4): map_active_site_residues(), Append a newly computed alignment to the on-disk cache file, under an exclusive, Performs a global sequence alignment between the canonical reference sequence, update_alignment_stats()

### Community 60 - "Community 60"
Cohesion: 0.50
Nodes (4): _build_hd1_line(), enforce_catalytic_protonation(), Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to, Give each catalytic residue the protonation its ROLE requires (CFG §15).      Ev

### Community 61 - "Community 61"
Cohesion: 0.50
Nodes (4): _add_label(), build_pptx(), _emu(), Assemble the editable PowerPoint deck from the Figure_*.svg already written next

## Knowledge Gaps
- **4 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `Community 6` to `Community 0`, `Community 1`, `Community 2`, `Community 3`, `Community 4`, `Community 5`, `Community 7`, `Community 8`, `Community 9`, `Community 10`, `Community 11`, `Community 13`, `Community 14`, `Community 15`, `Community 16`, `Community 17`, `Community 20`, `Community 22`, `Community 23`, `Community 26`, `Community 30`, `Community 37`, `Community 38`, `Community 39`, `Community 44`, `Community 54`, `Community 59`?**
  _High betweenness centrality (0.543) - this node is a cross-community bridge._
- **Why does `main()` connect `Community 8` to `Community 6`, `Community 53`, `Community 56`, `Community 27`, `Community 28`, `Community 29`?**
  _High betweenness centrality (0.097) - this node is a cross-community bridge._
- **Why does `_generate_comprehensive_figures_impl()` connect `Community 1` to `Community 2`, `Community 34`, `Community 4`, `Community 6`, `Community 7`, `Community 13`, `Community 15`?**
  _High betweenness centrality (0.074) - this node is a cross-community bridge._
- **Are the 75 inferred relationships involving `CFG` (e.g. with `generate_plots()` and `calculate_sn2_metrics()`) actually correct?**
  _`CFG` has 75 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `main()` (e.g. with `print_elapsed()` and `print_script_banner()`) actually correct?**
  _`main()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _522 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.046464646464646465 - nodes in this community are weakly interconnected._