# Graph Report - FAcDs pipeline  (2026-07-10)

## Corpus Check
- 13 files · ~156,974 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 605 nodes · 1222 edges · 38 communities (36 shown, 2 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 68 edges (avg confidence: 0.74)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Logging & Report Manager|Logging & Report Manager]]
- [[_COMMUNITY_02 Production Catalytic Geometry & SN2|02 Production: Catalytic Geometry & SN2]]
- [[_COMMUNITY_02 Production Role & Metric Derivation|02 Production: Role & Metric Derivation]]
- [[_COMMUNITY_02 Production Analysis Workers|02 Production: Analysis Workers]]
- [[_COMMUNITY_02 Production Model Analysis & PAE|02 Production: Model Analysis & PAE]]
- [[_COMMUNITY_05 Top-N CIF→PDB Preparation|05 Top-N: CIF→PDB Preparation]]
- [[_COMMUNITY_05 Top-N Interaction Maps|05 Top-N: Interaction Maps]]
- [[_COMMUNITY_04 Dendrogram & Phylogeny|04 Dendrogram & Phylogeny]]
- [[_COMMUNITY_05 Top-N PLIP Interaction Diagrams|05 Top-N: PLIP Interaction Diagrams]]
- [[_COMMUNITY_06 Prime MM-GBSA|06 Prime MM-GBSA]]
- [[_COMMUNITY_06 SID Job Processing|06 SID Job Processing]]
- [[_COMMUNITY_Project Config (CFG)|Project Config (CFG)]]
- [[_COMMUNITY_01 Merge & QC|01 Merge & QC]]
- [[_COMMUNITY_Utils Geometry & Console|Utils: Geometry & Console]]
- [[_COMMUNITY_07 Defluorination Figures|07 Defluorination Figures]]
- [[_COMMUNITY_05 Top-N Reference & Map Loading|05 Top-N: Reference & Map Loading]]
- [[_COMMUNITY_07 NAC & QSite Core|07 NAC & QSite Core]]
- [[_COMMUNITY_07 QSite Input Prep|07 QSite Input Prep]]
- [[_COMMUNITY_02 Production Resume & MSA|02 Production: Resume & MSA]]
- [[_COMMUNITY_07 Per-Job Analysis & Dashboards|07 Per-Job Analysis & Dashboards]]
- [[_COMMUNITY_Cluster 20|Cluster 20]]
- [[_COMMUNITY_Cluster 21|Cluster 21]]
- [[_COMMUNITY_Cluster 22|Cluster 22]]
- [[_COMMUNITY_Cluster 23|Cluster 23]]
- [[_COMMUNITY_Cluster 24|Cluster 24]]
- [[_COMMUNITY_Cluster 25|Cluster 25]]
- [[_COMMUNITY_Cluster 26|Cluster 26]]
- [[_COMMUNITY_Cluster 27|Cluster 27]]
- [[_COMMUNITY_Cluster 28|Cluster 28]]
- [[_COMMUNITY_Cluster 29|Cluster 29]]
- [[_COMMUNITY_Cluster 30|Cluster 30]]
- [[_COMMUNITY_Cluster 31|Cluster 31]]
- [[_COMMUNITY_Cluster 32|Cluster 32]]
- [[_COMMUNITY_Cluster 33|Cluster 33]]
- [[_COMMUNITY_Cluster 34|Cluster 34]]
- [[_COMMUNITY_Cluster 35|Cluster 35]]
- [[_COMMUNITY_Cluster 36|Cluster 36]]
- [[_COMMUNITY_Cluster 37|Cluster 37]]

## God Nodes (most connected - your core abstractions)
1. `CFG` - 44 edges
2. `main()` - 37 edges
3. `process_single_job()` - 23 edges
4. `process_single_job()` - 21 edges
5. `_generate_comprehensive_figures_impl()` - 21 edges
6. `prep_and_convert_phase()` - 18 edges
7. `check_catalytic_geometry()` - 17 edges
8. `ReportManager` - 16 edges
9. `run_figure_generation()` - 15 edges
10. `_echo()` - 15 edges

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

## Communities (38 total, 2 thin omitted)

### Community 0 - "Logging & Report Manager"
Cohesion: 0.07
Nodes (64): Simultaneous console + file logger shared by the pipeline steps.      The log pa, ReportManager, analyse_conflicts(), _auto_label_colour(), _aux_dir(), calculate_pareto_fronts(), console_info(), _fig23_multitarget() (+56 more)

### Community 1 - "02 Production: Catalytic Geometry & SN2"
Cohesion: 0.09
Nodes (34): analyse_candidate_structure(), analyse_pi_interactions(), calculate_sn2_metrics(), _canonical_resname(), check_catalytic_geometry(), classify_pair(), compute_pocket_fit(), generate_detailed_interactions() (+26 more)

### Community 2 - "02 Production: Role & Metric Derivation"
Cohesion: 0.07
Nodes (26): _derive_burgi_dunitz(), _derive_flippin_lodge(), format_control_mappings(), format_full_role_map(), generate_rich_justification(), get_cached_alignment_for_protein(), get_plane_normal(), init_gpu_reservations() (+18 more)

### Community 3 - "02 Production: Analysis Workers"
Cohesion: 0.11
Nodes (27): analysis_worker_loop(), console_info(), console_separator(), cpu_usage_summary(), extract_short_fasta_id(), flatten_job_result(), generate_scientific_ranking_csv(), load_cached_alignments() (+19 more)

### Community 4 - "02 Production: Model Analysis & PAE"
Cohesion: 0.08
Nodes (26): analyse_model_task(), check_job_status(), colabfold_a3m_path(), colabfold_meta_path(), compute_cross_interface_pae(), extract_3r3u_sequence(), heal_smiles_in_files(), load_structure_safe() (+18 more)

### Community 5 - "05 Top-N: CIF→PDB Preparation"
Cohesion: 0.10
Nodes (26): check_prep_needed(), _check_residue_identity_guard(), cif_to_pdb_gemmi(), collect_best_cifs(), generate_raw_step(), get_source_tag(), index_existing_files(), load_md_selected_jobs() (+18 more)

### Community 6 - "05 Top-N: Interaction Maps"
Cohesion: 0.12
Nodes (23): _append_auxiliary_log(), _im_don_acc(), _im_parse_pdb(), _im_project(), _im_render_diagram(), _im_ring(), _im_separate_atoms(), _parse_plip_xml() (+15 more)

### Community 7 - "04 Dendrogram & Phylogeny"
Cohesion: 0.17
Nodes (18): clean_id(), console_info(), console_separator(), generate_phylogenies(), generate_upgma_newick(), get_kmer_counts(), _load_module(), main() (+10 more)

### Community 8 - "05 Top-N: PLIP Interaction Diagrams"
Cohesion: 0.16
Nodes (13): _draw_interaction_diagram(), _fig_constants(), _fig_log(), Render a publication-quality 2D protein–ligand interaction map (no PyMOL)., Handles verification and automated installation of visual tools., Check if PLIP is installed as a Python module., Return (res_name, res_num, has_f) for the ligand in chain L.          Skips AA/U, Phase 2: Run PyMOL and PLIP rendering on all PDB files under ext_dir.     Called (+5 more)

### Community 9 - "06 Prime MM-GBSA"
Cohesion: 0.16
Nodes (19): _avg_dg(), _boltzmann_mean_dg(), _lookup_tiers(), main(), _mmgbsa_dg_series(), _natural_rank(), plot_mmgbsa_combined(), plot_mmgbsa_individual() (+11 more)

### Community 10 - "06 SID Job Processing"
Cohesion: 0.15
Nodes (20): is_eaf_complete(), _load_module(), out_eaf_frames(), _proc_alive(), process_jobs(), Path, _rank_of(), Count frames written into a SID-out.eaf (the analysed result vector).      Parse (+12 more)

### Community 11 - "Project Config (CFG)"
Cohesion: 0.11
Nodes (7): CFG, ===============================================================================, Central Configuration Repository — FAcDs Pipeline.      Attribute prefix → secti, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant, Single source of truth for the holistic mechanistic score (0–1) — the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the

### Community 12 - "01 Merge & QC"
Cohesion: 0.18
Nodes (18): apply_clean_spines(), clean_header(), clean_sequence_str(), generate_plots(), is_valid_protein(), _load_module(), main(), process_and_write() (+10 more)

### Community 13 - "Utils: Geometry & Console"
Cohesion: 0.11
Nodes (17): calculate_flippin_lodge(), ConsoleColours, find_nucleophile_od_fallback(), get_alignment_grade(), print_elapsed(), print_script_banner(), ===============================================================================, Print a coloured startup banner identifying which script is running.      Parame (+9 more)

### Community 14 - "07 Defluorination Figures"
Cohesion: 0.17
Nodes (17): clean_spines(), Academic-style axes: remove top/right spines, thin the remaining borders.      C, console_separator(), console_title(), format_job_label(), generate_defluorination_landscape(), generate_global_comparative_dashboard(), generate_viability_bar_chart() (+9 more)

### Community 15 - "05 Top-N: Reference & Map Loading"
Cohesion: 0.15
Nodes (17): console_info(), console_separator(), extract_chain_l_mol(), load_catalytic_anchor_map(), load_rank_map(), load_reference_data(), main(), prep_and_convert_phase() (+9 more)

### Community 16 - "07 NAC & QSite Core"
Cohesion: 0.14
Nodes (15): calculate_min_distance(), compute_wm_csv_stats(), console_qmm_ready(), extract_hybrid_smart_system(), _nac_dwell_stats(), parse_mapping(), _print_labeled(), Continuous-residence statistics for a per-frame strict-NAC boolean series. (+7 more)

### Community 17 - "07 QSite Input Prep"
Cohesion: 0.15
Nodes (13): find_eaf_file(), generate_qsite_inputs(), _load_module(), load_triad_mapping(), Path, Write an uncompressed QSite .mae trimmed to a solvation droplet: the full     pr, Write a valid QSite/Jaguar QM/MM relaxed-scan .in for the SN2     dehalogenation, Load a Python file as a module regardless of its filename. (+5 more)

### Community 18 - "02 Production: Resume & MSA"
Cohesion: 0.17
Nodes (12): deep_rename_job_folder(), fetch_msa_direct(), job_key_from_yaml(), Generates a stable hash for FASTA sequences to verify consistency during resumed, Extracts the stable resume key from an existing YAML definition file., Renames existing YAMLs and run folders to match the current input ordering durin, Recursively updates filenames within a job folder to reflect changes in identifi, Validates the A3M MSA file to detect corruption. Returns True if the file struct (+4 more)

### Community 19 - "07 Per-Job Analysis & Dashboards"
Cohesion: 0.21
Nodes (12): _blockade_vec(), _eaf_at(), generate_individual_dashboard(), load_eaf_scalar_series(), process_single_job(), DataFrame, ndarray, Vectorized water blockade — replaces check_water_blockade when positions     are (+4 more)

### Community 20 - "Cluster 20"
Cohesion: 0.24
Nodes (11): _draw_rama_background(), Any, _rama_classify(), _rama_stats(), Classify a phi/psi pair as Favored, Allowed, or Outlier., Calculate percentages of residues in favored, allowed, and outlier regions., Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes., Save a side-by-side comparison Ramachandran PNG. (+3 more)

### Community 21 - "Cluster 21"
Cohesion: 0.18
Nodes (11): Remove all ANSI/VT100 escape sequences from a string., _strip_ansi(), console_info(), extract_8residue_indices(), load_watermap_csv(), load_watermap_sites(), Launch the QSite/Jaguar executable on a generated .in, writing all output     in, Load WaterMap thermodynamic sites from a Maestro 'Analyse WaterMap'     CSV expo (+3 more)

### Community 22 - "Cluster 22"
Cohesion: 0.20
Nodes (10): calculate_angle(), calculate_burgi_dunitz(), calculate_improper_dihedral(), distance(), get_mic_vector(), MIC displacement vector from pos2 → pos1, supporting orthorhombic and     tricli, PBC-corrected or Euclidean distance between two Cartesian coordinates (Å)., Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.      Returns 0.0 if eit (+2 more)

### Community 23 - "Cluster 23"
Cohesion: 0.20
Nodes (10): append_rows_to_csv(), atomic_to_csv(), Scan every completed job's Best_Complex/ folder and copy each CIF into     the f, Appends a list of result rows to the master CSV, deduplicates, sorts by job_name, Single parallel scan over all *_summary.json files. Returns n_rows_written., Crash-safe CSV write: serialise to a sibling .tmp then os.replace() onto the, Carriage-return progress write. Always writes raw so ANSI colours pass     throu, rebuild_best_complexes_mirror() (+2 more)

### Community 24 - "Cluster 24"
Cohesion: 0.24
Nodes (4): Heartbeat, Periodically report a long Schrödinger step's progress by tailing its log., Largest last-match across all configured patterns (phase-tolerant)., Best-effort (done, total) for the Prime phase from a ``*-prime*.log`` in

### Community 25 - "Cluster 25"
Cohesion: 0.50
Nodes (7): _print_timing_table(), _run_all_steps(), run_step(), 00_00_run_pipeline_FAcDs.sh script, _sync_staging_log(), _tee(), WARN_EXIT_CODE

### Community 26 - "Cluster 26"
Cohesion: 0.31
Nodes (9): _box_diag(), calculate_min_distance(), _ensure_box_3x3(), mic_dists_2d(), ndarray, Minimum PBC-corrected distance between two sets of Cartesian positions., Coerce a periodic box to a (3,3) float64 matrix.     Handles: (3,3) arrays, flat, Orthorhombic fallback box lengths as a safe (3,) diagonal.      Used when the fu (+1 more)

### Community 27 - "Cluster 27"
Cohesion: 0.25
Nodes (9): compute_ligand_properties(), ensure_jobserver_on_working_disk(), _fs_device(), _jobserver_addresses(), Path, Compute per-ligand carbon count (nC), fluorine count (nF) and molecular     weig, st_dev of the nearest existing ancestor of `path` (identifies its mounted     fi, Registered local job-server addresses (e.g. ['localhost:40931']) parsed from (+1 more)

### Community 28 - "Cluster 28"
Cohesion: 0.22
Nodes (9): console_info(), console_separator(), console_title(), Logger, Bold section header — printed and optionally written to log file., Two-space-indented info line — printed and optionally written to log file., Horizontal rule — printed and optionally written to log file.      Parameters, Canonical file logger for all pipeline scripts.      Creates a file-only handler (+1 more)

### Community 29 - "Cluster 29"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), main(), Creates a fresh Conda environment from the exported PFAS.yml., Verifies installed packages in the current environment and prints a checklist., Exports the current active Conda environment to PFAS.yml and requirements.txt., verify_environment()

### Community 30 - "Cluster 30"
Cohesion: 0.29
Nodes (8): _im_angle(), _im_arom_hbond_check(), _im_find_contacts(), _im_hbond_check(), Angle (degrees) at vertex b for points a-b-c., True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.      Ret, Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).      Dista, Detect interactions using the Maestro criteria defined in CFG §3.      Priority

### Community 31 - "Cluster 31"
Cohesion: 0.32
Nodes (4): _echo(), OomdGuard, Print to terminal immediately (flush) — keeps live progress visible., Optionally mask systemd-oomd to prevent Out-Of-Memory kills during long SID runs

### Community 32 - "Cluster 32"
Cohesion: 0.25
Nodes (8): _extract_fluoride_charge_series(), _extract_scan_energies(), parse_qsite_barrier(), parse_qsite_profile(), Best-effort per-scan-point energy series (Hartree) from a Jaguar/QSite     relax, Best-effort ordered list of the most-negative Mulliken fluorine charge per     p, Parse the QM/MM SN2 relaxed scan into a full reaction profile: the barrier     a, Scalar barrier/reaction-energy/fluoride subset of parse_qsite_profile, for     t

### Community 33 - "Cluster 33"
Cohesion: 0.29
Nodes (6): _diagnose_mmgbsa_failure(), _mmgbsa_csv(), Locate a thermal_mmgbsa results CSV in `job_dir` (name varies by version).     R, Best-effort human-readable cause when MM-GBSA exits non-zero, read from the, Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:     ret, run_mmgbsa()

### Community 34 - "Cluster 34"
Cohesion: 0.33
Nodes (6): calculate_dihedral(), compute_ramachandran_angles(), _rama_get_atom_pos(), Find atom coordinates in a Gemmi residue structure., Extract (resname, resnum, phi, psi) for every residue that has both angles., Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.     Uses the Gr

### Community 36 - "Cluster 36"
Cohesion: 0.50
Nodes (4): Initialise the Top-N extraction log (separate handler from the prep log)., Initialises the preparation log via the shared utility., setup_logging(), setup_logging_extraction()

## Knowledge Gaps
- **2 isolated node(s):** `WARN_EXIT_CODE`, `git_push_FAcDs.sh script`
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `Project Config (CFG)` to `Logging & Report Manager`, `Cluster 33`, `02 Production: Catalytic Geometry & SN2`, `02 Production: Analysis Workers`, `Cluster 32`, `05 Top-N: CIF→PDB Preparation`, `04 Dendrogram & Phylogeny`, `06 Prime MM-GBSA`, `07 Defluorination Figures`, `05 Top-N: Reference & Map Loading`, `07 Per-Job Analysis & Dashboards`, `Cluster 21`?**
  _High betweenness centrality (0.574) - this node is a cross-community bridge._
- **Why does `main()` connect `02 Production: Analysis Workers` to `02 Production: Catalytic Geometry & SN2`, `02 Production: Role & Metric Derivation`, `02 Production: Model Analysis & PAE`, `Project Config (CFG)`, `Utils: Geometry & Console`, `02 Production: Resume & MSA`, `Cluster 23`?**
  _High betweenness centrality (0.165) - this node is a cross-community bridge._
- **Why does `topn_extraction_phase()` connect `05 Top-N: Reference & Map Loading` to `Cluster 34`, `Cluster 36`, `05 Top-N: Interaction Maps`, `05 Top-N: PLIP Interaction Diagrams`, `Project Config (CFG)`, `Cluster 20`?**
  _High betweenness centrality (0.098) - this node is a cross-community bridge._
- **Are the 30 inferred relationships involving `CFG` (e.g. with `calculate_sn2_metrics()` and `check_catalytic_geometry()`) actually correct?**
  _`CFG` has 30 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `main()` (e.g. with `CFG` and `safe_name()`) actually correct?**
  _`main()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `process_single_job()` (e.g. with `CFG` and `get_mic_vector()`) actually correct?**
  _`process_single_job()` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `===============================================================================`, `Central Configuration Repository — FAcDs Pipeline.      Attribute prefix → secti` to the rest of the system?**
  _248 weakly-connected nodes found - possible documentation gaps or missing edges._