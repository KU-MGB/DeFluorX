# Graph Report - .  (2026-08-05)

## Corpus Check
- 12 files · ~6,237 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1140 nodes · 2512 edges · 12 communities (11 shown, 1 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 174 edges (avg confidence: 0.72)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_00_00_run_pipeline|00_00_run_pipeline]]
- [[_COMMUNITY_00_01_Project_Config|00_01_Project_Config]]
- [[_COMMUNITY_00_02_Project_Utils|00_02_Project_Utils]]
- [[_COMMUNITY_00_03_Environment|00_03_Environment]]
- [[_COMMUNITY_01_Merge|01_Merge]]
- [[_COMMUNITY_02_Production|02_Production]]
- [[_COMMUNITY_03_Validation_Figures|03_Validation_Figures]]
- [[_COMMUNITY_04_Dendrogram|04_Dendrogram]]
- [[_COMMUNITY_05_TopN_and_PDB_Preparation|05_TopN_and_PDB_Preparation]]
- [[_COMMUNITY_06_Physics_Validation|06_Physics_Validation]]
- [[_COMMUNITY_07_MD_QMMM_Defluorination|07_MD_QMMM_Defluorination]]
- [[_COMMUNITY_git_push|git_push]]

## God Nodes (most connected - your core abstractions)
1. `03_Validation_Figures_DeFluorX.py` - 135 edges
2. `06_Physics_Validation_DeFluorX.py` - 114 edges
3. `07_MD_QMMM_Defluorination_DeFluorX.py` - 104 edges
4. `CFG` - 84 edges
5. `02_Production_DeFluorX.py` - 66 edges
6. `05_TopN_and_PDB_Preparation_DeFluorX.py` - 64 edges
7. `Path` - 59 edges
8. `Path` - 46 edges
9. `00_02_Project_Utils_DeFluorX.py` - 43 edges
10. `DataFrame` - 41 edges
11. `main()` - 36 edges
12. `main()` - 34 edges
13. `Path` - 34 edges
14. `process_single_job()` - 34 edges
15. `Path` - 32 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Import Cycles
- None detected.

## Communities (12 total, 1 thin omitted)

### Community 0 - "00_00_run_pipeline"
Cohesion: 0.22
Nodes (13): 00_00_run_pipeline_DeFluorX.sh script, WARN_EXIT_CODE, _print_step_menu(), _sync_staging_log(), _tee(), run_step(), _TT_W_NAME, _TT_W_TIME (+5 more)

### Community 1 - "00_01_Project_Config"
Cohesion: 0.10
Nodes (7): CFG, ===============================================================================, Central Configuration Repository - DeFluorX Pipeline.      Attribute prefix → se, Single source of truth for the holistic mechanistic score (0–1) - the raw geomet, Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic, Single source of truth for the gated continuous competence score (0–1), the, Internal-consistency guards (read-only; frozen-dataclass safe). Several constant

### Community 2 - "00_02_Project_Utils"
Cohesion: 0.03
Nodes (90): ConsoleColours, _ConsoleRuleFilter, install_console_rule_filter(), safe_name(), setup_logging(), _strip_ansi(), console_title(), console_info() (+82 more)

### Community 3 - "00_03_Environment"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), verify_environment(), main(), Exports the current active Conda environment to PFAS.yml and requirements.txt., Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing.

### Community 4 - "01_Merge"
Cohesion: 0.18
Nodes (18): _load_module(), setup_logger(), clean_sequence_str(), is_valid_protein(), clean_header(), process_and_write(), apply_clean_spines(), generate_plots() (+10 more)

### Community 5 - "02_Production"
Cohesion: 0.03
Nodes (129): _load_module(), _canonical_resname(), setup_logging(), console_info(), console_separator(), atomic_to_csv(), _tty_write(), extract_short_fasta_id() (+121 more)

### Community 6 - "03_Validation_Figures"
Cohesion: 0.02
Nodes (235): _load_module(), _tier_seps(), _stat_box(), _register_p(), _statistical_battery(), _write_statistical_tests(), _kruskal_by_tier(), _wilcoxon_ptm_iptm() (+227 more)

### Community 7 - "04_Dendrogram"
Cohesion: 0.18
Nodes (18): _load_module(), console_info(), console_separator(), _make_reporter(), clean_id(), get_kmer_counts(), generate_upgma_newick(), package_deployment() (+10 more)

### Community 8 - "05_TopN_and_PDB_Preparation"
Cohesion: 0.03
Nodes (121): _load_module(), setup_logging(), console_info(), console_separator(), index_existing_files(), collect_best_cifs(), load_rank_map(), load_md_selected_jobs() (+113 more)

### Community 9 - "06_Physics_Validation"
Cohesion: 0.02
Nodes (225): _load_module(), _open_step_log(), _progress_line(), _close_bar(), _worker_progress(), _echo(), _log(), _ok() (+217 more)

### Community 10 - "07_MD_QMMM_Defluorination"
Cohesion: 0.02
Nodes (198): LazyTrajectory, _load_module(), console_title(), console_info(), _mask_oomd_at_start(), console_separator(), _print_labeled(), console_qmm_ready() (+190 more)

## Knowledge Gaps
- **5 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `git_push_DeFluorX.sh script`
  These have ≤1 connection - possible missing edges or undocumented components.
- **1 thin communities (<3 nodes) omitted from report** - run `graphify query` to explore isolated nodes.