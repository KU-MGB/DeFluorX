# Graph Report - .  (2026-07-29)

## Corpus Check
- 12 files · ~0 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 1075 nodes · 2099 edges · 12 communities (11 shown, 1 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 150 edges (avg confidence: 0.71)
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
1. `CFG` - 91 edges
2. `main()` - 33 edges
3. `main()` - 33 edges
4. `_echo()` - 31 edges
5. `_generate_comprehensive_figures_impl()` - 30 edges
6. `process_single_job()` - 30 edges
7. `console_info()` - 27 edges
8. `prep_and_convert_phase()` - 25 edges
9. `run_mmgbsa_phase()` - 21 edges
10. `process_single_job()` - 19 edges

## Surprising Connections (you probably didn't know these)
- `calculate_sn2_metrics()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `check_catalytic_geometry()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `generate_scientific_ranking_csv()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `main()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py
- `map_active_site_residues()` --indirect_call--> `CFG`  [INFERRED]
  02_Production_DeFluorX.py → 00_01_Project_Config_DeFluorX.py

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
Nodes (84): ConsoleColours, _ConsoleRuleFilter, install_console_rule_filter(), safe_name(), setup_logging(), _strip_ansi(), console_title(), console_info() (+76 more)

### Community 3 - "00_03_Environment"
Cohesion: 0.33
Nodes (8): ConsoleColours, export_environment(), install_environment(), verify_environment(), main(), Exports the current active Conda environment to PFAS.yml and requirements.txt., Creates a fresh Conda environment from the exported PFAS.yml., Verify the installed pipeline packages and HALT if a mandatory one is missing.

### Community 4 - "01_Merge"
Cohesion: 0.17
Nodes (15): setup_logger(), clean_sequence_str(), is_valid_protein(), clean_header(), process_and_write(), apply_clean_spines(), generate_plots(), main() (+7 more)

### Community 5 - "02_Production"
Cohesion: 0.03
Nodes (126): _load_module(), _canonical_resname(), setup_logging(), console_info(), console_separator(), atomic_to_csv(), _tty_write(), extract_short_fasta_id() (+118 more)

### Community 6 - "03_Validation_Figures"
Cohesion: 0.02
Nodes (233): _tier_seps(), _stat_box(), _register_p(), _statistical_battery(), _write_statistical_tests(), _kruskal_by_tier(), _wilcoxon_ptm_iptm(), _umap_tier_separation() (+225 more)

### Community 7 - "04_Dendrogram"
Cohesion: 0.18
Nodes (15): console_info(), console_separator(), _make_reporter(), clean_id(), get_kmer_counts(), generate_upgma_newick(), package_deployment(), generate_phylogenies() (+7 more)

### Community 8 - "05_TopN_and_PDB_Preparation"
Cohesion: 0.03
Nodes (114): setup_logging(), console_info(), console_separator(), index_existing_files(), collect_best_cifs(), load_rank_map(), load_md_selected_jobs(), get_source_tag() (+106 more)

### Community 9 - "06_Physics_Validation"
Cohesion: 0.02
Nodes (221): _open_step_log(), _progress_line(), _close_bar(), _worker_progress(), _echo(), _log(), _ok(), _fail() (+213 more)

### Community 10 - "07_MD_QMMM_Defluorination"
Cohesion: 0.02
Nodes (156): LazyTrajectory, _load_module(), console_title(), console_info(), _mask_oomd_at_start(), console_separator(), _print_labeled(), console_qmm_ready() (+148 more)

## Knowledge Gaps
- **5 isolated node(s):** `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME`, `_TT_W_STAT`, `git_push_DeFluorX.sh script`
  These have ≤1 connection - possible missing edges or undocumented components.
- **1 thin communities (<3 nodes) omitted from report** - run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CFG` connect `00_01_Project_Config` to `00_00_run_pipeline`, `00_02_Project_Utils`, `00_03_Environment`, `01_Merge`, `02_Production`, `03_Validation_Figures`, `04_Dendrogram`, `05_TopN_and_PDB_Preparation`, `06_Physics_Validation`, `07_MD_QMMM_Defluorination`?**
  _High betweenness centrality (0.343) - this node is a cross-community bridge._
- **Why does `main()` connect `02_Production` to `00_01_Project_Config`, `00_02_Project_Utils`?**
  _High betweenness centrality (0.017) - this node is a cross-community bridge._
- **Are the 66 inferred relationships involving `CFG` (e.g. with `calculate_sn2_metrics()` and `check_catalytic_geometry()`) actually correct?**
  _`CFG` has 66 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `main()` (e.g. with `CFG` and `safe_name()`) actually correct?**
  _`main()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `main()` (e.g. with `print_elapsed()` and `print_script_banner()`) actually correct?**
  _`main()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 19 inferred relationships involving `_generate_comprehensive_figures_impl()` (e.g. with `CFG` and `_jf_circos()`) actually correct?**
  _`_generate_comprehensive_figures_impl()` has 19 INFERRED edges - model-reasoned connections that need verification._
- **What connects `WARN_EXIT_CODE`, `_TT_W_NAME`, `_TT_W_TIME` to the rest of the system?**
  _469 weakly-connected nodes found - possible documentation gaps or missing edges._