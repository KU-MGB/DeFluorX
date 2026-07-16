#!/usr/bin/env python3
"""Build the pushable code-knowledge-graph docs from the (local, un-pushed) graphify-out/.

Run this after every `/graphify` update so docs/ always holds the current graph:
    python docs/build_code_graph.py

Outputs (all under docs/, all committed):
  • code_graph.html         — self-contained interactive viewer (data embedded, no CDN, works offline)
  • code_graph_report.md    — copy of graphify's GRAPH_REPORT.md
  • code_graph.png          — a static snapshot of the viewer (rendered separately; see note below)

The PNG needs a browser render. If a headless Chrome/Chromium is on PATH this script makes it;
otherwise it prints how to snapshot code_graph.html manually (the .html and .md are always refreshed).
"""
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GOUT = ROOT / "graphify-out"
DOCS = ROOT / "docs"
TEMPLATE = DOCS / "code_graph_template.html"


def compact(graph: dict) -> dict:
    """graphify graph.json → the minimal {nodes, links} the viewer embeds."""
    nodes = graph["nodes"]
    links = graph.get("links", graph.get("edges", []))
    idx = {n["id"]: i for i, n in enumerate(nodes)}
    deg = Counter()
    for e in links:
        deg[e["source"]] += 1
        deg[e["target"]] += 1
    sf = lambda n: str(n.get("source_file", "")).split("/")[-1]
    out = {
        "nodes": [{"id": n["id"], "label": (n.get("label") or n["id"])[:48], "f": sf(n),
                   "t": n.get("file_type", n.get("node_type", "")), "d": int(deg[n["id"]])}
                  for n in nodes],
        "links": [{"s": idx[e["source"]], "t": idx[e["target"]], "r": e.get("relation", "")}
                  for e in links if e["source"] in idx and e["target"] in idx],
    }
    # Synthetic orchestration edges: 00_00 (bash) runs each pipeline step. The AST is python-import
    # based, so it never sees the shell's `python NN_*.py` calls and the orchestrator would float with
    # zero edges. Add them so the main orchestrator visibly connects to what it drives.
    mod = {n["f"]: i for i, n in enumerate(out["nodes"]) if n["label"] == n["f"]}
    orch = mod.get("00_00_run_pipeline_FAcDs.sh")
    steps = ["00_03_Environment_FAcDs.py", "01_Merge_FAcDs.py", "02_Production_FAcDs.py",
             "03_Validation_Figures_FAcDs.py", "04_Dendrogram_FAcDs.py",
             "05_TopN_and_PDB_Preparation_FAcDs.py", "06_Physics_Validation_FAcDs.py",
             "07_MD_QMMM_Defluorination_FAcDs.py"]
    def link(a, b, r):
        out["links"].append({"s": a, "t": b, "r": r})
        out["nodes"][a]["d"] += 1
        out["nodes"][b]["d"] += 1

    if orch is not None:
        for s in steps:
            if s in mod:
                link(orch, mod[s], "runs")
        if "git_push_FAcDs.sh" in mod:
            link(orch, mod["git_push_FAcDs.sh"], "runs")

    # Bridge the README doc-cluster into the code graph so it does not float as a separate component:
    # its semantic concepts have no edges to code, so without this the whole 100+ node README cluster
    # drifts off on its own. Anchor on the busiest README concept and link it to what it documents.
    readme = [i for i, n in enumerate(out["nodes"]) if n["f"] == "README.md"]
    if readme:
        anchor = max(readme, key=lambda i: out["nodes"][i]["d"])
        for s in steps + ["00_00_run_pipeline_FAcDs.sh", "00_01_Project_Config_FAcDs.py",
                          "00_02_Project_Utils_FAcDs.py"]:
            if s in mod:
                link(anchor, mod[s], "documents")
    return out


def build_html() -> Path:
    graph = json.loads((GOUT / "graph.json").read_text())
    data = compact(graph)
    html = TEMPLATE.read_text().replace("/*__DATA__*/", json.dumps(data))
    out = DOCS / "code_graph.html"
    out.write_text(html)
    print(f"  code_graph.html      — {len(data['nodes'])} nodes, {len(data['links'])} links, {round(len(html)/1024)} KB")
    return out


def copy_report() -> None:
    src = GOUT / "GRAPH_REPORT.md"
    if src.exists():
        shutil.copy2(src, DOCS / "code_graph_report.md")
        print("  code_graph_report.md — copied from graphify-out/GRAPH_REPORT.md")


def render_png(html: Path) -> None:
    png = DOCS / "code_graph.png"
    # 1) playwright headless chromium (preferred — deterministic, waits for the layout to settle)
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            br = p.chromium.launch()
            pg = br.new_page(viewport={"width": 1920, "height": 1200}, device_scale_factor=2)
            pg.goto(html.as_uri())
            pg.wait_for_timeout(6000)   # let the force layout settle + auto-fit fire
            pg.screenshot(path=str(png))
            br.close()
        print("  code_graph.png       — rendered with playwright chromium")
        return
    except Exception as e:
        print(f"  code_graph.png       — playwright unavailable ({str(e).splitlines()[0][:60]}); trying system chrome…")
    # 2) system chrome/chromium CLI fallback
    for b in ("chromium", "chromium-browser", "google-chrome", "google-chrome-stable", "chrome"):
        if shutil.which(b):
            png = DOCS / "code_graph.png"
            try:
                subprocess.run([b, "--headless=new", "--hide-scrollbars", "--force-device-scale-factor=1.5",
                                "--window-size=1920,1200", f"--screenshot={png}", "--virtual-time-budget=4000",
                                html.as_uri()], check=True, timeout=90, capture_output=True)
                print(f"  code_graph.png       — rendered with {b}")
                return
            except Exception as e:
                print(f"  code_graph.png       — {b} render failed: {e}")
    print("  code_graph.png       — SKIPPED (no headless Chrome/Chromium on PATH). "
          "Open docs/code_graph.html in a browser, 'Fit view', and save a screenshot as docs/code_graph.png.")


def main() -> None:
    if not (GOUT / "graph.json").exists():
        raise SystemExit("graphify-out/graph.json not found — run /graphify first.")
    print("Building code-graph docs from graphify-out/ …")
    html = build_html()
    copy_report()
    render_png(html)
    print("Done. docs/ now holds the current code graph (html + md + png).")


if __name__ == "__main__":
    main()
