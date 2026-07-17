# Pipeline architecture diagram

`pipeline_architecture.png` (embedded in the top-level README) is a post-processed render,
not a live Mermaid block — Mermaid cannot rotate subgraph labels or draw a numbered red
circle, and GitHub will not render `<foreignObject>` labels from an SVG loaded via `<img>`,
so the diagram is shipped as a PNG.

## Regenerate (after editing the diagram or updating per-script line counts)

1. Edit the source: `pipeline_architecture.mmd`.
2. Render + post-process + rasterise:

```bash
# needs: npx @mermaid-js/mermaid-cli, a headless Chrome, rsvg (optional)
mmdc -i pipeline_architecture.mmd -o base.svg -p pptr.json     # pptr.json = {"args":["--no-sandbox"]}
python3 pipeline_postprocess.py base.svg pipeline_final.svg    # rotate phase labels + red numbered circles
# screenshot the SVG through a browser so <foreignObject> labels render:
printf '<!doctype html><meta charset="utf-8"><style>body{margin:0;background:#fff}</style>' > wrap.html
cat pipeline_final.svg >> wrap.html
chrome --headless=new --force-device-scale-factor=2 --window-size=1854,2776 \
       --screenshot=pipeline_architecture.png file://$PWD/wrap.html
```
