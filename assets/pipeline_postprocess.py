import re, sys
svg=open(sys.argv[1]).read()

# ---------- 1. rotate phase labels to vertical, left-edge, centre-aligned ----------
def attr(blk,name):
    return float(re.search(rf'\b{name}="(-?[\d.]+)"',blk).group(1))

for pid in ("PHASE1","PHASE2","PHASE3"):
    blk=re.search(rf'<g class="cluster" id="my-svg-{pid}".*?</g>\s*</g>', svg, re.S).group(0)
    rect=re.search(r'<rect[^>]*/>',blk).group(0)
    x,y,w,h=attr(rect,'x'),attr(rect,'y'),attr(rect,'width'),attr(rect,'height')
    txt=re.sub(r'\s+',' ',re.search(r'<p>(.*?)</p>',blk,re.S).group(1)).strip()
    txt=re.sub(r'^[^\w]+','',txt).strip()   # drop leading emoji
    cm=re.search(r'color:(#[0-9a-fA-F]+)',blk); col=cm.group(1) if cm else '#e2e8f0'
    # drop original top label group
    blk2=re.sub(r'<g class="cluster-label".*?</g>\s*(?=</g>$)','',blk,flags=re.S)
    # fit font so rotated text stays within box height
    fs=min(17.0,(h-16)/max(1,len(txt))/0.5)
    fs=max(12.5,fs)
    cx=x+fs*0.72+5; cy=y+h/2
    label=(f'<text x="{cx:.1f}" y="{cy:.1f}" transform="rotate(-90 {cx:.1f} {cy:.1f})" '
           f'text-anchor="middle" dominant-baseline="central" '
           f'font-family="trebuchet ms,verdana,arial,sans-serif" font-size="{fs:.1f}" '
           f'font-weight="700" fill="#ffffff" letter-spacing="0.4">{txt}</text>')
    blk2=blk2[:-4]+label+'</g>'   # insert before final </g>
    svg=svg.replace(blk,blk2)

# ---------- 2. red numbered circles on each step node ----------
nums={'M1':'01','M2':'02','M3':'03','M4':'04','M5':'05','M6':'06','M7':'07'}
badges=[]
for node,nn in nums.items():
    m=re.search(rf'id="my-svg-flowchart-{node}-\d+"[^>]*transform="translate\(([\d.]+),\s*([\d.]+)\)"',svg)
    cx,cy=float(m.group(1)),float(m.group(2))
    r=re.search(r'<rect class="basic label-container"[^>]*x="(-?[\d.]+)" y="(-?[\d.]+)"',svg[m.end():m.end()+400])
    rx,ry=float(r.group(1)),float(r.group(2))
    bx,by=cx+rx+15,cy+ry+15   # top-left corner of node
    badges.append(
      f'<circle cx="{bx:.1f}" cy="{by:.1f}" r="14" fill="#dc2626" stroke="#ffffff" stroke-width="2"/>'
      f'<text x="{bx:.1f}" y="{by:.1f}" text-anchor="middle" dominant-baseline="central" '
      f'font-family="trebuchet ms,verdana,arial,sans-serif" font-size="14" font-weight="700" '
      f'fill="#ffffff">{nn}</text>')
svg=svg.replace('</svg>', ''.join(badges)+'</svg>')

open(sys.argv[2],'w').write(svg)
print("post-processed ->", sys.argv[2])
