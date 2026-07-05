"""
Cover image for "A Chiral-Order Glycation Channel for Synthetic Cells".
A homochiral scaffold (twin helical strands) crisp and single-handed at top,
fraying into disorder as D-glucose glycates it downward; chiral order U_chi -> 0,
information I = 1 - U_chi rises; a CD-readout motif and the governing relations.
Portrait, sized to the manuscript cover aspect (6.5 x 8.0625 in).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, RegularPolygon, FancyBboxPatch
from matplotlib.collections import LineCollection
import matplotlib.cm as cm

rng = np.random.default_rng(11)
W, H = 6.5, 8.0625
fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0,0,1,1]); ax.axis("off")
ax.set_xlim(0,10); ax.set_ylim(0,12.4)

# ---- background: soft vertical gradient (ordered green-teal -> warm disorder) ----
grad = np.linspace(0,1,400).reshape(-1,1)
ax.imshow(grad, extent=[0,10,0,12.4], origin="lower", aspect="auto",
          cmap=plt.cm.get_cmap("BrBG"), alpha=0.14, zorder=0)
ax.add_patch(FancyBboxPatch((0.0,0.0),10,12.4,boxstyle="square,pad=0",
            fc="white",ec="none",alpha=0.35,zorder=0))

# ---- title block ----
ax.text(5.0, 11.7, "A Chiral-Order Glycation Channel", ha="center", va="center",
        fontsize=21, weight="bold", color="#14332b")
ax.text(5.0, 11.15, "for Synthetic Cells", ha="center", va="center",
        fontsize=21, weight="bold", color="#14332b")
ax.text(5.0, 10.62, "Engineering — and Resisting — Sugar-Driven Information Loss",
        ha="center", va="center", fontsize=11.5, style="italic", color="#5a3a2a")

# ---- hero graphic: degrading twin helix ----
cx = 5.0
t = np.linspace(0, 1, 900)                 # 0 top ordered, 1 bottom disordered
y = 9.7 - 7.2*t                            # descend
amp = 1.35
disorder = t**1.7                          # grows downward
def strand(phase0):
    theta = 2*np.pi*4.3*t + phase0 + disorder*rng.normal(0,3.0,t.size).cumsum()*0.012
    r = amp*(1 - 0.15*disorder) + disorder*rng.normal(0,0.6,t.size)*0.25
    x = cx + r*np.sin(theta)
    return x, y
# color along strand: ordered (deep teal) -> disordered (warm red)
cmap = cm.get_cmap("Spectral_r")
def draw_strand(phase0, lw=3.2):
    x, yy = strand(phase0)
    pts = np.array([x, yy]).T.reshape(-1,1,2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=cmap, array=t, linewidth=lw, alpha=0.95, zorder=4)
    ax.add_collection(lc)
    return x, yy
x1,_ = draw_strand(0.0)
x2,_ = draw_strand(np.pi)
# faint rungs (base-pair-like) where still ordered
for i in range(0, 520, 26):
    ax.plot([x1[i], x2[i]], [y[i], y[i]], color="#4f8f86",
            lw=1.0, alpha=max(0, 0.5-0.5*t[i]), zorder=3)

# ---- D-glucose molecules impinging where fraying begins ----
def glucose(px, py, s=0.42, col="#c85b3c"):
    hexg = RegularPolygon((px,py), 6, radius=s, orientation=np.pi/6,
                          fc="none", ec=col, lw=2.0, zorder=6)
    ax.add_patch(hexg)
    for a in np.linspace(0, 2*np.pi, 6, endpoint=False):
        ax.plot([px+s*np.cos(a), px+(s+0.16)*np.cos(a)],
                [py+s*np.sin(a), py+(s+0.16)*np.sin(a)], color=col, lw=1.6, zorder=6)
for gx,gy in [(2.15,6.4),(7.7,5.6),(2.6,4.5),(7.3,3.8),(1.9,3.0)]:
    glucose(gx,gy)
    ax.add_patch(FancyArrowPatch((gx,gy),(cx+ (0.9 if gx<cx else -0.9), gy),
                 arrowstyle="-|>", mutation_scale=13, lw=1.6, color="#c85b3c",
                 alpha=0.8, zorder=5, connectionstyle="arc3,rad=0.12"))

# ---- annotations: order at top, disorder at bottom ----
ax.text(cx, 10.05, r"homochiral order   $U_\chi \approx 1$", ha="center",
        fontsize=11.5, color="#14332b", weight="bold")
ax.annotate("", xy=(cx-2.55, 5.0), xytext=(cx-2.55, 9.2),
            arrowprops=dict(arrowstyle="-|>", color="#8a5a3a", lw=2.2))
ax.text(cx-2.9, 7.1, "glycation", rotation=90, ha="center", va="center",
        fontsize=12, color="#8a5a3a", weight="bold", style="italic")
ax.text(cx, 2.15, r"information loss   $I = 1 - U_\chi$", ha="center",
        fontsize=12.5, color="#8a2f22", weight="bold")

# ---- CD-readout motif (upper-right inset) ----
axcd = fig.add_axes([0.70, 0.66, 0.24, 0.11]); axcd.set_facecolor("none")
w = np.linspace(0,1,200)
sig = np.sin(2*np.pi*w)*np.exp(-1.5*w)      # bisignate CD curve
axcd.plot(w, sig, color="#1a6f63", lw=2.2)
axcd.plot(w, 0.25*sig, color="#c0392b", lw=1.6, ls=(0,(3,2)))
axcd.axhline(0, color="#999", lw=0.6)
axcd.set_xticks([]); axcd.set_yticks([])
for sp in axcd.spines.values(): sp.set_color("#bbbbbb")
axcd.set_title("circular dichroism readout", fontsize=7.5, color="#333")

# ---- equations strip (bottom) ----
ax.add_patch(FancyBboxPatch((0.7,0.95),8.6,0.95,boxstyle="round,pad=0.08,rounding_size=0.12",
            fc="#f3f0ea", ec="#cdbfa8", lw=1.2, zorder=7))
ax.text(5.0, 1.42, r"$\dfrac{dI}{dt}=k_{glyc}\,G\,(1-I)-R\,I$"
        r"$\qquad\rho=\dfrac{R}{k_{glyc}\,G}\qquad I^{*}=\dfrac{1}{1+\rho}$",
        ha="center", va="center", fontsize=13, color="#2a2a2a", zorder=8)

# ---- author ----
ax.text(5.0, 0.45, "Leon Sandler  ·  Independent Researcher, Northbrook, Illinois",
        ha="center", va="center", fontsize=10.5, color="#333333")

plt.savefig("/home/claude/cover_chiral.png", dpi=150, facecolor="white")
print("cover saved")
