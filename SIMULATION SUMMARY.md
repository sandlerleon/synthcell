# Simulation Summary

**Project:** Information-loss / chiral-order modeling for synthetic-cell aging
**Author:** Leon Sandler (Independent Researcher, Northbrook, IL)
**Scope:** All simulations supporting (a) the SpudCell information-maintenance argument added to *Aging as Information Loss*, and (b) the new manuscript *A Chiral-Order Glycation Channel for Synthetic Cells*.
**Status:** All parameters illustrative and unfitted. Simulations establish internal coherence, coarse-graining validity, and mechanism — not measured effect sizes.

---

## 1. SpudCell information-maintenance model
**File:** `spudcell_maintenance.py` → `spudcell_maintenance.png`
**Role:** Adapts the scalar information-loss skeleton of the aging paper to the *maintenance-positive* regime, to answer whether a minimal synthetic cell can persist on its own.

**Governing law (generation-time):**
```
dI/dt = (r0 + L)(1 - I) - R·I
I* = L/(L + R)          # steady-state disorder
sustainable iff  R > L·(1/I_c - 1)   # for I_c = 0.5:  R > L
```
- `I` = fractional information disorder of the cell's specification (∈[0,1])
- `L` = per-generation loss pressure; `R` = energy-funded repair; `I_c` = error-catastrophe threshold
- Repair `R = R_max · A/(K_A+A) · m` — `A` = ATP availability (light + fed), `m` = repair-machinery availability

**Key illustrative parameters:** genome 90 kbp / 7 plasmids; `eps_base=1.5e-6`, `p_seg=0.03`, `r0=0.01`, `R_max=0.60`, `K_A=0.5`, `I_c=0.5`.

**Key results:**

| Loss channel | Value (per generation) |
|---|---|
| Replication copy errors `L_rep` | 0.045 |
| **7-plasmid missegregation `L_seg`** | **0.192 (dominant)** |
| Total loss pressure `r0 + L` | 0.247 |
| Repair needed to survive (`I_c=0.5`) | R > 0.247 |

| Scenario | `I*` | Verdict |
|---|---|---|
| Current SpudCell (m=0, no self-repair) | 1.000 | collapses, crosses `I_c` at ~2 generations |
| Fed + partial machinery (m=0.35) | 0.638 | collapses at ~4 generations |
| Light-powered + autopoietic (m=0.9) | 0.379 | **sustainable** |

**Takeaways:** (i) the dominant loss channel is plasmid mis-partitioning, not copy errors — consolidating the 7 plasmids to 1 is the highest-leverage move (independently matches the SpudCell team's own stated priority); (ii) with no self-repair the lineage cannot persist; (iii) light-derived ATP only lowers `I*` once repair machinery exists (m>0) — energy is the fuel, not the mechanic.

---

## 2. Chiral-order glycation channel — phenomenological model (Figure 1)
**File:** `chiral_glycation.py` → `chiral_glycation.png`
**Role:** Defines the channel and its coarse-grained dynamics for the new manuscript.

**Governing law:**
```
I = 1 - U_chi                          # information loss = loss of chiral order
k_glyc = k0·(1 - s·u_s)                # surface uniformity shields glycation sites
dI/dt = k_glyc·G·(1 - I) - R·I         # glycation loss vs scaffold turnover
I* = k_glyc·G / (k_glyc·G + R)
```
- `U_chi` ∈ [0,1] = chiral order parameter (CD-readable); `G` = glucose exposure
- `u_s` = geometric surface uniformity (resistance knob); `R` = turnover/repair

**Key illustrative parameters:** `k0=0.9`, `s=0.75`, `R0=0.15`, `I0=0.05`.

**Key results:**
- `k_glyc(u_s=0) = 0.900`; `k_glyc(u_s=0.8) = 0.360` — surface uniformity confers ~2.5× resistance.
- Panel A: chiral order decays faster at higher glucose (engineered vulnerability).
- Panel B: `k_glyc` falls linearly with `u_s`; steady-state `U_chi` rises (inset).
- Panel C: all-D vs all-L scaffolds diverge under D-glucose (illustrative; sign/magnitude are predictions, not data).

---

## 3. Chiral-order glycation channel — mechanistic lattice Monte Carlo (Figure 2)
**File:** `chiral_montecarlo.py` → `chiral_montecarlo.png`
**Role:** Independent, non-circular support. Chiral order, resistance, and the enantioselective split are **outputs** of site-level rules, not imposed. Also tests whether the Box 1 ODE is a valid reduction.

**Microscopic rules (imposed):** lattice of oriented chiral units; per-site glycation prob `p0·G·exposure·stereo`; each glycation kicks orientation (packing defect); adjacent glycated sites impose cooperative crosslink strain (AGE crosslinks); surface uniformity `u_s` sets the exposure distribution (uniform ⇒ tight, occluded, low-variance); stereochemical barrier `ε` makes D-glucose favor the matched-handed scaffold; optional turnover `R` replaces glycated units with fresh homochiral material.
**Order parameter (measured):** `U_chi = |⟨exp(i·θ_i)⟩|` over the lattice.

**Key illustrative parameters:** lattice 50×50 (2500 sites), 220 steps, 6-seed average; `p0=0.020`, `kappa_defect=0.95`, `kappa_strain=0.55`, `p_x=0.16`, `eps=0.35`, panel-C turnover `R_C=0.045`.

**Key results:**

| Test | Result | Supports |
|---|---|---|
| **A** — ODE reduction | measured `k/G` = 0.0182, 0.0151, 0.0117 for G=0.4/1.0/2.2 → constant to within ~30% over a 5× range | Box 1 ODE is a faithful coarse-graining; residual sub-linearity = crosslink cooperativity (the anticipated super-linear correction) |
| **B** — resistance | emergent `k_glyc` falls 0.0154 → 0.0068 as `u_s` 0→1 (ratio **2.27**) | resistance knob emerges from steric occlusion, not assumed |
| **C** — chiral signature | steady-state `U_chi`: all-D 0.663, racemic 0.715, all-L 0.781 (D/L order-loss ratio **1.54**) | persistent enantioselective split from a single barrier `ε` (prediction 4) |
| **D** — order field | snapshots show transition from uniform single-handed order → accumulated, crosslink-propagated defects | visual mechanism |

**Honesty notes:** parameters are illustrative; the sign/magnitude of the ε-split remain predictions. The resistance result is conditional on the physically-motivated assumption that surface uniformity lowers per-site accessibility — the *emergence* is real, the micro-rule is a modeling choice (stated in the manuscript).

---

## Reproducibility

- **Dependencies:** Python 3, NumPy, SciPy, Matplotlib.
- **Determinism:** all stochastic runs use fixed seeds; Monte Carlo curves are 6-seed averages.
- **Runtime:** each script runs in seconds on a single core.
- **Run:** `python3 <script>.py` — regenerates the corresponding PNG and prints the numerical summary above.

## File manifest

| File | Type | Produces |
|---|---|---|
| `spudcell_maintenance.py` | model | `spudcell_maintenance.png` (SpudCell autonomy) |
| `chiral_glycation.py` | model | `chiral_glycation.png` (Figure 1) |
| `chiral_montecarlo.py` | model | `chiral_montecarlo.png` (Figure 2) |

*Prepared to accompany the Zenodo/GitHub code deposit for the two manuscripts. DOIs to be inserted on deposit.*

---

## 4. Chiral order ↔ configurational entropy bridge (Figure 2, v3)
**File:** `entropy_bridge.py` → `entropy_bridge.png`
**Role:** Makes the "chiral order = information" claim exact (addresses Reviewer concern 1).

For the maximum-entropy (von Mises) distribution of scaffold orientations with concentration κ:
```
U_chi(κ) = I1(κ)/I0(κ)                       # mean resultant length = chiral order
S(κ)     = ln(2π I0(κ)) - κ·U_chi(κ)         # configurational entropy (nats)
H_chi    = S_max - S,  S_max = ln 2π         # chiral information content (nats)
```
Monotonic map U_chi ↔ S. Reference points: U_chi=0.5 → H_chi≈0.27 nats; U_chi=0.95 → H_chi≈1.56 nats. Loss of chiral order is exactly loss of configurational information.

## 5. Sensitivity & robustness + dimensionless collapse (Figure 5, v3)
**File:** `sensitivity_analysis.py` → `sensitivity_analysis.png`
**Role:** Robustness to parameters (Reviewer improvement 2) and dimensionless generality (improvement 3).

- **A — dimensionless collapse:** 400 random (G, R, u_s) triples fall on `I* = 1/(1+ρ)`, `ρ = R/(k_glyc·G)`; sustainability boundary `ρ*=1` (I_c=0.5).
- **B — stereochemical barrier:** enantioselective split `U_L−U_D` grows monotonically with ε (0.003 at ε=0 → 0.229 at ε=0.6); vanishes in the achiral control.
- **C — lattice size:** decay is size-independent (L=30/50/80) down to the `1/√N` floor (0.033/0.020/0.012).
- **D — turnover vs glucose:** steady-state U_χ rises with R, falls with G; catastrophe line crossed at `ρ=ρ*`.

**Single controlling group:** `ρ = R/[k0·G·(1−s·u_s)]` — conclusions depend on the turnover-to-glycation ratio, not absolute rates.

*Note:* `fig_overview.py` → `fig_overview.png` is the "why this matters" schematic (Figure 1, v3), a diagram rather than a simulation.

## Updated file manifest (v3)

| File | Type | Produces |
|---|---|---|
| `spudcell_maintenance.py` | model | SpudCell autonomy |
| `chiral_glycation.py` | model | Figure 3 (phenomenological channel) |
| `chiral_montecarlo.py` | model | Figure 4 (lattice Monte Carlo) |
| `entropy_bridge.py` | calculation | Figure 2 (entropy–order bridge) |
| `sensitivity_analysis.py` | model | Figure 5 (robustness + dimensionless collapse) |
| `fig_overview.py` | schematic | Figure 1 (overview diagram) |

---

## 6. Statistical-mechanics (site-diluted XY) reformulation (Figure 5, v4)
**File:** `xy_statmech.py` → `xy_statmech.png`
**Role:** Recasts the lattice model as a proper Hamiltonian system (Reviewer Critique 2 + Strategy C); gives the orientation "kick" and the stereochemical barrier ε a physical meaning.

**Hamiltonian (site-diluted 2D XY + pinning field):**
```
H = - Σ_<ij> J0 g_i g_j cos(θ_i-θ_j) - Σ_i h0 g_i cos(θ_i)
```
- `g_i ∈ [0,1]` order integrity (1 intact → g_gl≪1 on glycation): glycation locally destroys XY coupling J and pinning field h.
- Glycation = kinetic Monte Carlo reaction at rate `p0·G·exposure·exp(±ε/2)`, `ε = ΔΔG‡/kBT` (activation-barrier difference).
- Orientations relax by checkerboard Metropolis sweeps (kBT=1); the "kick" is the equilibrium response to lost coupling.

**Key results:** (A) reduces to the same Box 1 ODE — k/G = 0.0189/0.0182/0.0174 across G=0.5/1.0/2.0 (~8% constant); (B) coupling loss drives an order–disorder transition vs glycated fraction; (C) enantioselective split scales with ε (−0.02 at ε=0 → +0.29 at ε=1.6).

## Parameter grounding (v4, §3.1)
Parameters anchored to measured chemistry: glucose glycation of poly-L-lysine/model peptides (days–weeks, exponential; CD ellipticity change on glycation), lysine site-reactivity set by steric availability, peptide-hydrogel turnover (days–weeks). Both k_glyc·G and R fall in the inverse-days-to-weeks range → ρ = R/(k_glyc·G) ~ O(1) is experimentally reachable; HbA1c (glycation vs RBC turnover) is the in vivo analogue. Simulations reframed as "predictive bounds for an achievable design," not toy models.

## Updated file manifest (v4)

| File | Type | Produces |
|---|---|---|
| `spudcell_maintenance.py` | model | SpudCell autonomy |
| `chiral_glycation.py` | model | Figure 3 (phenomenological channel) |
| `chiral_montecarlo.py` | model | Figure 4 (lattice Monte Carlo) |
| `xy_statmech.py` | model | Figure 5 (site-diluted XY / kMC) |
| `sensitivity_analysis.py` | model | Figure 6 (robustness + dimensionless collapse) |
| `entropy_bridge.py` | calculation | Figure 2 (entropy–order bridge) |
| `fig_overview.py` | schematic | Figure 1 (overview diagram) |
