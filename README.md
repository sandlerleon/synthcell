# Sugar-Driven Information Loss in Minimal Synthetic Cells as a Model for Metabolic Aging

This repository contains the theoretical framework, mathematical models, and simulation scripts for investigating metabolic aging inside bottom-up minimal synthetic cells. 

Natural metabolic aging is highly confounded by organismal physiology, but this work engineers a standalone, fully controlled physical analog. We model a homochiral molecular scaffold embedded within a synthetic cell chassis where the geometric order parameter ($U_\chi$) serves as the information carrier. Non-enzymatic glycation by D-glucose acts as the defect-forming, stereoselective loss channel, which destroys configurational information in a quantifiable Shannon-Boltzmann sense.

This work serves as the direct experimental complement to the companion framework **"Aging as Information Loss: A Unified Dynamical Framework for Biological Aging."**

## 🔬 Repository Architecture

The codebase contains self-contained Python scripts that regenerate all figures and analytical boundaries presented in the manuscript:

* **Governing ODE Model**: Coarse-grained dynamical system balancing glycation flux against scaffold turnover.
* **Lattice Monte Carlo Simulation**: Code verifying emergent structural decay from local steric occlusion rules ($50\times50$ lattice setups).
* **Statistical Mechanics Model**: A site-diluted XY model with a pinning field mapping the loss of chiral order to an equilibrium thermal phase transition.
* **Sensitivity Analysis**: Scripts confirming parameter robustness across experimentally realistic timescales ($400$ randomly sampled parameter triples).

## 🔗 Repository Links
* **Zenodo Deposit:** https://zenodo.org

https://zenodo.org/uploads/21210277
