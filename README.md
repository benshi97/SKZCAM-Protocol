[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/benshi97/SKZCAM-Protocol/HEAD?labpath=generate_clusters.ipynb)

### SKZCAM-Protocol
Scripts to implement the SKZCAM protocol for systematically generating quantum clusters of increasing size in the electrostatic embedded cluster approach. More details can be found in the original paper: **General embedded cluster protocol for accurate modeling of oxygen vacancies in metal-oxides** by Benjamin X. Shi, Venkat Kapil, Andrea Zen, Ji Chen, Ali Alavi and Angelos Michaelides ([arXiv 2202.04633](https://arxiv.org/abs/2202.04633)).

The notebook `generate_clusters.ipynb` documents the Python function (just a single function!!) that needs to be used to generate the clusters from just one input file. These clusters can be generated and visualized interactively with [Binder](https://mybinder.org/v2/gh/benshi97/SKZCAM-Protocol/HEAD?labpath=generate_clusters.ipynb). If there are any problems with the scripts, please feel free to email me at: [bxs21@cam.ac.uk](mailto:bxs21@cam.ac.uk).

The core dependencies required for this notebook to work are:
* Python 3
* [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
* Matplotlib
* Numpy

