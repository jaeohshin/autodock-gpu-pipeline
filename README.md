# AutoDock-GPU Docking Pipeline

This repository contains a fully automated docking pipeline using [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU), [AutoGrid](https://github.com/ccsb-scripps/AutoGrid), [MGLTools](http://mgltools.scripps.edu/), and [Meeko](https://github.com/forlilab/Meeko), scripted in Python. It is designed to batch-process ligand–receptor docking jobs with high-performance GPU acceleration.

---

## 🔧 Features

- ✅ Automated receptor preparation using `prepare_receptor4.py` (MGLTools)
- ✅ Ligand preparation using [**Meeko**](https://github.com/forlilab/Meeko) with fallback to `prepare_ligand4.py` (MGLTools)
- ✅ Grid map generation using [`AutoGrid`](https://github.com/ccsb-scripps/AutoGrid) (`prepare_gpf4.py` and `autogrid4`)
- ✅ Docking using `autodock_gpu_128wi` (CUDA-enabled)
- ✅ Customizable grid center and box size
- ✅ Multi-ligand batch support
- ✅ Compatible with PyMOL for pose visualization

---

## 📁 Repository Structure


## ⚙️ Dependencies

- Python 3.x
- [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU)
- [AutoGrid](https://github.com/ccsb-scripps/AutoGrid)
- [MGLTools 1.5.7](http://mgltools.scripps.edu/)
- [RDKit](https://www.rdkit.org/)
- [Meeko](https://github.com/forlilab/Meeko)

### ✅ Installation

To install Meeko and RDKit:


conda install -c conda-forge meeko rdkit


<<<<<<< HEAD
### Ligands are first processed using Meeko. If Meeko fails (e.g., due to malformed PDB input), the pipeline automatically falls back to MGLTools' prepare_ligand4.py.
=======
🧠 Notes on Ligand Preparation
>>>>>>> 5c4d8a1 (update)

    Ligands are first processed using Meeko for modern chemical handling.

    If Meeko fails (e.g., due to malformed or minimal .pdb input), the pipeline automatically falls back to MGLTools' prepare_ligand4.py.

<<<<<<< HEAD

### This pipeline integrates tools developed by:
=======
    This ensures robustness across various ligand formats and chemical quality levels.

🧪 Receptor & Grid Preparation
>>>>>>> 5c4d8a1 (update)

    Receptor preparation and grid map generation are still handled by MGLTools and AutoGrid, which provide reliable .pdbqt, .gpf, and .fld generation compatible with AutoDock.

<<<<<<< HEAD
    The Forli Lab (Meeko)
=======
🤝 Acknowledgments

This pipeline integrates tools developed by:

    The Scripps Research Institute: AutoDock-GPU, AutoGrid, MGLTools

    The Forli Lab: Meeko
>>>>>>> 5c4d8a1 (update)
