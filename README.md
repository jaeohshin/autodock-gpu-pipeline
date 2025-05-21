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

