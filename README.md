# AutoDock-GPU Docking Pipeline

This repository contains a fully automated docking pipeline using [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU), [AutoGrid](https://github.com/ccsb-scripps/AutoGrid), and [MGLTools](http://mgltools.scripps.edu/), scripted in Python. It is designed to batch-process ligand–receptor docking jobs with high-performance GPU acceleration.

---

## 🔧 Features

- ✅ Automated receptor and ligand preparation (`prepare_receptor4.py`, `prepare_ligand4.py`)
- ✅ Grid map generation using [`AutoGrid`](https://github.com/ccsb-scripps/AutoGrid) (`prepare_gpf4.py` and `autogrid4`)
- ✅ Docking using `autodock_gpu_128wi` (CUDA-enabled)
- ✅ Customizable grid center and box size
- ✅ Multi-ligand batch support
- ✅ Compatible with PyMOL for pose visualization

---

## 📁 Repository Structure

