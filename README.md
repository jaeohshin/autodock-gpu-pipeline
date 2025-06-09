# AutoDock-GPU Docking Pipeline

This repository contains a fully automated virtual screening pipeline using [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU), [AutoGrid](https://github.com/ccsb-scripps/AutoGrid), [MGLTools](http://mgltools.scripps.edu/), and [Meeko](https://github.com/forlilab/Meeko). The pipeline enables batch docking of ligands against receptor ensembles using GPU acceleration.

---

## 🔧 Features

* ✅ Automated receptor preparation using `prepare_receptor4.py` (MGLTools)
* ✅ Ligand preparation via [Meeko](https://github.com/forlilab/Meeko), with fallback to `prepare_ligand4.py` (MGLTools)
* ✅ Grid map generation using `prepare_gpf4.py` + `autogrid4`
* ✅ Fast docking using `autodock_gpu_128wi` (CUDA)
* ✅ Customizable grid center and box size
* ✅ Supports batch ligand screening per target
* ✅ Compatible with PyMOL for visualizing poses

---

## 📁 Repository Structure

This modular structure supports scalable, resumable docking workflows, making it easy to manage thousands of ligands, protein conformations, and kinase targets.

```text
virtual_screening/
├── input/
│   ├── kinase_list.txt
│   ├── ligands/
│   │   └── abl1/
│   │       ├── actives/
│   │       └── decoys/
│   └── receptors/
│       └── abl1/
├── preprocessed/
│   ├── ligands_pdbqt/
│   ├── receptors_pdbqt/
│   ├── aligned_ligands/
│   └── grid_centers/
├── grids/
├── docking_output/
├── results/
├── logs/
├── scripts/
│   ├── run_vs.py
│   ├── align_ligands.py
│   └── prepare_grids.py
```

---

## ⚙️ Dependencies

* Python 3.x
* [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU)
* [AutoGrid](https://github.com/ccsb-scripps/AutoGrid)
* [MGLTools 1.5.7](http://mgltools.scripps.edu/)
* [RDKit](https://www.rdkit.org/)
* [Meeko](https://github.com/forlilab/Meeko)

---

## 🚀 Quickstart: Virtual Screening Workflow

Follow these steps to run the full docking pipeline:

1. **Prepare receptor structures**
   Generate receptor `.pdb` files (e.g., from BioEmu or other sources), then:

   ```bash
   bash scripts/prepare_dude_receptors.sh
   ```

2. **Prepare ligand files**
   Convert SDF ligands (actives/decoys) into `.pdbqt` using Meeko or MGLTools:

   ```bash
   bash scripts/prepare_dude_ligands.sh
   ```

3. **Generate grid maps**
   Compute grid centers (e.g., via PyMOL selection) and generate `.gpf` / `.fld` files:

   ```bash
   /usr/bin/python3 -m pymol -cq scripts/prepare_grids.py --kinase abl1
   ```

4. **Run virtual screening**
   Launch batch docking with AutoDock-GPU:

   ```bash
   python scripts/run_vs.py
   ```

Each step can be configured or modified via script arguments or by editing the corresponding files. Intermediate outputs (e.g., `.pdbqt`, `.fld`, `.dlg`) will be saved in the `preprocessed/`, `grids/`, and `docking_output/` directories.

---

## 📊 Receptor & Grid Preparation

Receptor preprocessing and grid generation are handled by MGLTools and AutoGrid to ensure compatibility with AutoDock's docking engine. The pipeline uses `.pdbqt`, `.gpf`, and `.fld` files for docking setup.

---

## 🤝 Acknowledgments

This pipeline integrates open-source tools developed by:

* The Scripps Research Institute: AutoDock-GPU, AutoGrid, MGLTools
* The Forli Lab: Meeko

