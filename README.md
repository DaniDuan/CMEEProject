# My CMEE Project Repository

## Description
This repository contains scripts for simulaton and analysis on **The Role of Metabolic Strategies in Determining Microbial Community Richness along Temperature Gradients**. 
## Languages
Python, Jupyter Notebook, R, bash
## Dependencies
NumPy, Matplotlib, SciPy, LaTex

## Installation
To use scripts in this repository, clone the repository and run.

```bash
git clone https://github.com/DaniDuan/CMEEProject.git
```

## Project structure and Usage 
### Code directory
- **temp_microbial_CUE.ipynb:** A demonstration of the trajectory of thought flows on simulation, analysis and progress throughout the project.
- **Bacteria_vector_modular.py:** Main script for the simulation of resource uptake and growth of microbial communities.
- **model_func.py:** ODEs for integration.
- **parameters&#46;py:** Returning matrices and vectors of uptake, respiration and excretion.
- **size_temp_funcs.py:** TPCs for uptake and respiration rates.
- **simrun&#46;py:** Plotting script for figures used in thesis.
- **TPC_processing.R:** A script for plotting the growth rate TPCs using data from Tom Smith's empirical observations.
### Thesis directory
- **thesis.tex:** The Latex file for thesis. 
- **CompileLaTeX&#46;sh:** The bash file for compiling thesis.
## Author name and contact

Name: Danica Duan

Email: d.duan20@imperial.ac.uk