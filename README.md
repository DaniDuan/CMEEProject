# My CMEE Project Repository

## Description
This repository contains python scripts for simulaton and analysis on the impact of temperature on bacterial communities. 
## Languages
Python, Jupyter Notebook
## Dependencies
NumPy, Matplotlib, SciPy

## Installation
To use scripts in this repository, clone the repository and run.

```bash
git clone https://github.com/DaniDuan/CMEEProject.git
```

## Project structure and Usage 
- **temp_microbial_CUE.ipynb:** A demonstration of the simulation, analysis and current progress on the project. 
- **Bacteria_vector_modular.py:** Main script for the simulation of resource uptake and growth of microbial communities.
- **model_func.py:** ODEs for integration.
- **parameters&#46;py:** Returning matrices and vectors of uptake, respiration and excretion.
- **size_temp_funcs.py:** Size and temperature dependence for uptake and repiration.
- **plots&#46;py:** Plotting time series of resource, species biomass and community richness.
- **rich_temp.py:** Returning community richness at different temperatures.
- **comp&#46;py:** Calculating the actual species resource uptake accounting competition, and returning the survival rate of species with the actual uptake values at different temperatures. 
## Author name and contact

Name: Danica Duan

Email: d.duan20@imperial.ac.uk