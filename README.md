# mushroom_harvesters

Code for paper "Market integration and cooperative resource harvesting" 

## R code

- **descriptive_analysis.R** takes the network data and generate plots and summary tables, and formats/exports data for multitensor
- **regressions_hurdle.R** contains the stan regression models, diagnostics, tables and plots. 
- **multitensor_results.R** reads in the csv files with the multitensor results and generates plots because I am too lazy to do it in Python

## Python code

- **yunnan_MT.py** runs multitensor for various layer combinations and exports to csv files

Multitensor was developed by Caterina de Bacco and colleagues at the MPI for Intelligent Systems. Details of the multitensor algorithm and instructions for installation can be found at https://github.com/MPI-IS/multitensor.

## Data availability

