# CPP-DIP: Multi-objective Coverage Path Planning for MAVs in Dispersed and Irregular Plantations

## Overview
CPP-DIP is a multi-objective coverage path planning framework designed for Micro Air Vehicles (MAVs) to efficiently navigate dispersed and irregular plantations. It integrates image-based tree detection, density-aware waypoint generation, target-based path planning, and object-based path optimization to minimize flight distance, turning angles, intersection paths, and redundant coverage.

Our paper is available on arXiv [[PDF]]().

---

## Repository Structure
This repository contains the core algorithms evaluated in the paper:

- **Greedy Heuristic Insertion.m**  
  Implementation of the Greedy Heuristic Insertion (GHI) algorithm for CPP-DIP framework.

- **Ant Colony Optimization.m**  
  Implementation of the Ant Colony Optimization (ACO) algorithm for CPP-DIP framework.

- **Monte Carlo Reinforcement Learning.m**  
  Implementation of the Monte Carlo Reinforcement Learning (MCRL) method for CPP-DIP framework.

- **Clustering_KDE_DBSCAN.m**  
  Script for waypoint clustering and density-aware tree grouping using KDE and DBSCAN.

- **palm_tree_coordinates.txt**  
  Example plantation coordinates used as input data.

---

## Sample Data
To ensure reproducibility, we provide:

- **Sample images**: Oil palm dataset: https://drive.google.com/file/d/1pqzwsjBEopTnlbHHbPO2hpJrrhZ_qyu1/view?usp=sharing 
- **Example waypoint files**: `palm_tree_coordinates.txt`  

## Running the Algorithms

To reproduce the path planning results, open MATLAB and run the algorithm corresponding to your chosen method.  

- **Greedy Heuristic Insertion (GHI)**: run `Greedy Heuristic Insertion.m`  
- **Ant Colony Optimization (ACO)**: run `Ant Colony Optimization.m`  
- **Monte Carlo Reinforcement Learning (MCRL)**: run `Monte Carlo Reinforcement Learning.m`  

Before running any algorithm, preprocess the input tree coordinates using `Clustering_KDE_DBSCAN.m`.  

The provided example input file is `palm_tree_coordinates.txt`. Running the scripts will generate the optimized coverage path, which is saved either as a MATLAB figure or a `.mat` file for further analysis.
