# Communication-Efficient Model Averaging Prediction for Massive Data with Asymptotic Optimality

This repository contains the code used for the simulation and empirical studies presented in the paper *"Communication-Efficient Model Averaging Prediction for Massive Data with Asymptotic Optimality."*

## Structure

### 1. Simulation Code (`simulation code for dmap`)

The folder `simulation code for dmap` includes the code for the simulations in the paper. It is divided into two subfolders, corresponding to **Example 1** and **Example 2** from the paper:

*   `sim DGP1`: Contains the code and results for Example 1.

*   `sim DGP2`: Contains the code and results for Example 2.

Each subfolder includes further subdirectories, which align with specific simulation steps in the paper:

*   **1 determine T and phi**: Code for determining the parameters (T) and (\phi).

*   **2 compare six method with fixed T and phi**: Comparison of six methods with fixed parameters (T) and (\phi).

*   **3 compare CPU time of DMAP-SA, DMAP-SL with gMAP**: Code for comparing CPU times of DMAP-SA, DMAP-SL, and gMAP.

*   **4 compare two initial values for DMAP-SL**: Comparison of two initial values for the DMAP-SL method.

### 2. Empirical Code (`real data code for dmap`)

The folder `real data code for dmap` contains the code used for the empirical analysis in the paper. It is organized into two subfolders:

*   `1 MSPEs and time for three types of candidate models`: Code for calculating the Mean Squared Prediction Errors (MSPEs) and computing time for for three sets of candidate models for
California housing prices dataset.

*   `2 Estimation of weighted coefficients`: Code for estimating the weighted coefficients.

## Usage

Each folder contains code and related data files for reproducing the results presented in the paper. Please refer to the specific folder corresponding to the example or empirical analysis you are interested in.
