**Course:** Neural Signal Processing
**Authors:** Muralidhar Rao M, Prateek Kumar, Shreeansh Hota

---

## Project Overview
This project investigates the "Shape of Stillness" hypothesis: the transition from a resting state to focused meditation is not merely signal amplification, but a collapse of degrees of freedom. Rather than treating the brain as a static array of independent channels or a "Black Box" pattern generator, this research aims to show that focus leads to the emergence of a low-dimensional neural manifold with high temporal stability.

## Dataset
* **Source:** Biswas et al.
* **Data Type:** High-density 64-channel EEG.
* **Subjects:** 30 Rajyoga practitioners and 30 Controls.
* **Conditions:** EO1 (Baseline) versus M1 (Meditation).

## Methodology
The analysis framework is divided into three main tiers:

### 1. Signal Subspace Extraction
* Utilizes Principal Component Analysis (PCA) combined with Horn's Parallel Analysis.
* **Purpose:** Collapses correlated electrodes to handle volume conduction and filters out high-frequency noise that does not contribute to the cognitive state.

### 2. Quantitative Metrics Analysis
Calculates three dynamical markers within the k-dimensional signal subspace:
* **Trajectory Velocity:** Measures the Euclidean distance between consecutive time points to proxy volatility.
* **State Space Volume:** Measures the volume of the hyper-ellipsoid enclosing the trajectory to quantify dynamic range.
* **Excursion Radius:** Measures the mean distance of trajectory points from the geometric centroid to assess spatial spread.

### 3. Intrinsic Complexity & Visualization
* **Two-Nearest Neighbors (TNN) Estimator:** Used to estimate the true non-linear degrees of freedom ($d$), as linear methods like PCA can overestimate dimensionality in curved manifolds.
* **Mapping Attractor Topology:** Projects trajectories onto PC1–PC3 for visual inspection of stable loops versus chaotic wandering.

## Expected Outcomes
* **Metric Reduction:** Meditators are expected to show reduced Trajectory Velocity and State Space Volume, confirming the "Shape of Stillness" hypothesis.
* **Lower Intrinsic Complexity:** A significant decrease in the Intrinsic Dimension ($d$) via TNN compared to the control group.
* **High Reliability:** A lower standard deviation of $d$ within the meditation group, indicating that long-term practice drives neural dynamics toward a shared, stable complexity regime.