# Gamma-Band State Space Analysis: Findings

**Project:** Characterizing the Neural State Space Dynamics of Long-Term Rajyoga Meditation
**Authors:** Muralidhar Rao M, Prateek Kumar, Shreeansh Hota
**Date:** April 2026
**Analysis Band:** Gamma (30-80 Hz bandpass, downsampled to 250 Hz)
**Dataset:** Biswas et al. — 30 Rajyoga meditator-control pairs, 64-channel EEG
**Conditions:** EO1 (eyes-open baseline, 5 min) vs M1 (eyes-open meditation, 15 min)

---

## Analysis Pipeline

1. **Preprocessing:** Raw segmented EEG trials were bandpass filtered (30-80 Hz, 4th-order Butterworth, zero-phase), downsampled to 250 Hz, and concatenated across good trials to form continuous electrode-time trajectories. Bad trials were rejected using a three-layer pipeline (eye position, RMS amplitude, PSD thresholding) inherited from the Biswas et al. preprocessing code. Minimum thresholds: 30 good trials and 10 good electrodes per subject-protocol.

2. **Signal Subspace Extraction (Tier 1):** PCA was applied to each subject's electrode-time matrix. Horn's Parallel Analysis (100 permutations, 95th percentile threshold) determined the individual embedding dimension k. A global standardized dimension K_std was set as the median k of Controls at EO1 baseline to ensure all subsequent metrics were computed in the same dimensional space.

3. **Dynamical Metrics (Tier 2):** Three markers were computed in the K_std-dimensional PCA space:
   - Trajectory Velocity: mean Euclidean step size between consecutive time points
   - State Space Volume: log hyper-ellipsoid volume from eigenvalues of the projected covariance
   - Excursion Radius: mean distance of trajectory points from the geometric centroid

4. **Intrinsic Complexity (Tier 3):** The Two-Nearest Neighbors (TNN) estimator (Facco et al., 2017) was applied to the full-dimensional data (not k-truncated) to estimate nonlinear intrinsic dimension d. Data was subsampled to 10,000 points via uniform random sampling when trajectory length exceeded this threshold.

5. **Statistics:** All between-group comparisons used paired Wilcoxon signed-rank tests (matched meditator-control pairs). Within-group EO1 vs M1 comparisons also used paired Wilcoxon signed-rank. Significance thresholds: * p<0.05, ** p<0.01, *** p<0.005.

---

## Key Parameters

| Parameter | Value |
|---|---|
| Frequency band | 30-80 Hz (gamma) |
| Filter | 4th-order Butterworth bandpass, zero-phase (filtfilt) |
| Downsampled Fs | 250 Hz |
| Trial handling | Concatenation (not averaging — spontaneous gamma is not phase-locked) |
| Parallel analysis permutations | 100 |
| K_std | 6 (median of Control-EO1 k values) |
| Control-EO1 k range | [1, 10], mean = 5.5 |
| TNN subsampling | 10,000 points (uniform random) |
| Valid pairs (EO1) | 29/30 (053DR missing EO1) |
| Valid pairs (M1) | 29/30 (019CKa had only 5 good trials) |

---

## Results

### Between-Group Comparison: Meditators vs Controls (Paired Wilcoxon)

| Metric | Condition | Meditators (mean +/- SEM) | Controls (mean +/- SEM) | p-value | Direction |
|---|---|---|---|---|---|
| k | EO1 | 6.62 +/- 0.51 | 5.59 +/- 0.46 | 0.165 | Med > Con (N.S.) |
| k | M1 | 6.62 +/- 0.48 | 6.10 +/- 0.44 | 0.490 | Med > Con (N.S.) |
| Velocity | EO1 | 26.21 +/- 2.04 | 20.86 +/- 1.83 | 0.061 | Med > Con (trend) |
| **Velocity** | **M1** | **30.06 +/- 2.23** | **21.87 +/- 1.62** | **0.007 ★★** | **Med > Con** |
| Log Volume | EO1 | 13.76 +/- 0.46 | 12.09 +/- 0.40 | 0.013 ★ | Med > Con |
| **Log Volume** | **M1** | **14.97 +/- 0.46** | **13.02 +/- 0.46** | **0.007 ★★** | **Med > Con** |
| Radius | EO1 | 22.11 +/- 1.68 | 17.59 +/- 1.54 | 0.048 ★ | Med > Con |
| **Radius** | **M1** | **25.23 +/- 1.87** | **18.46 +/- 1.35** | **0.009 ★★** | **Med > Con** |
| TNN d | EO1 | 18.36 +/- 0.38 | 18.51 +/- 0.53 | 0.871 | No difference |
| TNN d | M1 | 17.81 +/- 0.34 | 17.79 +/- 0.41 | 0.888 | No difference |

### Within-Group Comparison: EO1 vs M1 (Paired Wilcoxon)

| Metric | Group | EO1 (mean +/- SEM) | M1 (mean +/- SEM) | p-value | Direction |
|---|---|---|---|---|---|
| k | Meditators | 6.64 +/- 0.53 | 6.71 +/- 0.49 | 0.938 | No change |
| k | Controls | 5.47 +/- 0.46 | 6.03 +/- 0.43 | 0.147 | No change |
| Velocity | Meditators | 26.60 +/- 2.08 | 29.40 +/- 2.21 | 0.080 | Trending up |
| Velocity | Controls | 20.98 +/- 1.77 | 21.67 +/- 1.58 | 0.393 | No change |
| **Log Volume** | **Meditators** | **13.83 +/- 0.47** | **14.83 +/- 0.46** | **0.004 ★★** | **Increases** |
| Log Volume | Controls | 12.11 +/- 0.38 | 12.96 +/- 0.45 | 0.015 ★ | Increases |
| Radius | Meditators | 22.40 +/- 1.72 | 24.62 +/- 1.83 | 0.092 | Trending up |
| Radius | Controls | 17.69 +/- 1.49 | 18.29 +/- 1.31 | 0.393 | No change |
| TNN d | Meditators | 18.32 +/- 0.39 | 17.79 +/- 0.35 | 0.328 | No change |
| TNN d | Controls | 18.58 +/- 0.52 | 17.90 +/- 0.41 | 0.037 ★ | Decreases |

---

## Hypothesis Evaluation

The "Shape of Stillness" hypothesis predicted that meditation would induce a collapse of neural degrees of freedom: lower trajectory velocity, smaller state space volume, reduced excursion radius, and lower intrinsic dimensionality in meditators compared to controls, particularly during the M1 meditation condition.

### Predictions vs Observations

| # | Prediction | Observed | p-value | Verdict |
|---|---|---|---|---|
| 1 | Meditators show lower velocity at M1 | Meditators show significantly **higher** velocity | 0.007 | **Rejected** |
| 2 | Meditators show lower state space volume at M1 | Meditators show significantly **higher** volume | 0.007 | **Rejected** |
| 3 | Meditators show lower excursion radius at M1 | Meditators show significantly **higher** radius | 0.009 | **Rejected** |
| 4 | Meditators show lower TNN intrinsic dimension at M1 | No difference between groups | 0.888 | **Not supported** |
| 5 | Velocity decreases from EO1 to M1 in meditators | Trends upward, not downward | 0.080 | **Rejected** |
| 6 | TNN d decreases from EO1 to M1 in meditators | No significant change | 0.328 | **Not supported** |
| 7 | Lower within-group variance of TNN d in meditators | Both groups show similar variance | — | **Not supported** |

**Summary: The "Shape of Stillness" hypothesis is not supported in the gamma band. The dynamical metrics show a statistically robust pattern in the opposite direction — meditators exhibit larger, faster, and more expansive gamma-band state space trajectories than controls, and this difference is amplified during meditation.**

---

## Interpretation

1. **Meditators have more dynamic gamma activity, not less.** All three dynamical metrics (velocity, volume, radius) are significantly elevated in meditators at M1 (p < 0.01). This is consistent with the Biswas et al. finding that meditation produces broadband gamma power increases across the scalp, and that meditators show higher spontaneous gamma even at rest (EO1).

2. **The group difference is already present at baseline.** Velocity, volume, and radius all show Med > Con trends at EO1 (p = 0.013 to 0.061), suggesting this is a trait-level difference in experienced practitioners, not purely a state effect of meditating. Meditation amplifies an already-existing difference.

3. **State space volume increases during meditation in both groups.** The EO1-to-M1 increase in log volume is significant for both meditators (p = 0.004) and controls (p = 0.015). This suggests that the act of sitting with eyes open in a meditation-like posture expands the gamma state space regardless of expertise, but meditators expand more.

4. **TNN intrinsic dimension shows no group effect.** The nonlinear intrinsic dimensionality is virtually identical between meditators and controls (~18 in both groups, both conditions). This implies that the manifold's intrinsic geometry is the same — the groups differ in the amplitude and speed with which they traverse that manifold, not in its topological complexity.

5. **The "collapse" hypothesis may apply to other frequency bands.** The literature motivating the hypothesis (Yu et al., 2022) focused on alpha-band (8-15 Hz) manifolds. Gamma-band dynamics may behave fundamentally differently: gamma reflects active cortical computation, so more focused meditation may produce more organized yet larger-amplitude gamma patterns rather than a suppression of activity. The hypothesis could still hold for alpha or broadband analyses.

---

## Limitations

- Parallel analysis used 100 permutations (reduced from 1000 for computational feasibility); k estimates may have slightly higher variance than with 1000 permutations.
- TNN subsampling to 10,000 points introduces stochasticity; results may vary slightly across runs.
- The analysis is restricted to the gamma band (30-80 Hz). The hypothesis may be better suited to lower frequency bands or broadband signals.
- Two subject-protocol combinations were excluded (019CKa M1: insufficient trials; 053DR EO1: missing data).
