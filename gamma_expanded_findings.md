# Gamma-Band Expanded State Space Analysis

**Band:** Gamma (30–80 Hz)
**Data folder:** `savedDataGamma_AllStates`
**Protocols:** EO1 (eyes-open baseline), G1 (gamma-induction protocol), M2 (meditation, second session)
**Pipeline:** Phase Randomization surrogates, within-trial velocity, 1e-6 eigenvalue clamp, PR, ANCOVA
**K_std:** 5
**Significance markers:** ★ p<0.05, ★★ p<0.01, ★★★ p<0.005

---

## 1. Between-Group at M2: Meditators vs Controls

| Metric | Meditators (mean +/- SEM) | Controls (mean +/- SEM) | p-value |
|---|---|---|---|
| k | 6.13 +/- 0.49 | 4.87 +/- 0.47 | 0.1215 |
| **Velocity** | **30.53 +/- 2.44** | **22.30 +/- 1.52** | **0.0047 ★★★** |
| **Log Volume** | **13.57 +/- 0.40** | **11.54 +/- 0.36** | **0.0008 ★★★** |
| **Radius** | **25.79 +/- 2.07** | **18.90 +/- 1.28** | **0.0050 ★★★** |
| PR | 5.86 +/- 0.43 | 4.75 +/- 0.41 | 0.0571 |
| TNN d | 17.34 +/- 0.26 | 17.50 +/- 0.45 | 0.7813 |

## 2. ANCOVA at M2: LogVolume ~ Group + TotalVariance

| Predictor | Coefficient | SE | p-value |
|---|---|---|---|
| **Group (Med vs Con)** | **+1.5283** | 0.4704 | **0.0019 ★★★** |
| TotalVariance | +0.00039 | 0.00008 | 1.18×10⁻⁵ ★★★ |
| (R² = 0.427, n = 30 per group) | | | |

---

## 3. Within-Group Wilcoxon Signed-Rank: All Three Transitions

### 3.1 EO1 vs G1

| Metric | Group | EO1 (mean +/- SEM) | G1 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 6.41 +/- 0.50 | 6.37 +/- 0.47 | 0.6651 |
| k | Controls | 5.27 +/- 0.47 | 5.38 +/- 0.43 | 0.3338 |
| Velocity | Meditators | 25.27 +/- 1.96 | 24.41 +/- 1.66 | 0.3695 |
| Velocity | Controls | 20.41 +/- 1.78 | 19.13 +/- 1.59 | 0.3147 |
| Log Volume | Meditators | 12.20 +/- 0.38 | 12.33 +/- 0.34 | 0.8372 |
| Log Volume | Controls | 10.93 +/- 0.33 | 10.74 +/- 0.32 | 0.6114 |
| Radius | Meditators | 21.36 +/- 1.62 | 20.68 +/- 1.36 | 0.4051 |
| Radius | Controls | 17.23 +/- 1.50 | 16.14 +/- 1.33 | 0.3044 |
| PR | Meditators | 5.81 +/- 0.41 | 6.16 +/- 0.43 | 0.9914 |
| PR | Controls | 5.27 +/- 0.42 | 5.33 +/- 0.44 | 0.2748 |
| TNN d | Meditators | 18.31 +/- 0.39 | 18.29 +/- 0.41 | 0.7869 |
| TNN d | Controls | 18.64 +/- 0.53 | 18.90 +/- 0.43 | 0.5666 |

### 3.2 G1 vs M2

| Metric | Group | G1 (mean +/- SEM) | M2 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 6.37 +/- 0.46 | 6.13 +/- 0.49 | 0.6229 |
| k | Controls | 5.38 +/- 0.43 | 4.87 +/- 0.48 | 0.2182 |
| Velocity | Meditators | 24.41 +/- 1.64 | 30.53 +/- 2.44 | 0.0207 ★ |
| Velocity | Controls | 19.13 +/- 1.59 | 22.30 +/- 1.54 | 0.0201 ★ |
| Log Volume | Meditators | 12.33 +/- 0.34 | 13.57 +/- 0.40 | **0.0020 ★★★** |
| Log Volume | Controls | 10.74 +/- 0.32 | 11.54 +/- 0.37 | 0.0125 ★ |
| Radius | Meditators | 20.68 +/- 1.34 | 25.79 +/- 2.07 | 0.0243 ★ |
| Radius | Controls | 16.14 +/- 1.33 | 18.90 +/- 1.30 | 0.0190 ★ |
| PR | Meditators | 6.16 +/- 0.42 | 5.86 +/- 0.43 | 0.7036 |
| PR | Controls | 5.33 +/- 0.44 | 4.75 +/- 0.41 | 0.2218 |
| TNN d | Meditators | 18.29 +/- 0.40 | 17.34 +/- 0.26 | 0.0140 ★ |
| TNN d | Controls | 18.90 +/- 0.43 | 17.50 +/- 0.46 | **0.0011 ★★★** |

### 3.3 EO1 vs M2

| Metric | Group | EO1 (mean +/- SEM) | M2 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 6.41 +/- 0.50 | 6.13 +/- 0.49 | 0.4913 |
| k | Controls | 5.27 +/- 0.46 | 4.87 +/- 0.47 | 0.4276 |
| Velocity | Meditators | 25.27 +/- 1.96 | 30.53 +/- 2.49 | 0.0369 ★ |
| Velocity | Controls | 20.41 +/- 1.75 | 22.30 +/- 1.52 | 0.0786 |
| Log Volume | Meditators | 12.20 +/- 0.38 | 13.57 +/- 0.40 | **0.0008 ★★★** |
| Log Volume | Controls | 10.93 +/- 0.32 | 11.54 +/- 0.36 | 0.0285 ★ |
| Radius | Meditators | 21.36 +/- 1.62 | 25.79 +/- 2.11 | 0.0314 ★ |
| Radius | Controls | 17.23 +/- 1.48 | 18.90 +/- 1.28 | 0.0627 |
| PR | Meditators | 5.81 +/- 0.41 | 5.86 +/- 0.44 | 0.9569 |
| PR | Controls | 5.27 +/- 0.41 | 4.75 +/- 0.41 | 0.2895 |
| TNN d | Meditators | 18.31 +/- 0.39 | 17.34 +/- 0.26 | 0.0169 ★ |
| TNN d | Controls | 18.64 +/- 0.52 | 17.50 +/- 0.45 | **0.0060 ★★** |
