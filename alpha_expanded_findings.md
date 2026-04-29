# Alpha-Band Expanded State Space Analysis

**Band:** Alpha (8–13 Hz)
**Data folder:** `savedDataAlpha_AllStates`
**Protocols:** EO1 (eyes-open baseline), EC1 (eyes-closed), M1 (meditation, eyes-open)
**Pipeline:** Phase Randomization surrogates, within-trial velocity, 1e-6 eigenvalue clamp, PR, ANCOVA
**K_std:** 3
**Significance markers:** ★ p<0.05, ★★ p<0.01, ★★★ p<0.005

---

## 1. Between-Group at M1: Meditators vs Controls

| Metric | Meditators (mean +/- SEM) | Controls (mean +/- SEM) | p-value |
|---|---|---|---|
| k | 3.83 +/- 0.22 | 3.47 +/- 0.18 | 0.2306 |
| Velocity | 21.49 +/- 2.32 | 22.63 +/- 1.96 | 0.5521 |
| Log Volume | 8.32 +/- 0.27 | 8.65 +/- 0.29 | 0.4300 |
| Radius | 17.74 +/- 1.92 | 18.69 +/- 1.64 | 0.5235 |
| PR | 2.61 +/- 0.11 | 2.49 +/- 0.15 | 0.4051 |
| **TNN d** | **12.11 +/- 0.31** | **10.86 +/- 0.34** | **0.0479 ★** |

## 2. ANCOVA at M1: LogVolume ~ Group + TotalVariance

| Predictor | Coefficient | SE | p-value |
|---|---|---|---|
| Group (Med vs Con) | −0.1071 | 0.3196 | 0.7388 |
| TotalVariance | +0.00049 | 0.00008 | 1.59×10⁻⁷ ★★★ |
| (R² = 0.403, n = 29 per group) | | | |

---

## 3. Within-Group Wilcoxon Signed-Rank: All Three Transitions

### 3.1 EO1 vs EC1

| Metric | Group | EO1 (mean +/- SEM) | EC1 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 3.52 +/- 0.15 | 3.10 +/- 0.16 | 0.0220 ★ |
| k | Controls | 3.30 +/- 0.11 | 3.10 +/- 0.15 | 0.1758 |
| Velocity | Meditators | 20.44 +/- 2.24 | 24.81 +/- 2.34 | **0.0006 ★★★** |
| Velocity | Controls | 19.64 +/- 1.68 | 24.19 +/- 1.84 | **0.0029 ★★★** |
| Log Volume | Meditators | 7.98 +/- 0.28 | 8.68 +/- 0.28 | **0.0001 ★★★** |
| Log Volume | Controls | 8.14 +/- 0.25 | 8.68 +/- 0.25 | **0.0027 ★★★** |
| Radius | Meditators | 16.85 +/- 1.88 | 20.69 +/- 2.01 | **0.0005 ★★★** |
| Radius | Controls | 16.08 +/- 1.37 | 20.00 +/- 1.54 | **0.0027 ★★★** |
| PR | Meditators | 2.68 +/- 0.14 | 2.29 +/- 0.13 | **0.0024 ★★★** |
| PR | Controls | 2.58 +/- 0.13 | 2.24 +/- 0.10 | 0.0125 ★ |
| TNN d | Meditators | 10.78 +/- 0.24 | 10.28 +/- 0.19 | 0.0092 ★★ |
| TNN d | Controls | 10.25 +/- 0.23 | 10.17 +/- 0.24 | 0.3252 |

### 3.2 EC1 vs M1

| Metric | Group | EC1 (mean +/- SEM) | M1 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 3.10 +/- 0.16 | 3.83 +/- 0.22 | **0.0035 ★★★** |
| k | Controls | 3.10 +/- 0.15 | 3.47 +/- 0.18 | 0.0758 |
| Velocity | Meditators | 24.81 +/- 2.34 | 21.49 +/- 2.32 | 0.0677 |
| Velocity | Controls | 24.19 +/- 1.84 | 22.63 +/- 1.96 | 0.0981 |
| Log Volume | Meditators | 8.68 +/- 0.28 | 8.32 +/- 0.27 | 0.0432 ★ |
| Log Volume | Controls | 8.68 +/- 0.25 | 8.65 +/- 0.29 | 0.3469 |
| Radius | Meditators | 20.69 +/- 2.01 | 17.74 +/- 1.92 | 0.0585 |
| Radius | Controls | 20.00 +/- 1.54 | 18.69 +/- 1.64 | 0.1329 |
| PR | Meditators | 2.29 +/- 0.13 | 2.61 +/- 0.11 | 0.0086 ★★ |
| PR | Controls | 2.24 +/- 0.10 | 2.49 +/- 0.15 | 0.0817 |
| TNN d | Meditators | 10.28 +/- 0.19 | 12.11 +/- 0.31 | **<0.0001 ★★★** |
| TNN d | Controls | 10.17 +/- 0.24 | 10.86 +/- 0.34 | **0.0018 ★★★** |

### 3.3 EO1 vs M1

| Metric | Group | EO1 (mean +/- SEM) | M1 (mean +/- SEM) | p-value |
|---|---|---|---|---|
| k | Meditators | 3.52 +/- 0.16 | 3.83 +/- 0.23 | 0.1082 |
| k | Controls | 3.30 +/- 0.11 | 3.47 +/- 0.18 | 0.4850 |
| Velocity | Meditators | 20.44 +/- 2.28 | 21.49 +/- 2.36 | 0.0108 ★ |
| Velocity | Controls | 19.64 +/- 1.65 | 22.63 +/- 1.93 | **0.0002 ★★★** |
| Log Volume | Meditators | 7.98 +/- 0.29 | 8.32 +/- 0.28 | 0.0242 ★ |
| Log Volume | Controls | 8.14 +/- 0.25 | 8.65 +/- 0.29 | **0.0008 ★★★** |
| Radius | Meditators | 16.85 +/- 1.91 | 17.74 +/- 1.96 | 0.0077 ★★ |
| Radius | Controls | 16.08 +/- 1.35 | 18.69 +/- 1.61 | **0.0001 ★★★** |
| PR | Meditators | 2.68 +/- 0.15 | 2.61 +/- 0.11 | 0.4388 |
| PR | Controls | 2.58 +/- 0.13 | 2.49 +/- 0.15 | 0.3493 |
| TNN d | Meditators | 10.78 +/- 0.25 | 12.11 +/- 0.31 | **<0.0001 ★★★** |
| TNN d | Controls | 10.25 +/- 0.23 | 10.86 +/- 0.34 | 0.0104 ★ |
