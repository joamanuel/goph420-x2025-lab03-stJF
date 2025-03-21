# goph420-x2025-lab03-stJF
GOPH 420/699 – Inv. and Param. Est. for Geoph. (W2025) 

Lab Assignment #3

Root Finding Using the Newton-Raphson or Modified Secant Method.

A Python script is developed to nd multiple roots (corresponding to dierent modes) of the
dispersion equation over a range of frequencies. A key aspect of the algorithm is an eective initial
guess strategy based on the asymptotic behavior of the tangent function within the dispersion
equation, allowing for ecient and accurate identication of multiple modes

## 🧪 Testing Root-Finding Functions

You can run a simple test script to verify the implementation of both the Newton-Raphson and Modified Secant methods.

The test function used is:
\[
f(x) = x^2 - 2
\]
which has a known analytical root at \( \sqrt{2} \approx 1.4142135623730951 \).

### 📂 Script: `test_root_finding.py`

This script is located in the project root and tests both methods:

- `root_newton_raphson()` from `root_finding_functions.py`
- `root_secant_modified()` from `root_finding_functions.py`

### ▶️ How to Run the Tests

Make sure you're in goph420-x2025-lab03-stJF/root_finding_functions/, then run:

```bash
python test_root_finding.py


## 🌊 Love Wave Dispersion Modeling

The script `love_waves_modes.py` computes the dispersion curves for Love waves in a two-layer system using the **Modified Secant Method** for root finding.

It solves the dispersion equation numerically for a range of frequencies and identifies multiple modes (roots) where valid Love wave velocities exist. The results are visualized in two plots (They are stored in the Plots folder located at goph420-x2025-lab03-stJF/src/goph420_lab03/Plots/):

- **Love wave velocity \( c_L \) vs. frequency**
- **Wavelength \( \lambda_L = c_L / f \) vs. frequency**

---

### ▶️ How to Run the Script

From the goph420-x2025-lab03-stJF/src/goph420_lab03, run:

```bash
python love_waves_modes.py
