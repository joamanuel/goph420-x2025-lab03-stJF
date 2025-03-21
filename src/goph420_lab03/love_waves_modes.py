import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 8, 5
rcParams["figure.subplot.hspace"] = (0.5)
#from root_finding.modified_secant import root_secant_modified

def root_secant_modified(x0, dx, f, tol=1e-6, max_iter=100):
    errors = []
    x = x0
    for i in range(max_iter):
        fx = f(x)
        if np.isnan(fx):
            break
        dfdx_approx = (f(x + dx * x) - fx) / (dx * x)
        if dfdx_approx == 0 or np.isnan(dfdx_approx):
            break
        x_new = x - fx / dfdx_approx
        error = abs((x_new - x) / x_new)
        errors.append(error)
        if error < tol:
            return x_new, i+1, np.array(errors)
        x = x_new
    return x, max_iter, np.array(errors)


# Physical parameters
rho1 = 1800       # kg/m³
rho2 = 2500       # kg/m³
beta1 = 1900      # m/s
beta2 = 3200      # m/s
H = 4000          # m
modes = 5         # Number of dispersion modes to compute

def dispersion_function(phi, f):
    A = H**2 * (1 / beta1**2 - 1 / beta2**2)
    if phi == 0:
        return np.inf
    if A - phi**2 <= 0:
        return np.nan
    rhs = (rho2 / rho1) * np.sqrt(A - phi**2) / phi
    lhs = np.tan(2 * np.pi * f * phi)
    return lhs - rhs

def compute_cl(phi):
    inside = 1 / beta1**2 - (phi / H)**2
    if inside <= 0:
        return np.nan
    return 1 / np.sqrt(inside)

def main():
    frequencies = np.linspace(0.2, 2.5, 100)
    cl_modes = [[] for _ in range(modes)]
    lambda_modes = [[] for _ in range(modes)]
    freq_modes = [[] for _ in range(modes)]

    for f in frequencies:
        phi_roots = []
        for n in range(1, modes + 5):
            x0 = (n - 0.5) / (2 * f)
            g_phi = lambda phi: dispersion_function(phi, f)
            try:
                phi_root, _, _ = root_secant_modified(x0, dx=1e-6, f=g_phi)
                cL = compute_cl(phi_root)
                if (
                    not np.isnan(phi_root)
                    and not np.isnan(cL)
                    and beta1 < cL < beta2
                    and not any(np.isclose(phi_root, prev, atol=1e-3) for prev in phi_roots)
                ):
                    phi_roots.append(phi_root)
            except Exception:
                continue

        phi_roots = sorted(phi_roots)[:modes]
        for i, phi in enumerate(phi_roots):
            cL = compute_cl(phi)
            if not np.isnan(cL) and beta1 < cL < beta2:
                cl_modes[i].append(cL)
                lambda_modes[i].append(cL / f)
                freq_modes[i].append(f)

    # --- Plot Love Wave Velocity (c_L) vs Frequency ---
    plt.figure(figsize=(8, 5))
    for i in range(modes):
        plt.plot(freq_modes[i], cl_modes[i], label=f"Mode {i}")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Love Wave Velocity $c_L$ (m/s)")
    plt.title("Dispersion Curves: $c_L$ vs Frequency (Modified Secant)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Plots/cl_vs_f_secant.png")
    plt.show()

    # --- Plot Wavelength (lambda_L) vs Frequency ---
    plt.figure(figsize=(8, 5))
    for i in range(modes):
        plt.plot(freq_modes[i], lambda_modes[i], label=f"Mode {i}")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Wavelength $\\lambda_L$ (m)")
    plt.title("Wavelengths: $\\lambda_L$ vs Frequency (Modified Secant)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Plots/lambda_vs_f_secant.png")
    plt.show()

if __name__ == "__main__":
    main()
