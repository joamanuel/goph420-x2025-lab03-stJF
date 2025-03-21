import numpy as np
from root_finding_functions import root_newton_raphson
from root_finding_functions import root_secant_modified

# --- Test function: f(x) = x^2 - 2, root at sqrt(2) ---
def f(x):
    return x**2 - 2

def dfdx(x):
    return 2 * x

def print_function_description():
    print("=" * 50)
    print("Root-Finding Test Case")
    print("Function being tested: f(x) = x^2 - 2")
    print("Analytical root: sqrt(2) ≈ 1.4142135623730951")
    print("=" * 50)

def test_newton_raphson():
    x0 = 1.0
    root, iterations, errors = root_newton_raphson(x0, f, dfdx)
    print("\n[Newton-Raphson Method]")
    print(f"  Initial guess: {x0}")
    print(f"  Root: {root}")
    print(f"  Iterations: {iterations}")
    print(f"  Final relative error: {errors[-1]:.2e}")
    print(f"  Expected: {np.sqrt(2)}")

def test_modified_secant():
    x0 = 1.0
    dx = 1e-6
    root, iterations, errors = root_secant_modified(x0, dx, f)
    print("\n[Modified Secant Method]")
    print(f"  Initial guess: {x0}, dx: {dx}")
    print(f"  Root: {root}")
    print(f"  Iterations: {iterations}")
    print(f"  Final relative error: {errors[-1]:.2e}")
    print(f"  Expected: {np.sqrt(2)}")

def main():
    print_function_description()
    test_newton_raphson()
    test_modified_secant()

if __name__ == "__main__":
    main()
