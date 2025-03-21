import numpy as np

def root_newton_raphson(x0, f, dfdx, tol=1e-6, max_iter=100):
    errors = []
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfx = dfdx(x)
        if dfx == 0:
            raise ZeroDivisionError("Derivative is zero.")
        x_new = x - fx / dfx
        error = abs((x_new - x) / x_new)
        errors.append(error)
        if error < tol:
            break
        x = x_new
    return x, i+1, np.array(errors)



def root_secant_modified(x0, dx, f, tol=1e-6, max_iter=100):
    errors = []
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfdx_approx = (f(x + dx * x) - fx) / (dx * x)
        if dfdx_approx == 0:
            raise ZeroDivisionError("Approx derivative is zero.")
        x_new = x - fx / dfdx_approx
        error = abs((x_new - x) / x_new)
        errors.append(error)
        if error < tol:
            break
        x = x_new
    return x, i+1, np.array(errors)


