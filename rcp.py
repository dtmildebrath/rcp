import ctypes
import time

# Polynomial time version
lib = ctypes.CDLL("./librcp.so")

def _set_bindings():
    """ Set the I/O types for the C wrappers """
    lib.solve_rcp.restype = ctypes.c_int
    lib.solve_rcp.argtypes = (
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
    )

    # Brute force version
    lib.solve_rcp_brute_force.restype = ctypes.c_int
    lib.solve_rcp_brute_force.argtypes = (
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
    )

# Automatically set the I/O bindings for the C functions upon import
_set_bindings()

def solve_rcp(
    a, b, c, L, r, phi, tol=1e-10, display_timing=False, method="true", return_timing=False, call_set_bindings=False,
):
    """ Solve the RCP for the specified problem data.

    Args:
        a (ndarray):
            Coefficients of the mean of the random variable.
        b (ndarray):
            Coefficients of the variance of the random variable (assumed all
            strictly positive).
        c (ndarray):
            Objective function coefficients.
        L (ndarray):
            Upper bounds for each decision variable (assumed all strictly
            positive).
        r (float):
            Value of the right-hand side of the constraint.
        phi (float):
            Value corresponding to the risk tolerance in the chance constraint.
        tol (float, optional):
            Tolerance value used for multiple purposes in the algorithm.
            Default is 1e-10.
        display_timing (bool, optional):
            If True, prints timing information to stdout. This timing excludes
            the time needed to set C bindings and validate parts of the input.
        method (str, optional):
            Either `true` or `brute`. If `true`, solves the problem using the
            algorithm in the paper. If `brute` solves the problem using brute
            force. Default is `true`.
        return_timing (bool, optional):
            If True, returns wall clock execution time in seconds in addition
            to the optimal solution. Default is False.
        call_set_bindings (bool, optional):
            If True, sets the input/return types for all the C functions. This
            should be called automatically any time the rcp module is imported.
            This flag is included for "legacy" reasons. Default is False.

    Returns:
        list of float:
            Optimal solution vector.
        float (if return_timing=True):
            Total wall clock execution time in seconds. 
            
    """
    method = method.lower()
    if method not in ["true", "brute"]:
        raise ValueError("Method must be `true' or `brute'")
    n = len(a)
    soln = [0.0 for _ in range(n)]

    if call_set_bindings:
        _set_bindings()

    # Cast everything to C types
    c_args = (
        ctypes.c_int(n),
        (ctypes.c_double * n)(*a),
        (ctypes.c_double * n)(*b),
        (ctypes.c_double * n)(*c),
        (ctypes.c_double * n)(*L),
        ctypes.c_double(r),
        ctypes.c_double(phi),
        ctypes.c_double(tol),
        (ctypes.c_double * n)(*soln),
    )

    if display_timing or return_timing:
        t0 = time.time()

    if method == "true":
        status = lib.solve_rcp(*c_args)
    else:
        status = lib.solve_rcp_brute_force(*c_args)
        
    if status != 0:
        raise RuntimeError("Solve failed")

    if display_timing or return_timing:
        t1 = time.time()
        print(f"Time in C call: {t1-t0:.2f}s")

    soln = [c_args[-1][i] for i in range(n)]

    if return_timing:
        return soln, t1 - t0

    return soln
