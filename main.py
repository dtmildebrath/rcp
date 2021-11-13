import numpy as np
import rcp


def generate_random_instance(n, rng, round_data=False, remove_trivial_coords=False):
    """ Generate random problem instances.

    Args:
        n (int):
            Number of variables in the problem.
        rng (Generator):
            Numpy random number generator.
        round_data (bool, optional):
            If True, rounds the entries in a, b and c to two decimal places,
            and rounds r to one decimal place. Rounding is used to simplify
            debugging. Default is False.
        remove_trivial_coords (bool, optional):
            If True, ensures that no coordinates satisfy c_i <= 0 and a_i >= 0.
            Such coordinates are always 0 at optimality. Default is False.
    """
    a = rng.uniform(low=-1., high=1., size=n)
    b = rng.uniform(low=1e-2, high=1., size=n)
    c = rng.uniform(low=-1., high=1., size=n)
    phi = 1.0
    r = rng.uniform(low=0., high=0.5)
    L = np.ones(n)

    if round_data:
        a = np.around(a, decimals=2)
        b = np.around(b, decimals=2)
        c = np.around(c, decimals=2)
        r = np.around(r, decimals=1)

    if remove_trivial_coords:
        ix = np.where((c <= 0.) & (a >= 0.))[0]
        a[ix] *= -1.

    return a, b, c, L, r, phi


def brute_force_checker(n):
    """ Solve instances until we find one where the brute force solution and
    the algorithmic solution do not agree.

    Used for debugging and correction validation.

    Input:
        n (int):
            Size of instances to check. Cannot be too large or the brute force
            checker will fail.
    """
    rng = np.random.default_rng()

    while True:
        a, b, c, L, r, phi = generate_random_instance(n, rng)

        print("a =", list(a))
        print("b =", list(b))
        print("c =", list(c))
        print("L =", list(L))
        print("r =", r)
        print("phi =", phi)

        soln = rcp.solve_rcp(a, b, c, L, r, phi, display_timing=True, method="true")
        soln = np.array(soln)
        objval = np.dot(soln, c)

        soln_brute = rcp.solve_rcp(a, b, c, L, r, phi, display_timing=True, method="brute")
        soln_brute = np.array(soln_brute)
        objval_brute = np.dot(soln_brute, c)
        
        if np.abs(objval - objval_brute) > 1e-5:
            print("FOUND DIFFERENT SOLUTIONS")
            print("Brute force:")
            print("\tSolution:", soln_brute)
            print(f"\tObjective value: {objval_brute:.4f}")
            print("Algorithm:")
            print("\tSolution:", soln)
            print(f"\tObjective value: {objval:.4f}")
            break
        else:
            print("Brute:", soln_brute)
            print("Algor:", soln)


def generate_and_solve_one_instance(n, seed=None, solve_brute_too=False):
    if seed is None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    a, b, c, L, r, phi = generate_random_instance(n, rng, round_data=True)

    soln = rcp.solve_rcp(a, b, c, L, r, phi, display_timing=True, method="true")
    soln = np.array(soln)
    objval = np.dot(soln, c)

    if solve_brute_too:
        soln_brute = rcp.solve_rcp(a, b, c, L, r, phi, display_timing=True, method="brute")
        soln_brute = np.array(soln_brute)
        objval_brute = np.dot(soln_brute, c)

    print("Algo soln:", soln)
    if solve_brute_too:
        print("Brute: ", soln_brute)
    print("Algo: ", objval)
    if solve_brute_too:
        print("Brute: ", objval_brute)


def paper_numerical_test():
    NUM_TRIALS = 10
    n_values = (250, 500, 1000, 2000)
    seed = 1774  # For reproducibility

    rng = np.random.default_rng(seed)

    for n in n_values:
        avg_time = 0.0
        for t in range(NUM_TRIALS):
            instance = generate_random_instance(n, rng)
            _, time = rcp.solve_rcp(*instance, return_timing=True, display_timing=False)
            avg_time += time
            print(f"n = {n} instance {t+1} time {time:.6f}s")
        avg_time /= NUM_TRIALS
        print(f"\nn = {n} average time {time:.6f}s")
        print()


if __name__ == "__main__":
    paper_numerical_test()
