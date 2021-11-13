/*
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rcp.h"


/* Compute the roots of the quadratic ax^2 + bx + c
 * The values smaller_root and larger_root are set to NaN if there are no real
 * roots. Otherwise, computes the real roots in ascending order (returns the
 * same value twice if there is a single root with multiplicity 1).
 */
void solve_quadratic(double a, double b, double c, double* smaller_root, double* larger_root)
{
    static double tol = 1e-15;
    if(fabs(a) < tol) { // Function is linear
        *smaller_root = -1. * c / b;
        *larger_root = -1. * c / b;
        return;
    }
    double quad = b * b - 4 * a * c;
    if(quad < 0){ // No real roots
        *smaller_root = nan("1"); // Hope this works!
        *larger_root = nan("1");
        return;
    }
    double discrim = sqrt(quad);
    if(fabs(discrim) < tol) { // Single repeated root
        *smaller_root = -b / (2 * a);
        *larger_root = -b / (2 * a);
        return;
    }
    double x1, x2, tmp;
    if(b >= 0) {
        x1 = -(b + discrim) / (2 * a);
        x2 = -(2 * c) / (b + discrim);
    } else {
        x1 = (2 * c) / (-b + discrim);
        x2 = (-b + discrim) / (2 * a);
    }
    if(x1 > x2) {
        tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
    *smaller_root = x1;
    *larger_root = x2;
}


/* Compute all values of x for which f(x)=0 where
 * f(x) = z1 * x + z2 + phi * sqrt(z3 * x + z4) - r
 */
void get_roots(
    double z1, double z2, double z3, double z4, double phi, double r, double* smaller_root, double* larger_root
) {
    double aa, bb, cc;
    aa = z1 * z1;
    bb = -1. * (2 * (r - z2) * z1 + phi * phi * z3);
    cc = (r - z2) * (r - z2) - phi * phi * z4;
    solve_quadratic(aa, bb, cc, smaller_root, larger_root);
}


/* Solve the 1D optimization problem
 * 
 * max obj_coeff * x
 *  st f(x) <= 0
 *     0 <= x <= 1
 *
 * where f(x) = z1 * x + z2 + phi * sqrt(z3 * x + z4) - r
 *
 * Returns -1 if the problem is infeasible. Otherwise returns
 * the optimal value of x (_not_ the optimal objective value)
 */
double solve_single_edge_subproblem(
    double z1, double z2, double z3, double z4, double obj_coeff, double phi, double r, double tol
) {
    double fmn = z2 - (z4 * z1 / z3) - r;
    double x_opt; // -1 if the problem is infeasible
    double r0, r1, true_root;

    // f(x) > 0 for all u in the domain of f
    if((z1 >= -tol) && (fmn > tol)) {
        x_opt = -1.;
    // f has a single root r0, and f(x) <= 0 for all x in [mn, r0]
    } else if((z1 >= -tol) && (fmn <= tol)) {
        get_roots(z1, z2, z3, z4, phi, r, &r0, &r1);
        true_root = r0;

        if(true_root >= -tol) {
            if(obj_coeff < 0) { // Minimization problem
                x_opt = 0; 
            } else { // Maximization problem
                if(true_root < 0.) {
                    x_opt = 0.;
                } else if(true_root > 1.) {
                    x_opt = 1.;
                } else {
                    x_opt = true_root;
                }
            }
        } else {
            x_opt = -1.;
        }
    // f has a single root r0, and f(x) <= 0 for all x in [r0, inf)
    } else if((z1 < -tol) && (fmn > tol)) {
        get_roots(z1, z2, z3, z4, phi, r, &r0, &r1);
        true_root = r1;

        if(true_root <= 1. + tol) {
            if(obj_coeff < 0) {
                if(true_root < -tol) {
                    x_opt = 0;
                } else {
                    x_opt = true_root;
                }
            } else {
                x_opt = 1.;
            }
        } else {
            x_opt = -1.;
        }
    // f may have either zero or two roots
    } else if((z1 < -tol) && (fmn <= tol)) {
        double thing, arg_max, mx;

        thing = (phi * z3) / (2 * z1);
        arg_max = (thing * thing - z4) / z3;
        mx = z2 + z1 * arg_max + phi * sqrt(z4 + z3 * arg_max) - r;

        if(mx > tol) {
            // f has two roots r0, r1, and f(x) <= 0
            // for all x in the union of [mn, r0] and [r1, inf).
            get_roots(z1, z2, z3, z4, phi, r, &r0, &r1);

            if(obj_coeff < 0) {
                if(r0 > -tol) {
                    x_opt = 0.;
                } else {
                    if(r1 < tol) {
                        x_opt = 0.;
                    } else if(r1 > 1. + tol) {
                        x_opt = -1.;
                    } else {
                        x_opt = r1;
                    }
                }
            } else {
                if(r1 <= 1. + tol) {
                    x_opt = 1.;
                } else if(r0 < -tol) { // and r1 > 1 (because the first case)
                    x_opt = -1;
                } else { // r0 >= 0 and r1 > 1
                    if(1. < r0) {
                        x_opt = 1.;
                    } else {
                        x_opt = r0;
                    }
                }
            }
        } else {
            // f has no roots, and f(x) <= 0 on its entire domain.
            if(obj_coeff < 0) {
                x_opt = 0.;
            } else {
                x_opt = 1.;
            }
        }
    }

    return x_opt;
}


/* Comparator function for sorting floats */
int double_compare(const void *a, const void *b)
{
    /* Needed to use the qsort() function */
    /* https://stackoverflow.com/questions/1787996/c-library-function-to-perform-sort */
    double *x = (double *) a;
    double *y = (double *) b;
    if(*x < *y)
        return -1;
    else if(*x > *y)
        return 1;
    return 0;
}


/* Comparator function for sorting tuple structs */
int tuple_compare(const void *a, const void *b)
{
    /* Needed for the argsort */
    /* https://stackoverflow.com/questions/36714030/c-sort-float-array-while-keeping-track-of-indices */
    struct ListTuple *a1 = (struct ListTuple *) a;
    struct ListTuple *a2 = (struct ListTuple *) b;
    if((*a1).value < (*a2).value)
        return -1;
    else if ((*a1).value > (*a2).value)
        return 1;
    else
        return 0;
}


void build_Y1(struct ScaledInstance m, double tol, int *length, double *Y1)
{
    double v; 
    *length = 0;
    for (int i=0; i<m.n; i++) {
        if(fabs(m.a[i]) > tol) {
            v = 0.5 * m.phi * m.b[i] / m.a[i];
            v = v * v;
            if(v <= m.b_max) {
                Y1[*length] = v;
                *length = *length + 1;
            }
        }
    }
}


void build_Y2(struct ScaledInstance m, double tol, int *length, double *Y2)
{
    double v, v2, denom;
    *length = 0;
    for (int i=0; i<m.n; i++) {
        for (int j=i+1; j<m.n; j++) {
            denom = m.a[j] * m.c[i] - m.a[i] * m.c[j];
            if(fabs(denom) > tol) {
                v = 0.5 * m.phi * (m.b[i] * m.c[j] - m.b[j] * m.c[i]) / denom;
                v2 = v * v;
                if((v >= 0) && (v2 <= m.b_max)) {
                    Y2[*length] = v2;
                    *length = *length + 1;
                }
            }
        }
    }
}


int is_plus_coord(double a, double c)
{
    return ((c > 0) && (a > 0));
}


int is_minus_coord(double a, double c)
{
    return ((c < 0) && (a < 0));
}


int is_one_coord(double a, double c)
{
    return ((c >= 0) && (a <= 0));
}


int is_zero_coord(double a, double c)
{
    return ((c <= 0) && (a >= 0));
}


double get_knapsack_rhs(struct ScaledInstance m, double y)
{
    return m.r * sqrt(y) - 0.5 * m.phi * y;
}


void get_knapsack_coeffs(struct ScaledInstance m, double y, double *p)
{
    for (int i=0; i<m.n; i++) {
        p[i] = sqrt(y) * m.a[i] + 0.5 * m.phi * m.b[i];
    }
}


/* Solves a fractional knapsack problem max{cx | ax<=b, 0<=x<=1} with no
 * assumptions on the signs of a, b and c.
 * Sets the is_feasible flag to indicate infeasibility
 * Returns the objective value
 */
double solve_general_knapsack_problem(int n, double *a, double *c, double b, int *is_feasible, double *solution)
{
    int num_non_trivial_coords;
    double rhs, rhs_shift, objval;

    struct ListTuple *non_trivial_coords = (struct ListTuple *) malloc(n * sizeof(struct ListTuple));

    // Divide up the types of coordinates
    num_non_trivial_coords = 0;
    rhs_shift = 0.0;
    objval = 0.0;
    for (int i=0; i<n; i++) {
        if(is_one_coord(a[i], c[i])) {
            solution[i] = 1.0;
            objval += c[i];
            rhs_shift += fabs(a[i]);
        } else if(is_zero_coord(a[i], c[i])) {
            solution[i] = 0.0;
        } else if(is_plus_coord(a[i], c[i])) {
            non_trivial_coords[num_non_trivial_coords].index = i;
            non_trivial_coords[num_non_trivial_coords].value = a[i] / c[i];
            non_trivial_coords[num_non_trivial_coords].type = PLUS_COORD;
            num_non_trivial_coords++;
        } else { // is_minus_coord
            non_trivial_coords[num_non_trivial_coords].index = i;
            non_trivial_coords[num_non_trivial_coords].value = a[i] / c[i];
            non_trivial_coords[num_non_trivial_coords].type = MINUS_COORD;
            num_non_trivial_coords++;
            rhs_shift += fabs(a[i]);
        }
    }
    rhs = b + rhs_shift;

    if (rhs < 0) { // Problem is infeasible
        *is_feasible = 0;
        objval = -INFINITY;  // Just in case
        free(non_trivial_coords);
        return objval;
    }
    *is_feasible = 1;

    if (num_non_trivial_coords == 0) {
        // objval has already been set
        free(non_trivial_coords);
        return objval;
    }

    double lhs, ak, frac_amount;
    int frac_index;

    // Solve the knapsack problem with a sort
    qsort(non_trivial_coords, num_non_trivial_coords, sizeof(*non_trivial_coords), tuple_compare);

    lhs = 0.0;
    frac_index = 0;
    for (int i=0; i<num_non_trivial_coords; i++) {
        ak = fabs(a[non_trivial_coords[i].index]);
        if(lhs + ak <= rhs) {
            lhs += ak;
        } else {
            frac_index = i;
            frac_amount = (rhs - lhs) / ak;
            break;
        }
    }
    
    // Reconstruct the solution
    for (int i=0; i<frac_index; i++) {
        if(non_trivial_coords[i].type == PLUS_COORD) {
            objval += c[non_trivial_coords[i].index];
            solution[non_trivial_coords[i].index] = 1.0; 
        } else {
            solution[non_trivial_coords[i].index] = 0.0;
        }
    }
    if(non_trivial_coords[frac_index].type == PLUS_COORD) {
        objval += c[non_trivial_coords[frac_index].index] * frac_amount;
        solution[non_trivial_coords[frac_index].index] = frac_amount;
    } else {
        objval += c[non_trivial_coords[frac_index].index] * (1.0 - frac_amount);
        solution[non_trivial_coords[frac_index].index] = 1.0 - frac_amount;
    }
    for (int i=frac_index+1; i<num_non_trivial_coords; i++) {
        if(non_trivial_coords[i].type == PLUS_COORD) {
            solution[non_trivial_coords[i].index] = 0.0; 
        } else {
            objval += c[non_trivial_coords[i].index];
            solution[non_trivial_coords[i].index] = 1.0;
        }
    }

    free(non_trivial_coords);

    return objval;
}


/* Solves the particular parametric fractional knapsack problem presented in
 * the paper (parameterized by y>=0).
 * Sets the is_feasible flag to indicate feasibility
 * Returns objval
 */
double solve_knapsack_problem(struct ScaledInstance m, double y, double tol, int *is_feasible, double *solution)
{
    double *p;
    double rhs, objval;
    p = (double *) malloc(m.n * sizeof(double));

    get_knapsack_coeffs(m, y, p);
    rhs = get_knapsack_rhs(m, y);

    objval = solve_general_knapsack_problem(m.n, p, m.c, rhs, is_feasible, solution);

    free(p);

    return objval;
}


/* Algorithm 1 in the paper.
 * Returns the best objective value (or -1 if infeasible)
 * This function is the bulk of the algorithm so speed is good.
 * 
 * TODO: If n_ell is 0, then the edge sequence reduces to a single point, which
 * may or may not be feasible. This needs to be handled in both the code and
 * the paper.
 */
double solve_edge_sequence(struct ScaledInstance m, double y, double tol, double* solution)
{
    double a, b, c, alpha_max, z_max, alpha, z;
    int k_max, num_non_trivial_indices;
    double *p;

    p = (double *) malloc(m.n * sizeof(double));
    get_knapsack_coeffs(m, y, p);

    // Build the index sets
    struct ListTuple *non_trivial_indices = (struct ListTuple *) malloc(m.n * sizeof(struct ListTuple));
    num_non_trivial_indices = 0;
    for (int i=0; i<m.n; i++) {
        if (is_plus_coord(p[i], m.c[i])) {
            non_trivial_indices[num_non_trivial_indices].index = i;
            non_trivial_indices[num_non_trivial_indices].value = p[i] / m.c[i];
            non_trivial_indices[num_non_trivial_indices].type = PLUS_COORD;
            num_non_trivial_indices++;
        } else if (is_minus_coord(p[i], m.c[i])) {
            non_trivial_indices[num_non_trivial_indices].index = i;
            non_trivial_indices[num_non_trivial_indices].value = p[i] / m.c[i];
            non_trivial_indices[num_non_trivial_indices].type = MINUS_COORD;
            num_non_trivial_indices++;
        }
    }

    z_max = -1.0;  // Optimal objval will always be >= 0

    a = 0.0; b = 0.0; c = 0.0;
    for (int i=0; i<m.n; i++) {
        if (is_one_coord(p[i], m.c[i]) || is_minus_coord(p[i], m.c[i])) {
            a += m.a[i]; b += m.b[i]; c += m.c[i];
        }
    }

    // Build the permutation pi
    qsort(non_trivial_indices, num_non_trivial_indices, sizeof(*non_trivial_indices), tuple_compare);

    // Solve each edge subproblem
    for (int k=0; k<num_non_trivial_indices; k++) {
        if (non_trivial_indices[k].type == MINUS_COORD) {
            a -= m.a[non_trivial_indices[k].index];
            b -= m.b[non_trivial_indices[k].index];
            c -= m.c[non_trivial_indices[k].index];
        }
        alpha = solve_single_edge_subproblem(
            m.a[non_trivial_indices[k].index],
            a,
            m.b[non_trivial_indices[k].index],
            b,
            m.c[non_trivial_indices[k].index],
            m.phi,
            m.r,
            tol
        );
        if (alpha > -0.5) { // Subproblem was feasible
            z = c + alpha * m.c[non_trivial_indices[k].index];
            if (z > z_max) {
                z_max = z;
                alpha_max = alpha;
                k_max = k;
            }
        }
        if (non_trivial_indices[k].type == PLUS_COORD) {
            a += m.a[non_trivial_indices[k].index];
            b += m.b[non_trivial_indices[k].index];
            c += m.c[non_trivial_indices[k].index];
        }
    }

    // Reconstruct the best solution
    if (z_max > -0.5) {
        for (int i=0; i<m.n; i++) {
            if (is_one_coord(p[i], m.c[i]) || is_minus_coord(p[i], m.c[i])) {
                solution[i] = 1.0;
            } else {
                solution[i] = 0.0;
            }
        }
        for (int k=0; k<k_max; k++) {
            if (non_trivial_indices[k].type == PLUS_COORD) {
                solution[non_trivial_indices[k].index] = 1.0;
            } else { // Minus coord
                solution[non_trivial_indices[k].index] = 0.0;
            }
        }
        solution[non_trivial_indices[k_max].index] = alpha_max;
    }

    free(non_trivial_indices);
    free(p);

    return z_max;
}


/* Solve RCP with coordinates scaled so the upper bounds L_i=1 for all i.
 * Does not do any input checking.
 */
void solve_scaled_rcp(struct ScaledInstance m, double tol, double* solution)
{
    double *Y;
    double *candidate_solution;
    double z_max, z_candidate, midpoint;
    int length_Y1, length_Y2, max_length_Y, L, is_feasible;

    z_max = -INFINITY;

    max_length_Y = m.n + ((m.n * (m.n - 1)) / 2) + 2;

    Y = (double *) malloc(max_length_Y * sizeof(double));
    candidate_solution = (double *) malloc(m.n * sizeof(double));
    build_Y1(m, tol, &length_Y1, Y);

    for(int i=0; i<length_Y1; i++) {
        z_candidate = solve_knapsack_problem(m, Y[i], tol, &is_feasible, candidate_solution);
        if ((is_feasible) && (z_candidate > z_max)) { // Checks feasibility
            z_max = z_candidate;
            for (int j=0; j<m.n; j++) solution[j] = candidate_solution[j];
        }
    }

    build_Y2(m, tol, &length_Y2, &Y[length_Y1]);

    /*
    printf("Y1\n");
    for (int i=0; i<length_Y1; i++) {
        printf("\t%.4f\n", Y[i]);
    }
    printf("Y2\n");
    for (int i=length_Y1; i<length_Y1+length_Y2; i++) {
        printf("\t%.4f\n", Y[i]);
    }
    */
    L = length_Y1 + length_Y2;
    Y[L] = 0.0;
    Y[L + 1] = m.b_max;
    qsort(Y, L + 2, sizeof(*Y), double_compare);

    /*
    for (int i=0; i<L; i++) {
        printf("\t%.4f\n", Y[i]);
    }
    */

    for (int ell=0; ell<=L; ell++) {
        midpoint = 0.5 * (Y[ell] + Y[ell+1]); 
        z_candidate = solve_edge_sequence(m, midpoint, tol, candidate_solution);
        // Don't need to check infeasibility here
        // because z_candidate would be -1 if it were infeasible
        if (z_candidate > z_max) {
            z_max = z_candidate;
            for (int j=0; j<m.n; j++) solution[j] = candidate_solution[j];
        }
    }

    free(candidate_solution);
    free(Y);
}


/* Solve the general RCP.
 * Also performs an O(n) check on the input parameters.
 */
int solve_rcp(
    int n, double *a, double *b, double *c, double *L, double r, double phi, double tol, double *soln
) {
    // Check the input before allocating anything (input check is O(n))
    if(r < 0) {
        fprintf(stderr, "r must be nonnegative\n");
        return -1;
    }
    if(phi < 0) {
        fprintf(stderr, "phi must be nonnegative\n");
        return -1;
    }
    for(int i=0; i<n; i++) {
        if(b[i] <= 0) {
            fprintf(stderr, "All entries of b must be strictly positive\n");
            return -1;
        }
        if(L[i] <= 0) {
            fprintf(stderr, "All entries of L must be strictly positive\n");
            return -1;
        }
        /*
        if((c[i] <= 0) && (a[i] >= 0)) {
            fprintf(stderr, "Coordinate %d with c_i <= 0 and a_i >= 0 can be removed\n", i);
            return -1;
        }
        */
    }

    // Input is clear
    struct ScaledInstance m;

    m.n = n;
    m.a = (double *) malloc(n * sizeof(double));
    m.b = (double *) malloc(n * sizeof(double));
    m.c = (double *) malloc(n * sizeof(double));
    m.r = r;
    m.phi = phi;

    // Scale the input
    for(int i=0; i<m.n; i++) {
        m.a[i] = a[i] * L[i];
        m.b[i] = b[i] * L[i];
        m.c[i] = c[i] * L[i];
    }

    // Compute upper bound on b
    m.b_max = 0.0;
    for (int i=0; i<n; i++) {
        m.b_max += m.b[i];
    }

    solve_scaled_rcp(m, tol, soln);

    free(m.a);
    free(m.b);
    free(m.c);

    // Scale the solution back
    for(int i=0; i<m.n; i++)
        soln[i] *= L[i];

    return 0;
}


/* For debugging. Main usage is by invoking the python wrapper. */
int main(int argc, char **argv)
{
    int n = 3;
    double a[] = {-0.48, -0.4, 0.63};
    double b[] = {0.1, 0.6, 0.73};
    double c[] = {-0.62, -0.89, -0.45};
    double L[] = {1.0, 1.0, 1.0};
    double phi = 1.0;
    double r = 0.3;
    double *soln = (double *) malloc(n * sizeof(double));
    double objval;
    int status;

    status = solve_rcp(n, a, b, c, L, r, phi, 1e-10, soln);

    if(status != 0) {
        fprintf(stderr, "Solve failed\n");
        free(soln);
        return -1;
    }

    printf("Algorithmic solution\n");
    objval = 0.0;
    for(int i=0; i<n; i++) {
        printf("  x[%d] = %.4f\n", i, soln[i]);
        objval += soln[i] * c[i];
    }
    printf("z* = %.4f\n", objval);

    status = solve_rcp_brute_force(n, a, b, c, L, r, phi, 1e-10, soln);
    if(status != 0) {
        fprintf(stderr, "Solve failed\n");
        free(soln);
        return -1;
    }

    printf("Brute force solution\n");
    objval = 0.0;
    for(int i=0; i<n; i++) {
        printf("  x[%d] = %.4f\n", i, soln[i]);
        objval += soln[i] * c[i];
    }
    printf("z* = %.4f\n", objval);

    free(soln);

    return 0;
}
