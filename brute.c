#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rcp.h"

void solve_scaled_rcp_brute_force(
    struct ScaledInstance m, double tol, double* solution
) {
    int *vertex = (int *) malloc(m.n * sizeof(int));
    int *best_vertex = (int *) malloc(m.n * sizeof(int));
    int best_index, ub;
    double z2, z4, frac, objval, best_objval, best_frac;

    ub = 1;
    for(int i=0; i<m.n; i++) ub *= 2; 

    best_objval = -INFINITY;

    for(int i=0; i<ub-1; i++) { // ub - 1 to skip all-ones vertex
        for(int j=0; j<m.n; j++) {
            vertex[j] = (i >> j) & 1;  // j^th bit of i
        }
        for(int j=0; j<m.n; j++) {
            if(vertex[j] == 1)
                continue;
            z2 = 0;
            z4 = 0;
            for(int k=0; k<m.n; k++) {
                z2 += m.a[k] * vertex[k]; 
                z4 += m.b[k] * vertex[k]; 
            }
            frac = solve_single_edge_subproblem(
                m.a[j], z2, m.b[j], z4, m.c[j], m.phi, m.r, tol
            );
            if(frac < -0.5) // Infeasible subproblem
                continue;
            objval = 0.0;
            for(int k=0; k<m.n; k++) {
                if(k == j) {
                    objval += frac * m.c[j];
                } else if(vertex[k] == 1) {
                    objval += m.c[k];
                }
            }
            if(objval > best_objval) {
                best_objval = objval;
                best_frac = frac;
                best_index = j;
                for(int q=0; q<m.n; q++) {
                    best_vertex[q] = vertex[q];
                }
            }
        }
    }


    // Write the solution to the output buffer
    for(int i=0; i<m.n; i++) {
        if(best_vertex[i] == 1) {
            solution[i] = 1.0;
        } else {
            solution[i] = 0.0;
        }
    }
    solution[best_index] = best_frac;

    free(vertex);
    free(best_vertex);
}

int solve_rcp_brute_force(
    int n, double *a, double *b, double *c, double *L, double r, double phi, double tol, double *soln
) {
    // Check the input before allocating anything
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

    solve_scaled_rcp_brute_force(m, tol, soln);

    free(m.a);
    free(m.b);
    free(m.c);

    // Scale the solution back
    for(int i=0; i<m.n; i++)
        soln[i] *= L[i];

    return 0;
}

