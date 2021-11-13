#ifndef RCP_H
#define RCP_H

static const int ZERO_COORD = 0;
static const int ONE_COORD = 1;
static const int PLUS_COORD = 2;
static const int MINUS_COORD = 3;

/* This structure holds a scaled problem instance
 * (i.e. the upper bounds are all assumed to equal 1)
 */
struct ScaledInstance {
    int     n;
    double* a; 
    double* b; 
    double* c; 
    double  phi;
    double  r;
    double  b_max; // Sum of all entries of b
};

struct ListTuple {
    int    index;
    double value;
    int    type;
};

void solve_quadratic(double a, double b, double c, double* smaller_root, double* larger_root);
void get_roots(double z1, double z2, double z3, double z4, double phi, double r, double* smaller_root, double* larger_root);
double solve_single_edge_subproblem(double z1, double z2, double z3, double z4, double obj_coeff, double phi, double r, double tol);
int double_compare(const void *a, const void *b);
int pair_compare(const void *a, const void *b);
int tuple_compare(const void *a, const void *b);
void build_Y1(struct ScaledInstance m, double tol, int *length, double *Y1);
void build_Y2(struct ScaledInstance m, double tol, int *length, double *Y2);
int is_plus_coord(double a, double c);
int is_minus_coord(double a, double c);
int is_one_coord(double a, double c);
int is_zero_coord(double a, double c);
double get_knapsack_rhs(struct ScaledInstance m, double y);
void get_knapsack_coeffs(struct ScaledInstance m, double y, double *p);
double solve_general_knapsack_problem(int n, double *a, double *c, double b, int *is_feasible, double *solution);
double solve_knapsack_problem(struct ScaledInstance m, double y, double tol, int *is_feasible, double *solution);
double solve_edge_sequence(struct ScaledInstance m, double y, double tol, double* solution);
void solve_scaled_rcp(struct ScaledInstance m, double tol, double* solution);
int solve_rcp(int n, double *a, double *b, double *c, double *L, double r, double phi, double tol, double *soln);
void solve_scaled_rcp_brute_force(struct ScaledInstance m, double tol, double* solution);
int solve_rcp_brute_force(int n, double *a, double *b, double *c, double *L, double r, double phi, double tol, double *soln);

#endif
