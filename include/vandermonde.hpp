#ifndef VANDERMONDE
#define VANDERMONDE

#include <armadillo>
using namespace arma;

mat VanDer(const uword &order, const vec &Stuetz);
mat TransposeDiffMatrix(const uword &order,const vec &Stuetz);

// Gauss-Legendre nodes & weights on [-1,1] for an arbitrary order
// (Golub-Welsch). Lets the solver support any order, not just the
// hardcoded 1/2/5 cases. nodes are returned in ascending order; weights sum to 2.
void GaussLegendre(const uword &order, vec &nodes, vec &weights);

// Values of all `order` Lagrange basis polynomials (built on `Stuetz`)
// evaluated at a single point xi. Used to build the edge-trace operator psi.
vec LagrangeAtPoint(const vec &Stuetz, double xi);

#endif // VANDERMONDE_HPP
