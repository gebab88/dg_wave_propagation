#ifndef VANDERMONDE
#define VANDERMONDE

#include <armadillo>
using namespace arma;

mat VanDer(const uword &order, const vec &Stuetz);
mat TransposeDiffMatrix(const uword &order,const vec &Stuetz);

#endif // VANDERMONDE_HPP
