#include "vandermonde.hpp"

mat VanDer(const uword &order,const vec &Stuetz){
    mat A=zeros<mat>(order,order);
    for (uword i=0; i<order;i++){
        A.col(i)=pow(Stuetz,i);
    }
    //A.print("A= ");
    return A;
}

mat TransposeDiffMatrix(const uword &order,const vec &Stuetz){
    mat Dtilde=zeros<mat>(order,order);
    for (uword l=1;l<order;l++){
        Dtilde.col(l)=(l)*pow(Stuetz,l-1);
    }
    mat Vandermonde=VanDer(order,Stuetz);
    //Vandermonde.print(" Vander=");
    //Stuetz.print("Stuetz =");
    mat DT=trans(Dtilde*inv(Vandermonde));
    //DT.print("DT=");
    return DT;
}

void GaussLegendre(const uword &order, vec &nodes, vec &weights){
    // Golub-Welsch: the Gauss-Legendre nodes are the eigenvalues of the
    // symmetric tridiagonal Jacobi matrix J (zero diagonal, off-diagonals
    // beta_k = k / sqrt(4 k^2 - 1)); the weights are 2*(first eigenvector
    // component)^2. Works for any order without hardcoding radicals.
    mat J=zeros<mat>(order,order);
    for (uword k=1; k<order; ++k){
        double beta = double(k)/sqrt(4.0*double(k)*double(k) - 1.0);
        J(k,k-1)=beta;
        J(k-1,k)=beta;
    }
    vec eigval;
    mat eigvec;
    eig_sym(eigval,eigvec,J);          // eigenvalues ascending -> nodes
    nodes=eigval;
    weights=2.0*square(trans(eigvec.row(0)));
}

vec LagrangeAtPoint(const vec &Stuetz, double xi){
    // l_j(xi) = prod_{m != j} (xi - z_m)/(z_j - z_m).
    uword order=Stuetz.n_elem;
    vec L=ones<vec>(order);
    for (uword j=0; j<order; ++j)
        for (uword m=0; m<order; ++m)
            if (m!=j)
                L(j)*=(xi-Stuetz(m))/(Stuetz(j)-Stuetz(m));
    return L;
}
