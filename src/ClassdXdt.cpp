#include "ClassdXdt.hpp"
#include "vandermonde.hpp"
#include <omp.h>
#include <stdexcept>

ClassdXdt::ClassdXdt(uword N, uword order,double dx,mat Aplus, mat Aminus,double Z_0,
                     double omega,  mat Q,  vec invMdiag,  double deltaT, double Jac, vec X0)
{
    this->N=N;
    this->order=order;
    this->dx=dx;
    this->Aplus=Aplus;
    this->Aminus=Aminus;
    this->Z_0=Z_0;
    this->omega=omega;
    this->Q=Q;
    // The nodal Gauss-Legendre basis yields a DIAGONAL mass matrix, so applying
    // M^{-1} is just a per-DOF scaling. main() now provides that diagonal
    // directly instead of building the dense kron(eye(N), inv(M)).
    this->invMdiag = invMdiag;
    if (this->invMdiag.n_elem != 2*order*N)
        throw std::invalid_argument("invMdiag has the wrong size");
    this->deltaT=deltaT;
    //this->Jac=Jac;
    this->Xend=uword(X0.size()-1);
    this->invJac=1/Jac;

    //Runge-Kutta-Coefficients
    k1 = zeros<vec>(2*order*N);
    k2 = zeros<vec>(2*order*N);
    k3 = zeros<vec>(2*order*N);
    k4 = zeros<vec>(2*order*N);
    k5 = zeros<vec>(2*order*N);
    k6 = zeros<vec>(2*order*N);
    
    psi=zeros<mat>(2*order,4);
    // Note: the volume operator is applied block-wise per cell in FluidMatrix
    // (Q * reshape(X)), so we no longer form the dense kron(eye(N),Q).
    //Q.save("Q.dat",raw_ascii);

    switch (order) {
        case 1: {
            psi.cols(2,3)=eye<mat>(2,2);
            psi.cols(0,1)=eye<mat>(2,2);
            }break;

        case 2: {
            psi.cols(2,3)= mat({{(1-sqrt(3)), 0 },  {0, (1-sqrt(3))}, {(1+sqrt(3)), 0}, {0, (1+sqrt(3))}});
            psi.cols(0,1)=mat({{(1+sqrt(3)), 0},    {0, (1+sqrt(3))},{(1-sqrt(3)), 0},{  0,  (1-sqrt(3))}});
            psi=0.5*psi;
            }break;

        case 5:
            {
            double z1=-1/3.0*sqrt(5+2*sqrt(10.0/7.0));
            double z2=-1/3.0*sqrt(5-2*sqrt(10.0/7.0));
            double z3=0.0;
            double z4=1/3.0*sqrt(5-2*sqrt(10.0/7.0));
            double z5=1/3.0*sqrt(5+2*sqrt(10.0/7.0));

            auto l1 = [z1,z2,z3,z4,z5] (double z) ->double { return (z-z2)/(z1-z2)*(z-z3)/(z1-z3)*(z-z4)/(z1-z4)*(z-z5)/(z1-z5);};
            auto l2 = [z1,z2,z3,z4,z5] (double z) ->double { return (z-z1)/(z2-z1)*(z-z3)/(z2-z3)*(z-z4)/(z2-z4)*(z-z5)/(z2-z5);};
            auto l3 = [z1,z2,z3,z4,z5] (double z) ->double { return (z-z1)/(z3-z1)*(z-z2)/(z3-z2)*(z-z4)/(z3-z4)*(z-z5)/(z3-z5);};
            auto l4 = [z1,z2,z3,z4,z5] (double z) ->double { return (z-z1)/(z4-z1)*(z-z2)/(z4-z2)*(z-z3)/(z4-z3)*(z-z5)/(z4-z5);};
            auto l5 = [z1,z2,z3,z4,z5] (double z) ->double { return (z-z1)/(z5-z1)*(z-z2)/(z5-z2)*(z-z3)/(z5-z3)*(z-z4)/(z5-z4);};


            psi.cols(2,3)= mat({{l1(1), 0 },  {0, l1(1)}, {l2(1), 0}, {0, l2(1)},  {l3(1), 0}, {0, l3(1)},  {l4(1), 0}, {0, l4(1)},  {l5(1), 0}, {0, l5(1)} });
            psi.cols(0,1)= mat({{l1(-1), 0 },  {0, l1(-1)}, {l2(-1), 0}, {0, l2(-1)},  {l3(-1), 0}, {0, l3(-1)},  {l4(-1), 0}, {0, l4(-1)},  {l5(-1), 0}, {0, l5(-1)} });
            //psi.print("psi= ");
            } break;
        default:
            {
            // Generic edge-trace operator for any order (e.g. 8). psi columns hold
            // the Lagrange basis values at the cell edges, interleaved with the 2x2
            // identity (p,u) per node: cols 0-1 = left edge (z=-1), cols 2-3 = right
            // edge (z=+1). kron(L,eye(2)) reproduces the explicit order-5 layout.
            vec Stuetz, weights;
            GaussLegendre(order, Stuetz, weights);
            psi.cols(0,1)=kron(LagrangeAtPoint(Stuetz,-1.0), eye<mat>(2,2));
            psi.cols(2,3)=kron(LagrangeAtPoint(Stuetz, 1.0), eye<mat>(2,2));
            } break;
        }
    }

ClassdXdt::~ClassdXdt() {
    //dtor
    }

vec ClassdXdt::numericalFlux(const vec &Ul, const vec &Ur) {
    //Berechnung der Zustände am Zellenrand linkes Riemann-Problem
    // right/left are LOCAL (not members) so this is reentrant -> safe to call
    // concurrently from the OpenMP-parallelized cell loop in FluidMatrix.
    vec right = trans(psi.cols(0,1))*Ur;
    vec left  = trans(psi.cols(2,3))*Ul;
    // Berechnung des Flusses
    return Aminus*right + Aplus*left;
    }


void ClassdXdt::FluidMatrix(const vec &X, vec &dXdt, const double &t) {
    //Q.print("Q=");
    //invM.print("invM=");
    vec Flux = zeros<vec>(X.size());
    switch (order) {

        case 1: {
            //Berechnung des numerischen Flusses
            //Fluss am Rand bzw Randbedingungen
            Flux.subvec(0,1)=   -psi.cols(2,3)*numericalFlux(X.subvec(0,1), X.subvec(2,3) )
                                //+ psi.cols(0,1) * numericalFlux(trans(vec({1,-1}))* X.subvec(0,1), X.subvec(0,1) ); //linker Rand
                                +psi.cols(0,1)*numericalFlux(vec({1,1/Z_0})*sin(omega*t) , X.subvec(0,1));  //linker Rand Sinus


            Flux.subvec(Xend-1,Xend)=   -psi.cols(2,3)*numericalFlux(X.subvec(Xend-1,Xend),trans(vec({1,-1})*X.subvec(Xend-1,Xend)))
                                        + psi.cols(0,1)*numericalFlux(X.subvec(Xend-3,Xend-2) , X.subvec(Xend-1,Xend));

             //Innere Zellen
            #pragma omp parallel for if(N>=1024) schedule(static)
            for (uword k=1; k<=uword(N-2); k++) {
                Flux.subvec(2*k,2*(k+1)-1) = -psi.cols(2,3) * numericalFlux( X.subvec(2*k,2*(k+1)-1) ,  X.subvec(2*(k+1), 2*(k+2)-1   ))  //rechter Rand
                                             +psi.cols(0,1) * numericalFlux( X.subvec(2*(k-1),2*k-1) ,  X.subvec(2*k,2*(k+1)-1) )  ; //linker Rand
                }
            }break;

        case 2: {
            //Berechnung des numerischen Flusses
            //Fluss am Rand bzw Randbedingungen
            Flux.subvec(0,3)=   - psi.cols(2,3) * numericalFlux(X.subvec(0,3), X.subvec(4,7))
                                //+ psi.cols(0,1) * numericalFlux(kron(mat({{0,1},{1,0}}), mat({{1,0},{0,-1}}))* X.subvec(0,3), X.subvec(0,3) );//zero Gradient Condition linker globaler Rand;
                                + psi.cols(0,1) * numericalFlux(kron(ones<vec>(2),vec({1,1/Z_0}))*sin(omega*t) , X.subvec(0,3)); //von links kommende Sinuswelle

            //cout << "Xend=" << Xend << endl;
            Flux.subvec(Xend-3,Xend)=   -psi.cols(2,3)*numericalFlux(X.subvec(Xend-3,Xend) , kron(mat({{0,1},{1,0}}), mat({{1,0},{0,-1}}))*X.subvec(Xend-3,Xend))
                                        +psi.cols(0,1)*numericalFlux(X.subvec(Xend-7,Xend-4),X.subvec(Xend-3,Xend));

            //Innere Zellen
            #pragma omp parallel for if(N>=1024) schedule(static)
            for (uword k=1; k<=uword(N-2); k++) {
                Flux.subvec(4*k,4*(k+1)-1) = -psi.cols(2,3) * numericalFlux( X.subvec(4*k,4*(k+1)-1) ,  X.subvec(4*(k+1), 4*(k+2)-1   ))  //rechter Rand
                                             +psi.cols(0,1) * numericalFlux( X.subvec(4*(k-1),4*k-1) ,  X.subvec(4*k,4*(k+1)-1) )  ; //linker Rand
                }
            } break;

        default:
            {
            //Berechnung des numerischen Flusses
            //Fluss am Rand bzw Randbedingungen
            Flux.subvec(0,2*order-1)=   - psi.cols(2,3) * numericalFlux(X.subvec(0,2*order-1), X.subvec(2*order,4*order-1))
                                        //+ psi.cols(0,1) * numericalFlux(kron(fliplr(eye<mat>(order,order)), mat({{1,0},{0,-1}}))* X.subvec(0,2*order-1), X.subvec(0,2*order-1) );//zero Gradient Condition linker globaler Rand;
                                        + psi.cols(0,1) * numericalFlux(kron(ones<vec>(order),vec({1,1/Z_0}))*sin(omega*t) , X.subvec(0,2*order-1)); //von links kommende Sinuswelle

            //cout << "Xend=" << Xend << endl;
            Flux.subvec(Xend-2*order+1,Xend) =  - psi.cols(2,3) * numericalFlux(
                                                X.subvec(Xend-2*order+1,Xend) , kron(fliplr(eye<mat>(order,order)),
                                                mat({{1,0},{0,-1}}))*X.subvec(Xend-2*order+1,Xend))
                                                + psi.cols(0,1)*numericalFlux(X.subvec(Xend-4*order+1,Xend-2*order),X.subvec(Xend-2*order+1,Xend));

            //Innere Zellen
            #pragma omp parallel for if(N>=1024) schedule(static)
            for (uword k=1; k<=uword(N-2); k++) {
                Flux.subvec(2*order*k,2*order*(k+1)-1) =    -psi.cols(2,3) * numericalFlux( X.subvec(2*order*k,2*order*(k+1)-1) ,  X.subvec(2*order*(k+1), 2*order*(k+2)-1   ))  //rechter Rand
                                                            +psi.cols(0,1) * numericalFlux( X.subvec(2*order*(k-1),2*order*k-1) ,  X.subvec(2*order*k,2*order*(k+1)-1) )  ; //linker Rand
            }
            }break;
    }
    // Volume term applied block-wise: A = kron(eye(N),Q) is block-diagonal, so
    // A*X == Q applied to each cell's 2*order-block. Reshaping X to (2*order)xN
    // turns it into one small GEMM (~N x less work than the dense 500x500 matvec),
    // which Accelerate parallelizes across cells for large N.
    dXdt = invMdiag % (vectorise(Q * reshape(X, 2*order, N)) + Flux);
    dXdt*=invJac;
    //Flux.print("Flux=");
    //X.print("X= ");
    //dXdt.print("dXdt= ");
    //Flux.print("Flux =");
    return;
    }

void ClassdXdt::EulerExplicitClassic(vec &X, double &t) {
    //k1.fill(0.0);
    //k1.print("k1= ");
    FluidMatrix(X,k1,t);
    //X=X+deltaT*k1;
    X+=deltaT*k1;
    //k1.print("k1= ");
    t+=deltaT;
    }

void ClassdXdt::RungeKuttaIntegrationClassic(vec &X, double &t) {
    //k1.print("k1=");
    //k2.print("k2=");
    //k3.print("k3=");
    //k4.print("k4=");
    //vec X1=X;
    FluidMatrix(X,k1,t);
    //k1.print("k1=");
    t+= deltaT/2.0;
    FluidMatrix(X+deltaT/2.0*k1,k2,t);
    FluidMatrix(X+deltaT/2.0*k2,k3,t);

    t+= deltaT/2.0;
    FluidMatrix(X+deltaT*k3,k4,t);

    /*
    k1.print("k1=");
    k2.print("k2=");
    k3.print("k3=");
    k4.print("k4=");
    */

    X+=deltaT/6.0*(k1+2.0*(k2+k3)+k4);
    //vec Xdiff=X1-X;
    //Xdiff.print("Xdiff= ");
    //X.print("X= ");
    }

void ClassdXdt::RungeKuttaSecondOrderHeun(vec &X, double &t) {
    FluidMatrix(X,k1,t);
    t+=deltaT;
    FluidMatrix(X+deltaT*k1,k2,t);
    X+=deltaT/2.0*(k1+k2);

}


void ClassdXdt::RungeKuttaIntegrationSecondVersion(vec &X, double &t){
    FluidMatrix(X,k1,t);
    t+=deltaT/3.0;
    FluidMatrix(X+deltaT/3.0*k1,k2,t);
    t+=deltaT/3.0;
    FluidMatrix(X-deltaT/3*k1+deltaT*k2,k3,t);
    t+=deltaT/3.0;
    FluidMatrix(X+deltaT*(k1-k2+k3),k4,t);
    X+=deltaT/8.0*(k1+3*(k2+k3)+k4);
}

void ClassdXdt::ButcherIntegration(vec &X, double &t){
//stub
}

void ClassdXdt::PredictorCorrectorHeun(vec &X, double &t){   //Heun-Verfahren
    epsilon =1;
    this->Xold=X;
    FluidMatrix(X,k1,t);
    this->RungeKuttaIntegrationClassic(X,t);

    do {
        //cout << &Xold << endl << &X <<endl;
        //Predictor
        this->Xpredictor=X;
        //Corrector
        FluidMatrix(X,k2,t);
        //New timestep
        X=this->Xold+deltaT/2.0*(k1+k2);
        epsilon = sum(abs(X-Xpredictor));

    } while (epsilon >= 1e-12);
    cout << "Convergence occured in Time : " << t-deltaT << endl;
}

    /*
vec CheckConvergence(){

    return
}
*/






