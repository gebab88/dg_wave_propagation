#ifndef CLASSDXDT_H
#define CLASSDXDT_H

#include <string>
#include <armadillo>
#include <cmath>

using namespace arma;
class ClassdXdt
{
    public:
        ClassdXdt(uword N, uword order,double dx,mat Aplus, mat Aminus,double Z_0,
                     double omega,  mat Q,  vec invMdiag,  double deltaT, double Jac, vec X0,
                     uword OpenMPMinCells = 50, uword OpenMPThreads = 2);

        ~ClassdXdt();

	//Attributes:
        double dx,Z_0, omega,deltaT,invJac;
        double epsilon;
        uword N,order,OpenMPMinCells,OpenMPThreads;
        mat Q;
        mat Aplus, Aminus,psi,psiLeft,psiRight,rightWallReflection;
        vec leftInflowState;
        vec invMdiag;   // diagonal of the (block-)inverse mass matrix, per DOF
        uword Xend;
        vec k1,k2,k3,k4,k5,k6;
        vec Xold,Xpredictor;

        //Member Functions:
        void RungeKuttaIntegrationClassic(vec &X, double &t);
        void EulerExplicitClassic(vec &X, double &t);
        void ButcherIntegration(vec &X, double &t);
        void RungeKuttaGill(vec &X, double &t);
        void RungeKuttaIntegrationSecondVersion(vec &X, double &t);
        void PredictorCorrectorHeun(vec &X, double &t);
        void RungeKuttaSecondOrderHeun(vec &X,double &t);

    protected:

    private:
        vec numericalFlux(const vec &Ul, const vec &Ur);
        void FluidMatrix(const vec &X, vec &dXdt, const double &t);

};

#endif // CLASSDXDT_H
