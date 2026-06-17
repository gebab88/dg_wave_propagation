#ifndef INITIALANDBOUNDARY_H
#define INITIALANDBOUNDARY_H

#include <string>
#include <armadillo>

using namespace arma;

class InitialAndBoundary
{
    public:
        InitialAndBoundary(double height, double Z_0, double mu, double sigma, std::string BoundaryRightEnd);
        InitialAndBoundary(double height, double Z_0, double mu, double sigma);
        ~InitialAndBoundary();
        double height;
        double Z_0;
        double mu;
        double sigma;
        std::string BoundaryRightEnd;

        void U_init(const vec &xi, mat &ui); //initialisierung ohne Null
        void U_init(mat &ui, const uword &order); // Initialisierung mit Null
        static double BoundaryInitialization(vec Ue, vec Uw, std::string inflow);

        void pAndu_init(vec &p_init, vec &u_init, const vec &xi);
        void gaussian(vec &gaussian,const vec &xi);


    protected:

    private:
};

#endif // INITIALANDBOUNDARY_H
