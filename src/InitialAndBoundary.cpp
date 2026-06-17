#include "InitialAndBoundary.h"
#include <math.h>

InitialAndBoundary::InitialAndBoundary(double height, double Z_0, double mu, double sigma, std::string BoundaryRightEnd)
{
    this->height=height;
    this->Z_0=Z_0;
    this->mu=mu;
    this->sigma=sigma;
    this->BoundaryRightEnd=BoundaryRightEnd;
}

InitialAndBoundary::InitialAndBoundary(double height, double Z_0, double mu, double sigma)
{
    this->height=height;
    this->Z_0=Z_0;
    this->mu=mu;
    this->sigma=sigma;
    return;
}

InitialAndBoundary::~InitialAndBoundary()
{
    cout << "Calling dtor InitialAndBoundary Class" << endl;
}

void InitialAndBoundary::U_init(mat &ui, const uword &order){
    ui=zeros<mat>(2,order);
    return;
}


void InitialAndBoundary::U_init(const vec &xi, mat &ui)  //stub
{
    //mat U_init=zeros<mat>(2,size(xi));
    vec p_init=zeros<vec>(size(xi));
    vec u_init=zeros<vec>(size(xi));

    this->pAndu_init(p_init,u_init,xi);
    //U_init={{p_init(xi,height,mu,sigma)},{u_init(xi,height,mu,sigma,Z_0)}};
    ui=join_cols(trans(p_init),trans(u_init));
}

void InitialAndBoundary::pAndu_init(vec &p_init, vec &u_init, const vec &xi){
        gaussian(p_init,xi);
        u_init=p_init/Z_0;
        return;
    }

void InitialAndBoundary::gaussian(vec &gaussian,const vec &xi){
        gaussian = height*exp((-1.0)*pow(xi - mu, 2.) / (2 * pow(sigma, 2)));
        return;
    }




