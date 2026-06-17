#ifndef CONFIGDATA
#define CONFIGDATA

const std::string scheme = "DG";
const std::string KeyOscillator = "without";

const std::string TimeIntegration = TIMEINTSCHEME;
/*
TimeIntegration = "EulerExplicit"; //
TimeIntegration ="RungeKuttaClassic"; //
//TimeIntegration ="RungeKuttaSecondOrder";
TimeIntegration="PredictorCorrectorHeun";
*/
#undef TIMEINTSCHEME

// Order-independent fraction of the CFL limit (deltaT already scales the
// order-dependent 1/(2*order-1) factor in main.cpp). 0.45 is stable across
// orders; for order 5 it reproduces the verified dt = 0.05*dx/c.
float CFLReserve=0.45;

const uword N=50; //Anzahl der Zellen
const uword Npoints=N+1;

const uword order = ORDER;
#undef ORDER

//order=2; //Spatial Order Fluid Simulation
//order=1; //Spatial Order Fluid Simulation
//order=5;

const double t0=0.0,t1=2e-02;
const float x0=0.0,x1=0.5;

vec x =linspace<vec>(x0,x1,Npoints);
const double dx=x[1]-x[0];
const double Jac=dx/2.0;


// Define problem
//Parameters of fluid
const float c_0=300.0;          //Speed of Sound
const float rho_0=1.5;          //Density Air
const float Z_0=c_0*rho_0;      //Impedance

const double freq_i=3000.0;
const double omega=2*M_PI*freq_i;

//Parameters for initial Distribution of p
const double mu=0.1;
const double sigma=0.02;
const double height=1.0;

#endif // CONFIGDATA
