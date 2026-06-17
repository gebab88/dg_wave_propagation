#include <iostream>
#include <string>
#include <cmath>
#include <armadillo>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

#include "gnuplot-iostream.h"

#include "InitialAndBoundary.hpp"
#include "ClassdXdt.hpp"
#include "vandermonde.hpp"
#include "ClassPlot.hpp"

using namespace arma;
using namespace std;


int main( int argc, char *argv[] ){
    chrono::time_point<chrono::system_clock> StartTime, EndTimeCalc,EndTimePlot;
    StartTime = chrono::system_clock::now();

    #include "config.hpp"

    float lambda_1 = c_0;
    float lambda_2 = -c_0;

    mat R = {{Z_0 , -Z_0}, {1.0 ,1.0}};
    mat R_inv= { {pow(2 * Z_0, -1) ,  0.5}, { (-1) * pow(2 * Z_0, -1) ,  0.5}};
    mat R_inv_calc=inv(R);

    mat lambda_plus = {{lambda_1, 0.0},{0.0 ,0.0}};

    mat lambda_minus ={{0.0 ,  0.0},{0.0 ,  lambda_2}};

    // Jacobi Matrices
    mat Aplus = R * lambda_plus*R_inv;
    //Aplus.print("Aplus= ");

    mat Aminus = R * lambda_minus *R_inv;
    //Aminus.print("Aminus= ");

    mat Q,invM;
    mat U0;
    mat A=eye<mat>(2,2);

    vec xnodes={ };

    #if ZEITODERFREQ == 1
    {
    InitialAndBoundary InitialBoundaryObj(height,Z_0,mu,sigma);
    U0 = zeros<mat>(2,order*N);
    switch (order)
    {
    case 1:
        {
            vec ui(2,fill::zeros);
            vec xi(order,fill::zeros);
            for (uword i=0; i<N; i++) {
                vec xi={0.5*(x[i]+x[i+1])};
                xnodes=join_cols(xnodes,xi);
                InitialBoundaryObj.U_init(ui,order); //Initialisierung mit 0
                //InitialBoundaryObj.U_init(xi,ui); //Initialisierung mit GaussGlocke
                //ui.print("ui =");
                U0.col(i)=ui;
            }
        //U0.print("U0= ");
        }break;

    case 2:
        {
        mat ui(2,2,fill::zeros);
        vec xi(order,fill::zeros);
        vec Stuetz={-sqrt(1.0/3.0) ,  sqrt(1.0/3.0)};
        for (uword i=0; i<N; i++) {
            //Stuetz.print("Stuetz=");
            vec xi=0.5*(x[i]*ones(order)+x[i+1]*ones(order) + dx * Stuetz);
            //xi.print("xi=");
            xnodes=join_cols(xnodes,xi);
            InitialBoundaryObj.U_init(ui,order); //Initialisierung mit 0
            //InitialBoundaryObj.U_init(xi,ui); //Initialisierung mit GaussGlocke
            //ui.print("ui=");
            U0.cols(i*order,i*order+1)=ui;

            //ui.print("ui =");
            //ui.save("ui.dat",raw_ascii);
            //a.print("a=");
            //b.print("b=");
            //Stuetz.print("Stütz=");
            //xi.print("xi: ");
        }
            // Aufbau der Massenmatrix
            mat M=eye<mat>(4,4);
            //mat invM=inv(M); //Inverse der Massenmatrix
            invM=M; //Inverse der Massenmatrix für 2 Stuetzstellen Einheitsmatrix (assign to outer invM, no shadowing)
            invM=kron(eye<mat>(N,N),invM); //Expansion auf alle Zellen
            //invM.print("invM= ");
            //invM.save("invM.dat",raw_ascii);
//
//          Berechnung des Flussintegrals, Steifigkeitsmatrix
            Q = kron(mat({{-1, -1},{ 1, 1}}), eye<mat>(2,2));
            Q=sqrt(3.0)/2.0 * Q;

            mat Funktionalmat={{0 ,  Z_0*c_0},{pow(rho_0,-1), 0}};
            Funktionalmat=kron(eye<mat>(2,2),Funktionalmat);

            Q=Funktionalmat*Q;
            Q.print("Q=");
            //Q.save("Q.dat",raw_ascii);
            } break;
    case 5:
        {
        mat ui(2,order,fill::zeros);
        vec xi(order,fill::zeros);

        double z1=-1.0/3.0*sqrt(5+2*sqrt(10.0/7.0));
        double z2=-1.0/3.0*sqrt(5-2*sqrt(10.0/7.0));
        double z3=0.0;
        double z4=1.0/3.0*sqrt(5-2*sqrt(10.0/7.0));
        double z5=1.0/3.0*sqrt(5+2*sqrt(10.0/7.0));

        vec Stuetz={z1,z2,z3,z4,z5};
        //Stuetz.print("z= ");

        for (uword i=0; i<N; i++) {
            //Stuetz.print("Stuetz=");
            vec xi=0.5*(x[i]*ones(order)+x[i+1]*ones(order) + dx * Stuetz);
            //xi.print("xi=");
            xnodes=join_cols(xnodes,xi);
            InitialBoundaryObj.U_init(ui,order); //Initialisierung mit 0
            //InitialBoundaryObj.U_init(xi,ui); //Initialisierung mit GaussGlocke
            //ui.print("ui=");
            U0.cols(i*order,(i+1)*order-1)=ui;

            //ui.print("ui =");
            //ui.save("ui.dat",raw_ascii);
            //a.print("a=");
            //b.print("b=");
            //Stuetz.print("Stütz=");
            //xi.print("xi: ");
        }

//          // Aufbau der Massenmatrix
            vec alpha(5);
            alpha[0]= (322-13*sqrt(70.0))/900.0;
            alpha[1]= (322+13*sqrt(70))/900.0;
            alpha[2]=128.0/225.0;
            alpha[3]=alpha[1];
            alpha[4]=alpha[0];
            //alpha.print("alpha = " );

            vec Mdiag=kron(eye<mat>(order,order),ones<vec>(2))*alpha;
            mat M=diagmat(Mdiag);
            //M.print("M= ");
            invM=inv(M); //Inverse der Massenmatrix (assign to outer invM, no shadowing)
            invM=kron(eye<mat>(N,N),invM); //Expansion auf alle Zellen
            //invM.print("invM= ");
            //invM.save("invM.dat",raw_ascii);
//
//          Berechnung des Flussintegrals, Steifigkeitsmatrix
            // Berechnung der transpornierten Differentationsmatrix
            mat DT=TransposeDiffMatrix(order,Stuetz);
            DT=kron(DT,eye<mat>(2,2));
            // Berechnung der Steifigkeitsmatrix
            Q=DT*M;

            mat Funktionalmat={{0 ,  Z_0*c_0},{pow(rho_0,-1), 0}};
            Funktionalmat=kron(eye<mat>(order,order),Funktionalmat);

            Q=Funktionalmat*Q;
            //Q.print("Q= ");
        } break;

    default:
        {
        // General order (e.g. 8): Gauss-Legendre nodes/weights computed on the
        // fly, then the same construction as case 5 (diagonal mass matrix from
        // the weights, stiffness Q = Funktional * D^T M). Mirrors case 5 exactly
        // but parametrized by `order`.
        mat ui(2,order,fill::zeros);

        vec Stuetz, alpha;
        GaussLegendre(order, Stuetz, alpha);   // nodes & weights on [-1,1]

        for (uword i=0; i<N; i++) {
            vec xi=0.5*(x[i]*ones(order)+x[i+1]*ones(order) + dx * Stuetz);
            xnodes=join_cols(xnodes,xi);
            InitialBoundaryObj.U_init(ui,order); //Initialisierung mit 0
            U0.cols(i*order,(i+1)*order-1)=ui;
        }

        // Diagonal mass matrix: each node's weight is shared by both fields (p,u).
        vec Mdiag=kron(alpha,ones<vec>(2));
        mat M=diagmat(Mdiag);
        invM=inv(M);
        invM=kron(eye<mat>(N,N),invM);         //Expansion auf alle Zellen

        // Stiffness matrix Q = Funktionalmat * D^T * M, kron'd to the 2-field system.
        mat DT=TransposeDiffMatrix(order,Stuetz);
        DT=kron(DT,eye<mat>(2,2));
        Q=DT*M;

        mat Funktionalmat={{0 ,  Z_0*c_0},{pow(rho_0,-1), 0}};
        Funktionalmat=kron(eye<mat>(order,order),Funktionalmat);
        Q=Funktionalmat*Q;
        } break;

    }
    }    //delete InitialBoundaryObj and other Objects

    //U0.save("U0.dat",raw_ascii);
    U0=vectorise(U0);
    //U0=zeros<vec>(U0.size()); // zu Nullsetzen von U0

    // Define ODE
    // Order-aware CFL: high-order DG needs dt ~ dx/(c*(2p+1)) with p = order-1,
    // so the stable timestep shrinks like 1/(2*order-1). Folding that in here
    // makes CFLReserve an order-independent safety fraction of the CFL limit.
    const double deltaT = CFLReserve * dx / (abs(lambda_1) * (2.0*order - 1.0));
    //U0.print("U0=");
    //invM.print("invM=");
    //cout<<t1<<endl;
    ClassdXdt dXdtObj(N,order,dx,Aplus,Aminus,Z_0,omega,Q,invM,deltaT, Jac,U0); //Initialization Objekt FundamentalMatrix
    vec X=U0;
  //dXdtObj.FluidMatrix(U0, dXdt0, t1);
    double t= t0;
    //Integration with RungeKutta
    vec TimeSteps={t};
    mat solution={trans(X)};
    vec TimeStep(1);

    X=U0;
    t=t0;
    TimeSteps={t};
    solution={trans(X)};

    if (TimeIntegration=="EulerExplicit"){
        do {
            dXdtObj.EulerExplicitClassic(X, t);
            TimeStep[0]=t;
            solution=join_vert(solution,trans(X));
            TimeSteps=join_vert(TimeSteps,TimeStep);
            //X.print("X= ");
            cout<<"t=" << t<<endl;
            } while (t< t1);
    }
    else if(TimeIntegration=="RungeKuttaClassic"){
        do{
            //vec X1=X;
            dXdtObj.RungeKuttaIntegrationClassic(X, t);
            TimeStep[0]=t;
            //vec Xdiff=X1-X;
            //Xdiff.print("Xdiff= ");
            solution=join_vert(solution,trans(X));
            TimeSteps=join_vert(TimeSteps,TimeStep);
            //cout<<"t=" << t<<endl;
            //printf("t= %f0 \n",t);
        } while (t<t1);
    }
    else if(TimeIntegration=="RungeKuttaSecond") {
        do {
            //vec X1=X;
            dXdtObj.RungeKuttaIntegrationSecondVersion(X, t);
            TimeStep[0] = t;
            //vec Xdiff=X1-X;
            //Xdiff.print("Xdiff= ");
            solution = join_vert(solution, trans(X));
            TimeSteps = join_vert(TimeSteps, TimeStep);
            //printf("t= %f0 \n",t);
        } while (t < t1);
        //cout<<"Sorry, this is still a stub"<<endl;
        //exit(EXIT_FAILURE);
    }

    else if(TimeIntegration=="RungeKuttaSecondOrder"){
        do{
                dXdtObj.RungeKuttaSecondOrderHeun(X,t);
                TimeStep[0]=t;
                solution=join_vert(solution,trans(X));
                TimeSteps=join_vert(TimeSteps,TimeStep);
            } while  (t<t1);
    }

    else if(TimeIntegration=="PredictorCorrectorHeun"){
        do{
            //vec X1=X;
            dXdtObj.PredictorCorrectorHeun(X, t);
            TimeStep[0]=t;
            solution=join_vert(solution,trans(X));
            TimeSteps=join_vert(TimeSteps,TimeStep);
            //cout<<"t=" << t<<endl;
            //printf("t= %f0 \n",t);
        } while (t<t1);
    }
    else{
        cout << "The TimeIntegration variable isn't defined correctly!" << endl;
        exit(EXIT_FAILURE);
    }

    //X.save("X.dat",raw_ascii);
    //solution.print("Solution =");
    //solution.save("X.dat",raw_ascii);
    //TimeSteps.save("time.dat",raw_ascii);

    EndTimeCalc = chrono::system_clock::now();


    double elapsed_time =  chrono::duration_cast<chrono::microseconds>
                             (EndTimeCalc-StartTime).count();



    cout << "The calculations are finished!" << endl;
    cout << "elapsed time: " << elapsed_time / 1000.0 /1000.0 << "s\n";
    
    cout << "Plotting starts!" << endl;

    ClassPlot PlotObj(solution,TimeSteps,xnodes);

    PlotObj.Plot();

    cout << "Encoding video!" << endl;
    PlotObj.MakeVideo("output.mp4");

    EndTimePlot = chrono::system_clock::now();
    elapsed_time =  chrono::duration_cast<chrono::microseconds>
                             (EndTimePlot-EndTimeCalc).count();
    cout << "Plotting has finished!" <<endl;
    cout << "elapsed time: " << elapsed_time / 1000.0 /1000.0 << "s\n";

    cout << "Press enter to exit." << endl;
	cin.get();

    #else

    #endif  //ifdef ZEITINT endif

    return 0;
}
