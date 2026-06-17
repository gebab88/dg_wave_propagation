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
    DT.print("DT=");
    return DT;
}
