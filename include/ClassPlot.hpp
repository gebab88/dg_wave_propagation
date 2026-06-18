#ifndef CLASSPLOT_H
#define CLASSPLOT_H


#include <armadillo>
#include <string>
using namespace arma;

class ClassPlot
{
    public:
        ClassPlot(mat pressure, vec TimeSteps, const vec &xnodes);
        ~ClassPlot();
        void Plot();
        void MakeVideo(const std::string &filename = "output.mp4");

    protected:

    private:
        mat pressure;
        vec TimeSteps;
        vec xnodes;
};

#endif // CLASSPLOT_H
