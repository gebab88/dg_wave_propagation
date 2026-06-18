#include "ClassPlot.hpp"
#include <string>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <thread>
#include <vector>

ClassPlot::ClassPlot(mat pressure,vec TimeSteps, const vec &xnodes)
{
    this->pressure=pressure;
    this->TimeSteps=TimeSteps;
    this->xnodes=xnodes;
}

ClassPlot::~ClassPlot()
{
    //dtor
}

void ClassPlot::Plot()
{
    const uword nframes = TimeSteps.size();
    if (nframes == 0) return;

    if (pressure.n_rows != nframes || pressure.n_cols != xnodes.size()) {
        std::cout << "Plotting failed: pressure data has shape "
                  << pressure.n_rows << "x" << pressure.n_cols
                  << ", expected " << nframes << "x" << xnodes.size()
                  << std::endl;
        return;
    }

    unsigned nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 4;
    if (nthreads > nframes) nthreads = static_cast<unsigned>(nframes);
    std::cout << "Plotting " << nframes << " frames on " << nthreads
              << " threads" << std::endl;

    // Each worker writes ONE self-contained gnuplot script (with inline data)
    // for its strided share of frames (tid, tid+nthreads, ...), then runs
    // gnuplot on that file via system(). Driving gnuplot through a script file
    // rather than a popen pipe avoids the macOS concurrent-popen deadlock, where
    // sibling gnuplot processes inherit each other's pipe write-ends (no
    // close-on-exec) and never see EOF. Frames are disjoint, so every PNG is
    // written exactly once.
    auto scriptPath = [](unsigned tid) {
        return "/tmp/dg_plot_" + std::to_string(tid) + ".gp";
    };

    // Phase 1 (parallel, in-process): each thread writes one self-contained
    // gnuplot script with inline data for its strided share of frames. This is
    // the CPU-bound part and scales well across threads.
    auto builder = [&](unsigned tid)
    {
        std::ofstream gp(scriptPath(tid));
        gp << std::scientific << std::setprecision(10);
        gp << "set terminal png size 800,600\n";
        gp << "set grid\n";
        gp << "set xrange [0:0.5]\nset yrange [-2:2]\n";

        std::string str;
        for (uword i = tid; i < nframes; i += nthreads)
        {
            vec p = trans(pressure.row(i));

            str = std::to_string(i);
            str.insert(str.begin(), 8 - str.length(), '0');
            gp << "set output '" << str << ".png'\n";
            gp << "plot '-' with linespoints ls 1\n";
            for (uword k = 0; k < xnodes.size(); ++k)
                gp << xnodes(k) << " " << p(k) << "\n";
            gp << "e\n";
        }
    };

    {
        std::vector<std::thread> pool;
        pool.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t)
            pool.emplace_back(builder, t);
        for (auto &th : pool)
            th.join();
    }

    // Phase 2: render all scripts concurrently. We launch the whole batch from a
    // SINGLE shell (`gnuplot a & gnuplot b & ... ; wait`) rather than calling
    // system() from each worker thread. Forking repeatedly from this heavily
    // multi-threaded process is pathologically slow on macOS; forking once into
    // a small single-threaded shell that then spawns the gnuplot fleet is ~10x
    // faster and gives the expected parallel speedup.
    std::string cmd;
    for (unsigned t = 0; t < nthreads; ++t)
        cmd += "gnuplot '" + scriptPath(t) + "' & ";
    cmd += "wait";
    int ret = std::system(cmd.c_str());
    if (ret != 0)
        std::cout << "Plotting failed (gnuplot returned " << ret
                  << "). Is gnuplot on the PATH?" << std::endl;

    for (unsigned t = 0; t < nthreads; ++t)
        std::remove(scriptPath(t).c_str());
}

void ClassPlot::MakeVideo(const std::string &filename)
{
    // Encode the %08d.png frames written by Plot() into a video via ffmpeg.
    // 60 input fps -> 30 output fps; yuv420p keeps the result widely playable.
    std::string cmd = "ffmpeg -y -framerate 60 -i %08d.png "
                      "-c:v libx264 -pix_fmt yuv420p -r 30 " + filename;
    std::cout << "Encoding video: " << cmd << std::endl;
    int ret = std::system(cmd.c_str());
    if (ret != 0)
        std::cout << "Video encoding failed (ffmpeg returned " << ret
                  << "). Is ffmpeg on the PATH?" << std::endl;
    else
        std::cout << "Video written to " << filename << std::endl;
}
