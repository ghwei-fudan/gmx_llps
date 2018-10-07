#include <iostream>
#include <gromacs/trajectoryanalysis.h>
#include "assembly.h"

int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<assembly>(argc,argv);
}