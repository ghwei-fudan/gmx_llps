//
// Created by Yiming Tang on 2018/10/6.
//

#ifndef GMX_ASSEMBLY_ASSEMBLY_H
#define GMX_ASSEMBLY_ASSEMBLY_H

#include <gromacs/trajectoryanalysis.h>
#include <iostream>

using namespace gmx;

class assembly : gmx::TrajectoryAnalysisModule {
public:
    assembly();

    virtual void initOptions(IOptionsContainer *options,
                             TrajectoryAnalysisSettings *settings);

    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation &top);

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);

    virtual void finishAnalysis(int nframes);

    virtual void writeOutput();

private:

    class ModuleData;

    // File names for output control

    // Each line of this output file is the number of clusters
    // in a particular frame.
    std::string fnClusterCount_;

    // Each line of this output file is the number of molecules
    // in the largest cluster in a particular frame.
    std::string fnLargestCluster_;

    // Input parameters

    // Cutoff for smallest number of molecules in a "Cluster"
    double cutoffClusterSize_;

    // Cutoff for intra-cluster interactions
    double cutoffSpace_;

    // Global parameters used in software

    // Number of molecules in this system
    int molCount_;

    Selection sel_;
    AnalysisNeighborhood nb_;

    AnalysisData dataClusterCount_;
    AnalysisData dataLargestCluster_;

    t_topology *top_;
    t_atoms atoms_;

};


#endif //GMX_ASSEMBLY_ASSEMBLY_H
