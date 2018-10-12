//
// Created by Yiming Tang on 2018/10/6.
//

#ifndef GMX_ASSEMBLY_ASSEMBLY_H
#define GMX_ASSEMBLY_ASSEMBLY_H

#include <gromacs/trajectoryanalysis.h>
#include <iostream>
#include <vector>

using namespace gmx;

class assembly : public gmx::TrajectoryAnalysisModule {
public:
    assembly();

    virtual void initOptions(IOptionsContainer *options,
                             TrajectoryAnalysisSettings *settings) override;

    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation &top) override;

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata) override;

    virtual void finishAnalysis(int nframes) override;

    virtual void writeOutput() override;

private:

    class ModuleData;

    // File names for output control

    // Each line of this output file is the number of clusters
    // in a particular frame.
    std::string fnClusterCount_;

    // Each line of this output file is the number of molecules
    // in the largest cluster in a particular frame.
    std::string fnLargestCluster_;

    // This is a PDB File containing the largest cluster.
    std::string fnLargestClusterPDB_;

    // Input parameters

    // Cutoff for smallest number of molecules in a "Cluster"
    double cutoffClusterSize_;

    // Cutoff for intra-cluster interactions
    double cutoffSpace_;

    double time_pdb_;

    // Global parameters used in software

    // Number of molecules in this system
    int molCount_;

    // Molecule ID -> atom ID mapping

    std::vector<std::vector<int>> idMap;

    Selection sel_;
    AnalysisNeighborhood nb_;

    AnalysisData dataClusterCount_;
    AnalysisData dataLargestCluster_;

    t_topology *top_;
    t_atoms atoms_;

};


#endif //GMX_ASSEMBLY_ASSEMBLY_H
