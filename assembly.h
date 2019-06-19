//
// Created by Yiming Tang on 2018/10/6.
//

#ifndef GMX_ASSEMBLY_ASSEMBLY_H
#define GMX_ASSEMBLY_ASSEMBLY_H

#include <gromacs/trajectoryanalysis.h>
#include <iostream>
#include <vector>
#include <set>

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

    // Each lien of this output file is the number of molecules that
    // are included in clusters in a particular frame.
    std::string fnMoleculesInCluster_;

    // This is a PDB File containing the largest cluster.
    std::string fnLargestClusterPDB_;

    // This file contains the liquidity of each frame
    // defined by fraction of aggregated molecules in last frame that remain aggregated in this frame.
    std::string fnLiquidity_;

    // These variables are for use of liquidity only, which only work under serial operations.

    std::set<int> cluster_last_frame_;
    float time_last_frame_;

    // This file contains density for aggregated and dispersed phases

    std::string fnDensity_;

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
    AnalysisNeighborhood nbDensity_;
    double density_nb_multiplier_;

    AnalysisData dataClusterCount_;
    AnalysisData dataLargestCluster_;
    AnalysisData dataMoleculesInCluster_;
    AnalysisData dataLiquidity_;
    AnalysisData dataDensity_;

    t_topology *top_;
    t_atoms atoms_;

};


#endif //GMX_ASSEMBLY_ASSEMBLY_H
