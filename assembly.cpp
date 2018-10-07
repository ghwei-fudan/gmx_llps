//
// Created by Yiming Tang on 2018/10/6.
//

#include "assembly.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace gmx;


assembly::assembly() : cutoffSpace_(0.40)
{
    registerAnalysisDataset(&dataClusterCount_, "cluster count");
    registerAnalysisDataset(&dataLargestCluster_, "largest cluster size");
}

void assembly::initOptions(IOptionsContainer          *options,
                           TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
            {
                    "Assembly Analyzer by Yiming Tang @ Fudan.\n",
                    "For now the only function of this program is to read in a trajectory ",
                    "and show the number of clusters, and the size of the largest cluster ",
                    "as a function of simulation time."
            };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("nb")
                               .filetype(eftPlot).outputFile().required()
                               .store(&fnClusterCount_).defaultBasename("number")
                               .description("Number of clusters as a function of time"));

    options->addOption(FileNameOption("sz")
                               .filetype(eftPlot).outputFile().required()
                               .store(&fnLargestCluster_).defaultBasename("size")
                               .description("Size of the largest cluster as a function of time"));

    options->addOption(SelectionOption("select")
                               .store(&sel_).required()
                               .description("Group that contain your molecules"));

    options->addOption(DoubleOption("cutoff_cz").store(&cutoffClusterSize_)
                               .defaultValue(20).required()
                               .description("Minimum number of molecules in a so-called cluster"));

    options->addOption(DoubleOption("cutoff_space").store(&cutoffSpace_)
                               .defaultValue(0.4).required()
                               .description("Maximum atom-wise distance to insure an interaction"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void assembly::initAnalysis(const TrajectoryAnalysisSettings &settings,
                            const TopologyInformation &top) {
    nb_.setCutoff(cutoffSpace_);
    dataLargestCluster_.setColumnCount(0, 1);
    dataClusterCount_.setColumnCount(0, 1);

    if (!fnClusterCount_.empty()) {
        AnalysisDataPlotModulePointer plotClusterCount(
                new AnalysisDataPlotModule(settings.plotSettings())
        );
        plotClusterCount->setFileName(fnClusterCount_);
        plotClusterCount->setTitle("Number of Clusters");
        plotClusterCount->setXAxisIsTime();
        plotClusterCount->setYLabel("Number of Clusters");
        dataClusterCount_.addModule(plotClusterCount);
    }

    if (!fnLargestCluster_.empty()) {
        AnalysisDataPlotModulePointer plotLargestCluster(
                new AnalysisDataPlotModule(settings.plotSettings())
        );
        plotLargestCluster->setFileName(fnLargestCluster_);
        plotLargestCluster->setTitle("Size of Largest Cluster");
        plotLargestCluster->setXAxisIsTime();
        plotLargestCluster->setYLabel("Size");
        dataLargestCluster_.addModule(plotLargestCluster);
    }

    this->top_ = top.topology();
    this->atoms_ = top.topology()->atoms;

    sel_.initOriginalIdsToGroup(top_, INDEX_MOL);
    molCount_ = sel_.mappedIds()[sel_.atomCount() - 1];
    cout << "Total number of molecules in this system: " << molCount_ << endl;

    // Now we are going to generate "molecule ID -> atom ID" map

    this->idMap = new vector<int>[molCount_];

    for(int i = 0; i < sel_.atomCount(); i++)
    {
        this->idMap[sel_.mappedIds()[i]].insert(this->idMap[sel_.mappedIds()[i]].end(), i);
    }

}

void assembly::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                            TrajectoryAnalysisModuleData *pdata) {
    AnalysisDataHandle dhClusterCount = pdata->dataHandle(dataClusterCount_);
    AnalysisDataHandle dhClusterSize = pdata->dataHandle(dataLargestCluster_);
    const Selection &sel = sel_;

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel_);

    dhClusterCount.startFrame(frnr, fr.time);
    dhClusterSize.startFrame(frnr, fr.time);

    // We first create two vector containing molecules that have been / not have been mapped.
    vector<int> mappedMolecule;
    vector<int> unmappedMolecule;
    for (int i = 0; i < molCount_; i++) {
        mappedMolecule.insert(unmappedMolecule.end(), i);
    }

    // We also create a cluster-vector containing vectors of molecule ids in this cluster.
    vector<vector<int>> clusterList;

    // For each molecule in unmappedMolecule, we build up a cluster based on it.
    //
    // Because in each loop, we delete molecules in that cluster from the vector
    // "unmappedMolecule", it is not wise to use for-loop. Instead we use a while-
    // loop.
    //
    while (!unmappedMolecule.empty()) {
        // We create a vector containing only one molecule and a pointer to it.
        vector<int> tempCluster = {unmappedMolecule[0]};

        unmappedMolecule.erase(unmappedMolecule.begin());
        mappedMolecule.insert(mappedMolecule.end(), tempCluster[0]);

        int pointer_processing = 0;

        // We now try to expand this cluster by neighborhood searching

        // We first perform neighborhood searching on "this" molecule.
        // There may be several atoms belonging to "this" molecule so here comes a for loop.

        bool completed_nb = false;

        while (!completed_nb) {
            bool foundNew = false;

            for (int atomIndex : idMap[tempCluster[pointer_processing]]) {
                AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(
                        sel.coordinates()[atomIndex]
                );
                AnalysisNeighborhoodPair pair;
                while (pairSearch.findNextPair(&pair)) {
                    if (find(tempCluster.begin(), tempCluster.end(), sel.mappedIds()[pair.refIndex()])
                        != tempCluster.end()) {

                        foundNew = true;

                        tempCluster.insert(tempCluster.end(), sel.mappedIds()[pair.refIndex()]);

                        // Now as we have assigned the cluster of refIndex, we delete it from
                        // unmapped molecule and insert it into mapped molecule

                        unmappedMolecule.erase(find(
                                tempCluster.begin(), tempCluster.end(), sel.mappedIds()[pair.refIndex()]));

                        mappedMolecule.insert(mappedMolecule.end(), sel.mappedIds()[pair.refIndex()]);
                    }
                }
            }

            pointer_processing += 1;

            if (!foundNew and pointer_processing == tempCluster.size()) {
                completed_nb = true;
            }
        }

        // Now tempCluster contains a new cluster.
        // We are now going to determine whether it is valid and put it into cluster list
        // if it is valid.

        if (tempCluster.size() >= cutoffClusterSize_) {
            clusterList.insert(clusterList.end(), tempCluster);
        }


    }

    int max_cluster_size = 0;

    for (vector<int> cluster : clusterList) {
        max_cluster_size = static_cast<int>(cluster.size() > max_cluster_size ?
                                            cluster.size() : max_cluster_size);
    }

    dhClusterCount.setPoint(0, clusterList.size());
    dhClusterSize.setPoint(0, max_cluster_size);

}

void assembly::finishAnalysis(int /*nframes*/) {}

void assembly::writeOutput() {}

