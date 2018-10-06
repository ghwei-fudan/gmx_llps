//
// Created by Yiming Tang on 2018/10/6.
//

#include "assembly.h"
#include <iostream>
#include <string>
#include <vector>

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





}

void assembly::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                            TrajectoryAnalysisModuleData *pdata) {
    AnalysisDataHandle dhClusterCount = pdata->dataHandle(dataClusterCount_);
    AnalysisDataHandle dhClusterSize = pdata->dataHandle(dataLargestCluster_);
    const Selection &sel = sel_;

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel_);

    dhClusterCount.startFrame(frnr, fr.time);
    dhClusterSize.startFrame(frnr, fr.time);























}