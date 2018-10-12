//
// Created by Yiming Tang on 2018/10/6.
//

#include "assembly.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gromacs/fileio/pdbio.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/pbcutil/rmpbc.h>
#include <gromacs/selection/nbsearch.h>



using namespace std;
using namespace gmx;


assembly::assembly() : cutoffSpace_(0.40)
{
    registerAnalysisDataset(&dataClusterCount_, "clusterCount");
    registerAnalysisDataset(&dataLargestCluster_, "clusterMax");
}

void assembly::initOptions(IOptionsContainer          *options,
                           TrajectoryAnalysisSettings *settings) {
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

    options->addOption(FileNameOption("pdb")
                               .filetype(eftPDB).outputFile()
                               .store(&fnLargestClusterPDB_)
                               .description("The Largest Cluster of the Last Frame"));

    options->addOption(SelectionOption("select")
                               .store(&sel_).required()
                               .description("Group that contain your molecules"));

    options->addOption(DoubleOption("cutoff_cz").store(&cutoffClusterSize_)
                               .defaultValue(20).required()
                               .description("Minimum number of molecules in a so-called cluster"));

    options->addOption(DoubleOption("cutoff_space").store(&cutoffSpace_)
                               .defaultValue(0.4).required()
                               .description("Maximum atom-wise distance to insure an interaction"));

    options->addOption(DoubleOption("pdb_time").timeValue().store(&time_pdb_)
                               .description(
                                       "Largest cluster of this time will be written to pdb if \"pdb\" specified."));


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

    molCount_ = sel_.mappedIds()[sel_.atomCount() - 1] + 1;
    cout << "Total number of molecules in this system: " << molCount_ << endl;

    // Now we are going to generate "molecule ID -> atom ID" map

    // this->idMap = new vector<int>[molCount_];

    for (int i = 0; i < molCount_; i++)
    {
        vector<int> tempVector;
        idMap.insert(idMap.end(), tempVector);
    }

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

    //cout << "DEBUG: We've just initialed neighborhood searching" << endl;

    dhClusterCount.startFrame(frnr, fr.time);

    //cout << "DEBUG: We've just started frame for dhClusterCount" << endl;

    dhClusterSize.startFrame(frnr, fr.time);

    //cout << "DEBUG: We've just started frame for dhClusterSize" << endl;

    // We first create two vector containing molecules that have been / not have been mapped.
    vector<int> mappedMolecule;
    vector<int> unmappedMolecule;
    for (int i = 0; i < molCount_; i++) {
        unmappedMolecule.insert(unmappedMolecule.end(), i);
    }

    unsigned long max_cluster_size = 0;
    vector<int> maxCluster;

/*
    cout << "The unmapped map: "<< endl;
    for (int i=0; i<unmappedMolecule.size(); i++)
    {
        cout << unmappedMolecule[i] << ',' ;
    }
    cout << endl;
*/
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

        //cout << "DEBUG: Now Process Molecule: " << tempCluster[0] << endl;


        unmappedMolecule.erase(unmappedMolecule.begin());
        //cout << "DEBUG: The next will be:     " << unmappedMolecule[0] << endl;

        mappedMolecule.insert(mappedMolecule.end(), tempCluster[0]);

        int pointer_processing = 0;

        // We now try to expand this cluster by neighborhood searching

        // We first perform neighborhood searching on "this" molecule.
        // There may be several atoms belonging to "this" molecule so here comes a for loop.

        bool completed_nb = false;

        while (!completed_nb) {
            bool foundNew = false;

            //////////////////////////DEBUG/////////////////////////////
            //cout << "DEBUG: We are in loop" << endl;

            //cout << "DEBUG: The Molecule Number is :" << tempCluster[pointer_processing] << endl;

            //cout << "DEBUG: The Atoms are :" ;
            //for (int i : idMap[tempCluster[pointer_processing]]){
            //    cout << i << ' ';
            //}
            //cout << endl;
            //////////////////////////DEBUG/////////////////////////////


            for (int atomIndex : idMap[tempCluster[pointer_processing]]) {


                AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(
                        sel.coordinates()[atomIndex]
                );

                //////////////////////////DEBUG/////////////////////////////

                //cout << "DEBUG: Now doing neighborhood searching for atom: " << atomIndex;
                //cout << " with (x,y,z) of " << sel.coordinates()[atomIndex][0] << ", "  ;
                //cout << sel.coordinates()[atomIndex][1] << ", " ;
                //cout << sel.coordinates()[atomIndex][2] << endl;



                //////////////////////////DEBUG/////////////////////////////

                AnalysisNeighborhoodPair pair;
                while (pairSearch.findNextPair(&pair)) {
                    if (find(tempCluster.begin(), tempCluster.end(), sel.mappedIds()[pair.refIndex()])
                        == tempCluster.end()) {


                        foundNew = true;

                        tempCluster.insert(tempCluster.end(), sel.mappedIds()[pair.refIndex()]);

                        // Now as we have assigned the cluster of refIndex, we delete it from
                        // unmapped molecule and insert it into mapped molecule

                        //cout << "DEBUG: Now found " << *find(
                        //        tempCluster.begin(), tempCluster.end(), sel.mappedIds()[pair.refIndex()]) << endl;
/*
                        cout << "The unmapped map: "<< endl;
                        for (int i=0; i<unmappedMolecule.size(); i++)
                        {
                            cout << unmappedMolecule[i] << ',' ;
                        }
                        cout << endl;
*/
                        //cout << "We are going to erase from unmapped list: " << sel.mappedIds()[pair.refIndex()] << endl;



                        unmappedMolecule.erase(find(unmappedMolecule.begin(), unmappedMolecule.end(),
                                                    sel.mappedIds()[pair.refIndex()]));

                        /*
                        cout << "The unmapped map: "<< endl;
                        for (int i=0; i<unmappedMolecule.size(); i++)
                        {
                            cout << unmappedMolecule[i] << ',' ;
                        }
                        cout << endl << endl;
*/

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

        if (tempCluster.size() > max_cluster_size) {
            max_cluster_size = tempCluster.size();
            maxCluster = tempCluster;
        }

        if (tempCluster.size() >= cutoffClusterSize_) {
            clusterList.insert(clusterList.end(), tempCluster);
        }


    }



    //cout << "hhh" << endl;


    dhClusterCount.setPoint(0, clusterList.size());
    dhClusterSize.setPoint(0, max_cluster_size);

    // Now it is time to write pdb

    if (!fnLargestClusterPDB_.empty() && fr.time == time_pdb_) {


        ////////////////////// REMOVING PERIODIC BOUNDARY CONDITIONS /////////////////////////

        // We first determine the number of atoms in this cluster.

        int atomNumber = 0;
        for(int moleculeIndex : maxCluster){
            atomNumber += idMap[moleculeIndex].size();
        }

        cout << "Number of all atoms are: " << atomNumber << endl;

        // We now initialize the atom coordinate list

        rvec *tempClusterAtoms = new rvec[atomNumber];

        atomNumber = 0;

        for (int moleculeIndex : maxCluster) {
            for (int atomIndex : idMap[moleculeIndex]) {
                tempClusterAtoms[atomNumber][0] = sel.coordinates()[atomIndex][0];
                tempClusterAtoms[atomNumber][1] = sel.coordinates()[atomIndex][1];
                tempClusterAtoms[atomNumber][2] = sel.coordinates()[atomIndex][2];
                atomNumber += 1;
            }
        }

        const rvec *clusterAtoms = tempClusterAtoms;

        cout << "Number of all atoms are: " << atomNumber << endl;

        // We now initial a neighborhood searching based on the grid generated by the above coordinate list.

        AnalysisNeighborhoodPositions nPLargestCluster(clusterAtoms, atomNumber);
        AnalysisNeighborhood nbLargestCluster;
        nbLargestCluster.setCutoff(static_cast<real>(cutoffSpace_));

        AnalysisNeighborhoodSearch nbSearchCluster = nbLargestCluster.initSearch(pbc, nPLargestCluster);

        // We now do neighborhood searching one by one and get the outcome coordinates.

        rvec *outputClusterAtoms = new rvec[atomNumber];

        outputClusterAtoms[0][0] = clusterAtoms[0][0];
        outputClusterAtoms[0][1] = clusterAtoms[0][1];
        outputClusterAtoms[0][2] = clusterAtoms[0][2];

        vector<int> unmappedAtoms;
        vector<int> mapped_unprocessedAtoms;
        vector<int> mapped_processedAtoms;

        for(int i=1; i<atomNumber; i++){
            unmappedAtoms.insert(unmappedAtoms.end(), i);
        }

        mapped_unprocessedAtoms.insert(mapped_unprocessedAtoms.end(), 0);


        while(!unmappedAtoms.empty())
        {
            int processingIndex = mapped_unprocessedAtoms[0];

            // cout << "DEBUG: Now Processing: " << processingIndex << endl;

            AnalysisNeighborhoodPairSearch pairSearchCluster = nbSearchCluster.startPairSearch(
                    clusterAtoms[processingIndex]
                    );

            AnalysisNeighborhoodPair pairCluster;
            while (pairSearchCluster.findNextPair(&pairCluster)) {

                if(find(mapped_processedAtoms.begin(), mapped_processedAtoms.end(), pairCluster.refIndex())\
                == mapped_processedAtoms.end()){

                    // If in this condition, then the atom is not fully processed, which needs further processing.

                    if(find(mapped_unprocessedAtoms.begin(), mapped_unprocessedAtoms.end(),pairCluster.refIndex())\
                    == mapped_unprocessedAtoms.end()) {
                        // If in this condition, then the atom has not been mapped. Now we map it.
                        mapped_unprocessedAtoms.insert(mapped_unprocessedAtoms.end(), pairCluster.refIndex());
                        unmappedAtoms.erase(find(unmappedAtoms.begin(), unmappedAtoms.end(), pairCluster.refIndex()));

                        outputClusterAtoms[pairCluster.refIndex()][0] = outputClusterAtoms[processingIndex][0] +
                                                                        pairCluster.dx()[0];
                        outputClusterAtoms[pairCluster.refIndex()][1] = outputClusterAtoms[processingIndex][1] +
                                                                        pairCluster.dx()[1];
                        outputClusterAtoms[pairCluster.refIndex()][2] = outputClusterAtoms[processingIndex][2] +
                                                                        pairCluster.dx()[2];
                    }

                    mapped_processedAtoms.insert(mapped_processedAtoms.end(), pairCluster.refIndex());
                }

            }
            // Anyway now this atom is fully processed, so we move its position.
            mapped_processedAtoms.insert(mapped_processedAtoms.end(), processingIndex);
            mapped_unprocessedAtoms.erase(mapped_unprocessedAtoms.begin());

        }

        //////////////////////////////////////////////////////////////////////////////////////


        //////////////////////////////// CENTERING THE CLUSTER ///////////////////////////////

        real centroid_x = 0.0, centroid_y = 0.0, centroid_z = 0.0;

        for(int i=0; i<atomNumber; i++){
            centroid_x += outputClusterAtoms[i][0];
            centroid_y += outputClusterAtoms[i][1];
            centroid_z += outputClusterAtoms[i][2];
        }

        centroid_x /= atomNumber;
        centroid_y /= atomNumber;
        centroid_z /= atomNumber;

        for(int i=0; i<atomNumber; i++){
            outputClusterAtoms[i][0] -= centroid_x;
            outputClusterAtoms[i][1] -= centroid_y;
            outputClusterAtoms[i][2] -= centroid_z;
        }


        /////////////////////////////////////////////////////////////////////////////////////


        // cout << endl << "Now We are going to write a pdb file to: " << fnLargestClusterPDB_ << endl;

        FILE *fpPDB = fopen(fnLargestClusterPDB_.c_str(), "w");

        // cout << "DEBUG: We've opened the file." << endl;

        //gmx_write_pdb_box(fpPDB, pbc->ePBC, pbc->box);

        // cout << "DEBUG: We've written pbc information." << endl;

        int fragmentCount = 1;
        int atomCount = 1;

        // cout << "DEBUG: Residue Count: " << atoms_.resinfo->rtp << endl;

        for (int moleculeIndex : maxCluster) {
            for (int atomIndex : idMap[moleculeIndex]) {


                gmx_fprintf_pdb_atomline(fpPDB, epdbATOM, atomCount, *(atoms_.atomname[atomIndex]), \
                /*atoms_.pdbinfo->altloc*/ ' ', *(top_->atoms.resinfo[0].name), atoms_.resinfo->chainid, \
                atoms_.atom[atomIndex].resind,
                        /*atoms_.resinfo->nr,*/ (char) atoms_.resinfo->ic, \

                                         outputClusterAtoms[atomCount - 1][0] * 10,
                                         outputClusterAtoms[atomCount - 1][1] * 10,
                                         outputClusterAtoms[atomCount - 1][2] * 10,
                        /*
                        sel.coordinates()[atomIndex][0]*10, sel.coordinates()[atomIndex][1]*10, \
                        sel.coordinates()[atomIndex][2]*10, atoms_.pdbinfo->occup, atoms_.pdbinfo->bfac*/0, 0, \
                *(atoms_.atomname[atomIndex]));
                atomCount += 1;

                //atoms_.pdbinfo->

            }
            fragmentCount += 1;
        }

        fclose(fpPDB);
    }



    dhClusterCount.finishFrame();
    dhClusterSize.finishFrame();

}

void assembly::finishAnalysis(int /*nframes*/) {}

void assembly::writeOutput() {}

