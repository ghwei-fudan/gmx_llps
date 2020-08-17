//
// Created by Yiming Tang on 2018/10/6.
//

#include "assembly.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <gromacs/fileio/pdbio.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/pbcutil/rmpbc.h>
#include <gromacs/selection/nbsearch.h>



using namespace std;
using namespace gmx;

bool compareClusterSize(vector<int> &cluster1, vector<int> &cluster2)
{
    return cluster1.size() > cluster2.size();
}

assembly::assembly() : cutoffSpace_(0.40)
{
    registerAnalysisDataset(&dataClusterCount_, "clusterCount");
    registerAnalysisDataset(&dataLargestCluster_, "clusterMax");
    registerAnalysisDataset(&dataMoleculesInCluster_, "moleculesInCluster");
    registerAnalysisDataset(&dataLiquidity_, "liquidity");
    registerAnalysisDataset(&dataDensity_, "density");
}

void assembly::initOptions(IOptionsContainer          *options,
                           TrajectoryAnalysisSettings *settings) {
    static const char *const desc[] =
            {
                    "Assembly Analyzer by Yiming Tang @ Fudan.\n",
                    "For peptide-self/co assembly studies, the system we meet usually contains ",
                    "hundreds of molecules. After they self/co-assembly into clusters it is ",
                    "very difficult to treat the periodic boundary conditions and ",
                    "cluster them into one or several big clusters. This program is thus constructred",
                    "to serve the purpose of clustering and analyzing these systems.\n",
                    "This program always read in a trajectory and a topology (tpr) file, which ",
                    "are both crucial to the analysis. You CANNOT use a .gro or .pdb file instead ",
                    "which has no periodic boundary condition information.\n",
                    "An interaction is defined as two molecules with minimum distance larger than a ",
                    "cut-off distance defined by -cutoff_space (default: 0.4 nm). A cluster is defined ",
                    "as an aggregate larger than a certain number of molecules forming interactions, whereas ",
                    "the minimum number of molecules is defined by -cutoff_cz (default: 20). ",
                    "The program calculates the number of big clusters and the size of ",
                    "the largest cluster as functions of time. File name can be specified",
                    "by -nb and -sz\n",
                    "A PDB file containing the largest cluster at a certain time point (defined by ",
                    "-pdb_time is also written if a PDB file name is provided in \"-pdb\" flag. ",
                    "Periodic boundary condition of this cluster is well-treated to avoid across-boundary ",
                    "instances and the aggregate is put at the center of simulation box.",
                    "This program will be further modified without notice at need of my own academic researches. ",
                    "If you have further question, feel free to contact tym@tymworld.com"
            };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("nb")
                               .filetype(eftPlot).outputFile()/*.required()*/
                               .store(&fnClusterCount_).defaultBasename("number")
                               .description("Number of clusters as a function of time"));

    options->addOption(FileNameOption("sz")
                               .filetype(eftPlot).outputFile()/*required()*/
                               .store(&fnLargestCluster_).defaultBasename("size")
                               .description("Size of the largest cluster as a function of time"));

    options->addOption(FileNameOption("molnumber")
                               .filetype(eftPlot).outputFile()
                               .store(&fnMoleculesInCluster_).defaultBasename("molecules_in_clusters")
                               .description("Number of molecules that are included in clusters"));

    options->addOption(FileNameOption("pdb")
                               .filetype(eftPDB).outputFile()
                               .store(&fnLargestClusterPDB_)
                               .description("The Largest Cluster of the Last Frame"));

    options->addOption(FileNameOption("liquidity")
                               .filetype(eftPlot).outputFile()
                               .store(&fnLiquidity_)
                               .description("Fraction of aggregated molecules in last frame that remain aggregated in this frame."));

    options->addOption(FileNameOption("density")
                               .filetype(eftPlot).outputFile()
                               .store(&fnDensity_)
                               .description("Density of aggregated / dispersed phase"));

    options->addOption(DoubleOption("cutoff_multi")
                               .store(&density_nb_multiplier_).defaultValue(2).required()
                               .description("Radius of Density Probe devided by cutoff_space"));

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

    options->addOption(BooleanOption("all_clusters").store(&write_all_clusters_)
                               .defaultValue(FALSE)
                               .description(
                                       "If TRUE, all cluster will be written to pdb file. If FALSE, only the largest will. Must use with -pdb"
                                       ));


    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void assembly::initAnalysis(const TrajectoryAnalysisSettings &settings,
                            const TopologyInformation &top) {

    //cout << "miao"<< endl;

    nb_.setCutoff(cutoffSpace_);
    nbDensity_.setCutoff(cutoffSpace_ * density_nb_multiplier_);
    dataLiquidity_.setColumnCount(0, 5);
    //cout << dataLiquidity_.dataSetCount();
    dataLargestCluster_.setColumnCount(0, 1);
    //cout << dataLiquidity_.columnCount(0);
    dataClusterCount_.setColumnCount(0, 1);
    dataMoleculesInCluster_.setColumnCount(0, 1);
    dataDensity_.setColumnCount(0, 2);


    if(!fnLiquidity_.empty()){
        AnalysisDataPlotModulePointer plotLiquidity(
                new AnalysisDataPlotModule(settings.plotSettings())
                );
        this->time_last_frame_ = -1;
        plotLiquidity->setFileName(fnLiquidity_);
        plotLiquidity->setTitle("Fraction of preservation");
        plotLiquidity->setXAxisIsTime();
        plotLiquidity->setYLabel("Fraction");
        plotLiquidity->appendLegend("Fraction of preservation");
        plotLiquidity->appendLegend("Fraction of cluster growth");
        plotLiquidity->appendLegend("Fraction of cluster shrink");
        plotLiquidity->appendLegend("Fraction of aggregation");
        plotLiquidity->appendLegend("Preservation of contacts");
        plotLiquidity->setYFormat(1, 7);
        dataLiquidity_.addModule(plotLiquidity);
    }

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

    if(!fnMoleculesInCluster_.empty()) {
        AnalysisDataPlotModulePointer plotMoleculesInCluster(
                new AnalysisDataPlotModule(settings.plotSettings())
        );
        plotMoleculesInCluster->setFileName(fnMoleculesInCluster_);
        plotMoleculesInCluster->setTitle("Number of Molecules in Clusters");
        plotMoleculesInCluster->setXAxisIsTime();
        plotMoleculesInCluster->setYLabel("Number of Molecules");
        dataMoleculesInCluster_.addModule(plotMoleculesInCluster);
    }

    if(!fnDensity_.empty())
    {
        AnalysisDataPlotModulePointer plotDensity(
                new AnalysisDataPlotModule(settings.plotSettings())
                );
        plotDensity->setFileName(fnDensity_);
        plotDensity->setTitle("Density");
        plotDensity->setXAxisIsTime();
        plotDensity->setYLabel("Density");
        dataDensity_.addModule(plotDensity);
    }

    this->top_ = top.topology();
    this->atoms_ = top.topology()->atoms;

    sel_.initOriginalIdsToGroup(top.mtop(), INDEX_MOL);

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

    //cout << "DEBUG: Column count for dataset 0 in liquidity is:" << dataLiquidity_.columnCount(0) << endl;
    //cout << "We've finished initializing." << endl;


}

void assembly::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                            TrajectoryAnalysisModuleData *pdata) {

    //cout << "We are analyzing Frames" << endl;

    AnalysisDataHandle dhClusterCount = pdata->dataHandle(dataClusterCount_);

    //cout << "DEBUG: We've just initialed datahandle for cluster count" << endl;

    AnalysisDataHandle dhClusterSize = pdata->dataHandle(dataLargestCluster_);

    //cout << "DEBUG: We've just initialed datahandle for largest cluster" << endl;

    AnalysisDataHandle dhLiquidity = pdata->dataHandle(dataLiquidity_);

    //cout << "DEBUG: We've just initialed datahandle for liquidity" << endl;

    AnalysisDataHandle dhMoleculesInCluster = pdata->dataHandle(dataMoleculesInCluster_);

    //cout << "DEBUG: We've just initialed datahandle for molecules in cluster" << endl;

    AnalysisDataHandle dhDensity = pdata->dataHandle(dataDensity_);

    //cout << "FRAME" << endl;



    const Selection &sel = sel_;

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel_);

    //cout << "DEBUG: We've just initialed neighborhood searching" << endl;

    dhClusterCount.startFrame(frnr, fr.time);

    //cout << "DEBUG: We've just started frame for dhClusterCount" << endl;

    dhClusterSize.startFrame(frnr, fr.time);

    //cout << "DEBUG: We've just started frame for dhClusterSize" << endl;

    dhMoleculesInCluster.startFrame(frnr, fr.time);

    dhLiquidity.startFrame(frnr, fr.time);

    dhDensity.startFrame(frnr, fr.time);

    //cout << "DEBUG: We've just started frame for dhMoleculesInCluster" << endl;

    // We first create two vector containing molecules that have been / not have been mapped.
    vector<int> mappedMolecule;
    vector<int> unmappedMolecule;
    for (int i = 0; i < molCount_; i++) {
        unmappedMolecule.insert(unmappedMolecule.end(), i);
    }

    unsigned long max_cluster_size = 0;
    unsigned long molecules_in_clusters = 0;
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

    // We now create a set that contains MOL ID of aggregated molecules

    set<int> cluster_this_frame;
    set<pair<int,int>> interaction_pairs;
    set<pair<int,int>> interaction_pairs_within_one_cluster;

    while (!unmappedMolecule.empty()) {
        // We create a vector containing only one molecule and a pointer to it.
        vector<int> tempCluster = {unmappedMolecule[0]};

        interaction_pairs_within_one_cluster.clear();

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

                std::pair<int,int> temp_pair;



                while (pairSearch.findNextPair(&pair)) {

                    if(sel.mappedIds()[pair.refIndex()] < sel.mappedIds()[atomIndex])
                    {
                        temp_pair.first = sel.mappedIds()[pair.refIndex()];
                        temp_pair.second = sel.mappedIds()[atomIndex];
                    }
                    else
                    {
                        temp_pair.first = sel.mappedIds()[atomIndex];
                        temp_pair.second = sel.mappedIds()[pair.refIndex()];
                    }

                    interaction_pairs_within_one_cluster.insert(temp_pair);

                    if (find(tempCluster.begin(), tempCluster.end(), sel.mappedIds()[pair.refIndex()])
                        == tempCluster.end()) {


                        foundNew = true;

                        tempCluster.insert(tempCluster.end(), sel.mappedIds()[pair.refIndex()]);
                        cluster_this_frame.insert(sel.mappedIds()[pair.refIndex()]);

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
            molecules_in_clusters += tempCluster.size();
            interaction_pairs.insert(interaction_pairs_within_one_cluster.begin(), interaction_pairs_within_one_cluster.end());
        }


    }



    //cout << "hhh" << endl;


    dhClusterCount.setPoint(0, clusterList.size());
    dhClusterSize.setPoint(0, max_cluster_size);
    dhMoleculesInCluster.setPoint(0, molecules_in_clusters);

    // Now it is time to calculate liquidity

    if(!fnLiquidity_.empty()){

        if(this->time_last_frame_ == -1 )
        {
            dhLiquidity.setPoint(0, 0);
            dhLiquidity.setPoint(1, float(cluster_this_frame.size()) / this->molCount_);
            dhLiquidity.setPoint(2, 0);
            dhLiquidity.setPoint(3, float(cluster_this_frame.size()) / this->molCount_);
            dhLiquidity.setPoint(4, 0);
        }
        else
        {
            if(this->cluster_last_frame_.size() == 0)
            {
                dhLiquidity.setPoint(0, 0);
                dhLiquidity.setPoint(1, float(cluster_this_frame.size()) / this->molCount_);
                dhLiquidity.setPoint(2, 0);
                dhLiquidity.setPoint(3, float(cluster_this_frame.size()) / this->molCount_);
            }
            if(this->interaction_pairs_last_frame_.size() == 0)
            {
                dhLiquidity.setPoint(4, 0);
            }
            else
            {
                //set<int> intersection;
                //set<int>::iterator p;

                std::vector<int> cluster_intersection(this->molCount_);

                std::vector<int>::iterator it;

                it = std::set_intersection(this->cluster_last_frame_.begin(), this->cluster_last_frame_.end(),
                                           cluster_this_frame.begin(), cluster_this_frame.end(), cluster_intersection.begin());

                cluster_intersection.resize(it - cluster_intersection.begin());

                dhLiquidity.setPoint(0, float(cluster_intersection.size()) / float(cluster_last_frame_.size()));
                dhLiquidity.setPoint(1, float(cluster_this_frame.size() - cluster_intersection.size()) / this->molCount_ );
                dhLiquidity.setPoint(2, float(cluster_last_frame_.size() - cluster_intersection.size()) / this->molCount_);
                dhLiquidity.setPoint(3, float(cluster_this_frame.size()) / this->molCount_);

                std::vector<std::pair<int,int>> intersection_interactions_pairs((int)max(interaction_pairs.size(), interaction_pairs_last_frame_.size()));
                std::vector<std::pair<int,int>>::iterator it_interactions_pairs;

                it_interactions_pairs = std::set_intersection(this->interaction_pairs_last_frame_.begin(),
                        this->interaction_pairs_last_frame_.end(), interaction_pairs.begin(), interaction_pairs.end(),
                        intersection_interactions_pairs.begin());

                intersection_interactions_pairs.resize(it_interactions_pairs - intersection_interactions_pairs.begin());

                dhLiquidity.setPoint(4, float(intersection_interactions_pairs.size()) / interaction_pairs_last_frame_.size());


            }
        }



        this->time_last_frame_ = fr.time;
        this->cluster_last_frame_ = cluster_this_frame;
        this->interaction_pairs_last_frame_ = interaction_pairs;

    }

    // Now let's calculate density. This usually takes time.

    if(!fnDensity_.empty()) {
        /// First we get list of aggregated atoms. These will serve as reference points for building nbsearch mesh.
        /// We also calculate mass here.

        long double aggregated_mass = 0;
        long double dispersed_mass = 0;

        std::vector<RVec> aggregated_atom;

        gmx::ArrayRef<float const[3]>::iterator iter_coordinate;
        gmx::ArrayRef<const int>::iterator iter_index;
        gmx::ArrayRef<const float>::iterator iter_mass;

        for (iter_coordinate = sel.coordinates().begin(), iter_index = sel.atomIndices().begin(), iter_mass = sel.masses().begin();
             iter_coordinate != sel.coordinates().end() && iter_index != sel.atomIndices().end() &&
             iter_mass != sel.masses().end();
             ++iter_coordinate, ++iter_index, ++iter_mass) {
            if (find(cluster_this_frame.begin(), cluster_this_frame.end(),
                     sel.mappedIds()[*iter_index]) != cluster_this_frame.end()) {
                aggregated_atom.insert(aggregated_atom.end(), *iter_coordinate);
                aggregated_mass += *iter_mass;
            } else {
                dispersed_mass += *iter_mass;
            }
        }

        /// We now initialize the nbsearch algorithm

        const std::vector<RVec> const_aggregated_atom = aggregated_atom;
        const AnalysisNeighborhoodPositions aggregated_positions(const_aggregated_atom);
        AnalysisNeighborhoodSearch nbDensitySearch = nbDensity_.initSearch(pbc, aggregated_positions);

        // We now start searching to get volume informations.

        double aggregated_point_number = 0;
        double dispersed_point_number = 0;

        //cerr << "DEBUG: box vector is " << pbc->box[0][0] << ", " << pbc->box[1][1] << ", " << pbc->box[2][2] << endl;

        for (real ix = 0; ix < pbc->box[0][0]; ix += 0.1) {
            for (real iy = 0; iy < pbc->box[1][1]; iy += 0.1) {
                for (real iz = 0; iz < pbc->box[2][2]; iz += 0.1) {
                    const rvec test_position = {ix, iy, iz};
                    AnalysisNeighborhoodPairSearch pairDensitySearch = nbDensitySearch.startPairSearch(test_position);
                    AnalysisNeighborhoodPair pair;
                    if (pairDensitySearch.findNextPair(&pair)) {
                        aggregated_point_number += 1;
                    } else {
                        dispersed_point_number += 1;
                    }
                }
            }
        }

        //cout << "Point Number: " << aggregated_point_number << ", " << dispersed_point_number << endl;

        double aggregated_volume =
                aggregated_point_number / 1000;
        double dispersed_volume =
                dispersed_point_number / 1000;

        // Now we calculate density

        //cout << "Mass:    " << aggregated_mass << ", " << dispersed_mass << endl;
        //cout << "Volume:  " << aggregated_volume << ", " << dispersed_volume << endl;
        //cout << "Density: " << aggregated_mass / aggregated_volume / (0.602) << ", " << dispersed_mass / dispersed_volume / (0.602) << endl;

        dhDensity.setPoint(0, (real)(aggregated_mass / aggregated_volume / 0.602));
        dhDensity.setPoint(1, (real)(dispersed_mass / dispersed_volume / 0.602));

    }








    // Now it is time to write pdb

    if (!fnLargestClusterPDB_.empty() && fr.time == time_pdb_) {

        // The clusterList variable contains all clusters.
        // Unfortunately they are unsorted.
        // We first sort all clusters in clusterList

        //cout << "We'll now write the pdb file" << endl;
        sort(clusterList.begin(), clusterList.end(), compareClusterSize);

        // We now open a pdb file

        FILE *fpPDB = fopen(fnLargestClusterPDB_.c_str(), "w");

        //cout << "We've open the PDB file" << endl;

        gmx_write_pdb_box(fpPDB, pbc->ePBC, pbc->box);

        // We now process all clusters one by one

        bool writen_largest_cluster = FALSE;

        for(auto tempCluster : clusterList) {

            if (not write_all_clusters_ and writen_largest_cluster) {
                break;
            }


            ////////////////////// REMOVING PERIODIC BOUNDARY CONDITIONS /////////////////////////

            // We first determine the number of atoms in this cluster.

            int atomNumber = 0;
            for (int moleculeIndex : tempCluster) {
                atomNumber += idMap[moleculeIndex].size();
            }

            cout << "Number of all atoms in this cluster are: " << atomNumber << endl;

            // We now initialize the atom coordinate list

            rvec *tempClusterAtoms = new rvec[atomNumber];

            atomNumber = 0;

            for (int moleculeIndex : tempCluster) {
                for (int atomIndex : idMap[moleculeIndex]) {
                    tempClusterAtoms[atomNumber][0] = sel.coordinates()[atomIndex][0];
                    tempClusterAtoms[atomNumber][1] = sel.coordinates()[atomIndex][1];
                    tempClusterAtoms[atomNumber][2] = sel.coordinates()[atomIndex][2];
                    atomNumber += 1;
                }
            }

            const rvec *clusterAtoms = tempClusterAtoms;

            // cout << "Number of all atoms are: " << atomNumber << endl;

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

            for (int i = 1; i < atomNumber; i++) {
                unmappedAtoms.insert(unmappedAtoms.end(), i);
            }

            mapped_unprocessedAtoms.insert(mapped_unprocessedAtoms.end(), 0);


            while (!unmappedAtoms.empty()) {
                int processingIndex = mapped_unprocessedAtoms[0];

                // cout << "DEBUG: Now Processing: " << processingIndex << endl;

                AnalysisNeighborhoodPairSearch pairSearchCluster = nbSearchCluster.startPairSearch(
                        clusterAtoms[processingIndex]
                );

                AnalysisNeighborhoodPair pairCluster;
                while (pairSearchCluster.findNextPair(&pairCluster)) {

                    if (find(mapped_processedAtoms.begin(), mapped_processedAtoms.end(), pairCluster.refIndex())\
 == mapped_processedAtoms.end()) {

                        // If in this condition, then the atom is not fully processed, which needs further processing.

                        if (find(mapped_unprocessedAtoms.begin(), mapped_unprocessedAtoms.end(), pairCluster.refIndex())\
 == mapped_unprocessedAtoms.end()) {
                            // If in this condition, then the atom has not been mapped. Now we map it.
                            mapped_unprocessedAtoms.insert(mapped_unprocessedAtoms.end(), pairCluster.refIndex());
                            unmappedAtoms.erase(
                                    find(unmappedAtoms.begin(), unmappedAtoms.end(), pairCluster.refIndex()));

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

            for (int i = 0; i < atomNumber; i++) {
                centroid_x += outputClusterAtoms[i][0];
                centroid_y += outputClusterAtoms[i][1];
                centroid_z += outputClusterAtoms[i][2];
            }

            centroid_x /= atomNumber;
            centroid_y /= atomNumber;
            centroid_z /= atomNumber;

            for (int i = 0; i < atomNumber; i++) {
                outputClusterAtoms[i][0] -= (centroid_x - pbc->box[0][0] / 2);
                outputClusterAtoms[i][1] -= (centroid_y - pbc->box[1][1] / 2);
                outputClusterAtoms[i][2] -= (centroid_z - pbc->box[2][2] / 2);


            }


            /////////////////////////////////////////////////////////////////////////////////////



            int fragmentCount = 1;
            int atomCount = 1;

            // cout << "DEBUG: Residue Count: " << atoms_.resinfo->rtp << endl;

            for (int moleculeIndex : tempCluster) {
                for (int atomIndex : idMap[moleculeIndex]) {


                    gmx_fprintf_pdb_atomline(fpPDB, epdbATOM, atomCount, *(atoms_.atomname[atomIndex]), \
                    /*atoms_.pdbinfo->altloc*/ ' ', *(top_->atoms.resinfo[atoms_.atom[atomIndex].resind].name), \
                    ' '/*top_->atoms.resinfo[atoms_.atom[atomIndex].resind].chainid*/ , \
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
            fprintf(fpPDB, "ENDMDL\n");
            writen_largest_cluster = TRUE;

        }

        fclose(fpPDB);
    }



    dhClusterCount.finishFrame();
    dhClusterSize.finishFrame();
    dhMoleculesInCluster.finishFrame();
    dhLiquidity.finishFrame();
    dhDensity.finishFrame();

}

void assembly::finishAnalysis(int /*nframes*/) {}

void assembly::writeOutput() {}

