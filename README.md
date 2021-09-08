# Citation

Please cite the following literature:

Tang et al., Prediction and characterization of liquid-liquid phase separation of minimalistic peptides, Cell Reports Physical Science (2021), https://doi.org/10.1016/j.xcrp.2021.100579


# Update 0.1: Fraction of aggregation

With "-molnumber" the program prints number of molecules that 
are included in clusters.

# Update 0.2: Liquidity

This program now calculates the liquidity, defined by the 
fraction of aggregated molecules that are still aggregated
in the next frame. Use '-liquidity' flag!

# Update 0.3: Density

This program now calculates the density of aggregated and solute state. Use "-density" and "-cutoff_multi" (optional) flags.


# Initial Release

# Introduction

Assembly Analyzer by Yiming Tang @ Fudan.

For peptide-self/co assembly studies, the system we meet usually 
contains hundreds of molecules. After they self/co-assembly into 
clusters it is very difficult to treat the periodic boundary 
conditions and cluster them into one or several big clusters. 
This program is thus constructed to serve the purpose of clustering 
and analyzing these systems.

# Input

This program always read in a trajectory and a topology (tpr) file, 
which are both crucial to the analysis. 

You CANNOT use a .gro or .pdb file instead which has no periodic 
boundary condition information.

An interaction is defined as two molecules with minimum distance 
larger than a cut-off distance defined by -cutoff_space 
(default: 0.4 nm). A cluster is defined as an aggregate larger 
than a certain number of molecules forming interactions, whereas 
the minimum number of molecules is defined by -cutoff_cz 
(default: 20). 

# Standard Output

The program calculates the number of big clusters and the size of 
the largest cluster as functions of time. File name can be specified
by -nb and -sz

# PDB Output

A PDB file containing the largest cluster at a certain time point 
(defined by -pdb_time is also written if a PDB file name is 
provided in \"-pdb\" flag. 

Periodic boundary condition of this cluster is well-treated to 
avoid across-boundary instances and the aggregate is put at 
the center of simulation box.

