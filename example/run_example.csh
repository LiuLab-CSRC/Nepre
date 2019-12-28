#!/bin/csh -f 
#
set pdb_folder = './for_stat_pdb/'  # the pdb folder name needs to end with '/'
set cutoff=10
set pdblist='pdblist.txt'

python ../tools/preprocess.py $pdb_folder $pdblist $cutoff

echo $pdb_folder $pdblist $cutoff

echo 'local coordinate mapping completed.\n data saved in coordinate.txt'

echo 'now compile the energy functions'

set n_theta=20
set n_phi=20
set pair_freq='../tools/cc-10-Energy.txt'
set coord_file = 'coordinate.txt'

python ../tools/ExtractEnergyDistribute.py $coord_file $cutoff $n_theta $n_phi $pair_freq

echo 'energy profile saved in ' $cutoff.npy
