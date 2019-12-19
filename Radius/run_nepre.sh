cutoff=6
pdb_folder=../example/test_decoy_set/
pdb_name='native.pdb'

echo 'Calculating energy for a single PDB structure with radius cutoff'

python Nepre_R.py -s $pdb_folder/$pdb_name 


echo 'Calculating energy for PDB structures in a folder with radius cutoff'
python Nepre_R.py -m $pdb_folder/ 
