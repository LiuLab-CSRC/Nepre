cutoff=6
pdb_folder=../example/test_decoy_set/
pdb_name='native.pdb'

echo 'Calculating energy for a single PDB structure with cutoff=6'

python Nepre_F.py -s $pdb_folder/$pdb_name $cutoff


echo 'Calculating energy for PDB structures in a folder with cutoff=6'
python Nepre_F.py -m $pdb_folder/ $cutoff
