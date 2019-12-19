# Nepre-R
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* Nepre_R.py (Main program)
* AminoAcid.py (Class for establish amino acid)
* radius.npy (Energy matrix)

A statistic of radius of different kind of amino acid have been done by us. We use gaussian distribute function to 
fit them and use the **gaussian mean** data as the default value. See details in mean_radius.txt.

You can see help information by typing:
<pre><code>
Nepre@liulab:~$ python Nepre_R.py -h
usage: Nepre_R.py [-h] [-s | -m] [-o] pdb_filename

Nepre-R Scoring Function Created by CSRC

positional arguments:
  path_filename        Path to PDB file 

optional arguments:
  -h, --help    show this help message and exit
  -s, --single  calculate single PDB
  -m, --multi   calculate a series of PDB
  -o, --output  save the results as a text file in running folder
</code></pre>

For **single** protein potential energy calculate, go to linux shell and type:
<pre><code>
#Not save the results in a text file
Nepre@liulab:~$ python Nepre_R.py -s ../example/test_decoy_set/native.pdb

#Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre_R.py -s -o ../example/test_decoy_set/native.pdb
</code></pre>

For **multi-object** job you can type:
<pre><code>
#Not save the results in a text file
Nepre@liulab:~$ python Nepre_R.py -m ../example/test_decoy_set/

#Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre_R.py -m -o ../example/test_decoy_set/
</code></pre>

You can also using **Nepre_R.py** as a module:

<pre></code>
import Nepre_R

#select a protein
path = "./example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre_R.load_EnergyMatrix()

#calculate Nepre potential energy
E = Nepre_R.calculate_Energy(f,Matrix)
</code></pre>

Extensions
----------
Nepre module also provide some useful function:
* Calculate the pearson coefficient correlation.
* Extract data from standard PDB file.
<pre><code>
"""
Pearson Coefficient
"""
import Nepre_R
x = [1,2,3,4]
y = [1,2,3,4]
p = Nepre_R.Pearson(x,y)

"""
Extract Data
"""
import Nepre_R
f = open("./example.pdb")
res = []
for line in f.readlines():
    res.append(Nepre_R.extract_Data(line))
</code></pre>
