# Nepre-F
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* Nepre_F.py (Main program)
* AminoAcid.py (Class for establish amino acid)
* cutoff.npy (Energy matrix)

We provide **7** cutoff options with cutoff between **4** angstrom and **10** angstrom.
You can see help information by typing:
<pre><code>
Nepre@liulab:~$ python Nepre_F.py -h
usage: Nepre_F.py [-h] [-s | -m] [-o] path cutoff

Nepre-F Scoring Function Created by CSRC

positional arguments:
  path          PDB file path of folder path
  cutoff        cutoff parameter for Nepre-F

optional arguments:
  -h, --help    show this help message and exit
  -s, --single  calculate single PDB
  -m, --multi   calculate a series of PDB
  -o, --output  save the results as a text file in running folder
</code></pre>

For **single** protein potential energy calculate, choose a cutoff (6 angstrom e.g) and type:

#Not save the results in a text file
```shell
Nepre@liulab:~$ python Nepre_F.py -s ../example/pdb/native.pdb 6

Nepre Potential Energy
Using Cutoff: 6
../example/pdb/native.pdb -74.65268164490764
```

#Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre_F.py -s -o ../example/pdb/native.pdb 6
</code></pre>

For **multi-object** calculation, you can type:
<pre><code>
#Not save the results in a text file
Nepre@liulab:~$ python Nepre_F.py -m ./pdb_folder/ 6

#Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre_F.py -m -o ./pdb_folder/ 6
</code></pre>

You can also use **Nepre_F.py** as a module if you want to use the calculation results for other purposes:
<pre><code>
import Nepre_F

#choose a cutoff
cutoff = 6

#select a protein
path = "./example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre_F.load_EnergyMatrix(cutoff)

#calculate Nepre potential energy
E = Nepre_F.calculate_Energy(f,Matrix,cutoff)
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
import Nepre_F
x = [1,2,3,4]
y = [1,2,3,4]
p = Nepre_F.Pearson(x,y)

"""
Extract Data
"""
import Nepre_F
f = open("./example.pdb")
res = []
for line in f.readlines():
    res.append(Nepre_F.extract_Data(line))
</code></pre>
