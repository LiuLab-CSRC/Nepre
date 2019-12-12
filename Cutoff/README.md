# Nepre-F
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* Nepre_F.py (Main program)
* AminoAcid.py (Class for establish amino acid)
* {cutoff}.npy (Energy matrix)

where cutoff is a number from 4 to 10, describing the distance threshold for the neighborhood. We provide **7** cutoff options with cutoff between **4** angstrom and **10** angstrom.

You can see help information by typing:
```Shell
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
```

For **single** protein potential energy calculate, choose a cutoff (6 angstrom e.g) and turn on **-s** flag:

## Print results to the terminal
```shell
Nepre@liulab:~$ python Nepre_F.py -s ../example/pdb/native.pdb 6
```
The calculation results are shown as:
```
Nepre Potential Energy
Using Cutoff: 6
../example/pdb/native.pdb -74.65268164490764
```

## Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
```Shell
Nepre@liulab:~$ python Nepre_F.py -s -o ../example/pdb/native.pdb 6
```

For **multi-object** calculation, you can use **-m** flag:

## Print results to the terminal
```Shell
Nepre@liulab:~$ python Nepre_F.py -m ../example/pdb/ 6
```
The results are:
```
Nepre Potential Energy
Using Cutoff: 6
decoy1_10.pdb 	-66.87601505181286
decoy1_38.pdb 	-69.27243871698501
decoy1_9.pdb 	-73.25992941774311
decoy1_12.pdb 	-72.71724455052468
decoy1_13.pdb 	-59.75656264882953
decoy1_8.pdb 	-67.00973379132571
decoy1_17.pdb 	-66.47342460888504
decoy1_16.pdb 	-72.80176813916913
decoy1_14.pdb 	-70.55491570484827
decoy1_15.pdb 	-73.71455055593965
decoy1_3.pdb 	-68.16973959578826
decoy1_18.pdb 	-67.44919368439325
decoy1_19.pdb 	-69.57480105848336
decoy1_27.pdb 	-68.63085067776933
decoy1_26.pdb 	-67.23543088201664
decoy1_1.pdb 	-73.25169055504813
decoy1_5.pdb 	-71.45508750375832
decoy1_36.pdb 	-77.04372881925707
decoy1_4.pdb 	-69.80589852136566
decoy1_6.pdb 	-69.44655719981684
decoy1_35.pdb 	-69.38003347950509
native.pdb 	-74.65268164490764
decoy1_43.pdb 	-68.45272155017172
```

## Save the results in a text file(Same folder with Nepre.py with name "latest_results.txt")
```Shell
Nepre@liulab:~$ python Nepre_F.py -m -o ../example/pdb/ 6
```

You can also use **Nepre_F.py** as a module if you want to use the calculation results for other purposes:
```Python
import Nepre_F

#choose a cutoff
cutoff = 6

#select a protein
path = "example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre_F.load_EnergyMatrix(cutoff)

#calculate Nepre potential energy
E = Nepre_F.calculate_Energy(f,Matrix,cutoff)
```

Extensions
----------
Nepre module also provide some useful function:
* Calculate the pearson coefficient correlation.
* Extract data from standard PDB file.
```Python
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
```
