In this folder:
```shell
cc-10-Energy.txt is the pairwise energy derived based on observed frequency over expected frequency.
```

```shell
preprocess.py converts a set of PDB files to local coordinates centered around each amino acid, and save the data to coordinate.txt
```

```shell
ExtractEnergyDistribute.py compiles Nepre statistical potentials based on the data in coordinate.txt file.
```

There is an example script in the example folder. Note that a sufficient number of structures that with low redundancy are required to compile a good potential energy.
