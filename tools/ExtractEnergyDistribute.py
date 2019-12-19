import os
import math
import numpy as np
import sys

ccDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}


cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}


if __name__ == "__main__":
    
    args = sys.argv[1:]
    if(len(args) != 5):
      print "\n usage: python ExtractEnergyDistribute.py <local_coord.txt> <cutoff> <N_theta> <N_phi> <pairing_frequency_file>"
      exit()

    coord_file = args[0]
    cutoff = int(args[1])
    N_theta = int(args[2])
    N_phi = int(args[3])
    neighboring_frequency_data = args[4]
   
    List = ccDict.keys()
    List.sort()
    fr = open(neighboring_frequency_data)
    i = 0
    for line in fr.readlines():
        line=line.strip().split()
        for j in range(20):
            ccDict[List[i]][List[j]]=float(line[j])
        i += 1
    fr.close()

    
    emptyarea = {}
   
    print("Start to read coordinate.txt")
    fr=open(coord_file, 'r')
    amino1=''
    amino2=''
    
    #line = fr.readline()
    for line in fr.readlines():
        line=line.rstrip('\n')
        if line=='':
            continue

        if len(line)==3:
            if amino2=='VAL' or amino1=='':
                amino1=line
                amino2=''
            else:
                amino2=line
                print amino1,amino2
        else:
            line=line.split()
            if amino2 not in cdDict[amino1]:
                cdDict[amino1][amino2]=[[float(line[0]),float(line[1]),float(line[2])]]
            else:
                cdDict[amino1][amino2].append([float(line[0]),float(line[1]),float(line[2])])
            
    fr.close()

    for amino1 in cdDict:
        for amino2 in cdDict[amino1]:
            emp = 0
            dualArray=np.zeros((N_theta,N_phi))
            print amino1,amino2
            for coor in cdDict[amino1][amino2]:
                coor = np.array(coor)
                rho = sum(coor**2)**0.5
                theta = np.arccos(coor[2]/rho)
                phi = np.arctan2(coor[1],coor[0])
                theta = min(int(theta*N_theta/np.pi),N_theta-1)
                phi = min(int((phi+np.pi)/(2*np.pi/N_phi)),N_phi-1)
                #rho = min(int(rho),10)
                
                if(rho < cutoff):
                    dualArray[theta][phi] += 1.0

            print("Start to find min data")
            mindata = None
            for theta in range(N_theta):
                for phi in range(N_phi):
                    if(mindata == None and dualArray[theta][phi] != 0):
                        mindata = dualArray[theta][phi]
                    if(mindata > dualArray[theta][phi] and dualArray[theta][phi] != 0):
                        mindata = dualArray[theta][phi]
            print("Start to check empty area")
            for theta in range(N_theta):
                for phi in range(N_phi):
                    if dualArray[theta][phi] == 0:
                        dualArray[theta][phi] = mindata
                        emp += 1
            dualArray = dualArray / dualArray.sum()

            Integral = np.ones((N_theta,N_phi))
            for j in range(N_theta):
                Integral[j] = Integral[j]*(np.cos(j*np.pi/N_theta)-np.cos((j+1)*np.pi/N_theta))

            Integral = Integral / Integral.sum()
            dualArray = dualArray / Integral

            cdDict[amino1][amino2] = ccDict[amino1][amino2]-np.log(dualArray)
            emptyarea[amino1+'-'+amino2] = emp

    f = file(str(cutoff) + ".npy",'wb')
    for amino1 in List:
        for amino2 in List:
            np.save(f, cdDict[amino1][amino2])

