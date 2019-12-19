import os
import sys
import numpy as np

aminoDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
           "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
           "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
           "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

aDict={"ALA":0,"VAL":0,"LEU":0,"ILE":0,"PHE":0,\
       "TRP":0,"MET":0,"PRO":0,"GLY":0,"SER":0,\
       "THR":0,"CYS":0,"TYR":0,"ASN":0,"GLN":0,\
       "HIS":0,"LYS":0,"ARG":0,"ASP":0,"GLU":0,}

ccDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

def EuclideanDistances(A, B):
    BT = B.transpose()
    vecProd = A * BT
    SqA =  A.getA()**2
    sumSqA = np.matrix(np.sum(SqA, axis=1))
    sumSqAEx = np.tile(sumSqA.transpose(), (1, vecProd.shape[1]))
    SqB = B.getA()**2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))
    SqED = sumSqBEx + sumSqAEx - 2*vecProd
    return np.matrix(SqED)

def processPDB2LocalCoord(fileNameList,cutoff, pdb_folder_path):
    for string in fileNameList:
        #string=string.lower()
        filename=pdb_folder_path + string
        try:
            fr=open(filename)
            #each DNA or RNA center
            center=[]
            coordinate=[]
            for line in fr.readlines():
                line = line.strip()
                if line=='ENDMDL':
                    break
                if  line[0:4]=='ATOM' and line[16]!='B' and \
                    line[17:20].split()[0] in aDict and line[12:16].split()[0]!='OXT': 
                    res=line[17:20].split()[0]
                    chain=line[21]
                    index=line[22:26].split()[0]
                    x=float(line[30:38].split()[0])
                    y=float(line[38:46].split()[0])
                    z=float(line[46:54].split()[0])
                    if center!=[] and res==center[-1][0] and \
                       chain==center[-1][1] and index==center[-1][2]:
                       coordinate[-1].append([x,y,z])
                       if res=='GLY' or line[12:16].split()[0] not in {'N','CA','C','O','O1','02'}:
                           center[-1][3]+=x
                           center[-1][4]+=y
                           center[-1][5]+=z
                           center[-1][6]+=1
                    else:
                        if len(center)>=1 and center[-1][-1]==0 or len(coordinate)>=1 and len(coordinate[-1])<4:
                            del center[-1],coordinate[-1]
                        center.append([res,chain,index])
                        center[-1].extend([0.0,0.0,0.0,0])
                        aDict[res]+=1
                        coordinate.append([[x,y,z]])
                        if res=='GLY' or line[12:16].split()[0] not in {'N','CA','C','O','O1','02'}:
                            center[-1][3]+=x
                            center[-1][4]+=y
                            center[-1][5]+=z
                            center[-1][6]+=1
            if len(center)>=1 and center[-1][-1]==0 or len(coordinate)>=1 and len(coordinate[-1])<4:
                del center[-1],coordinate[-1]
            #search ionic around amino acid
            if center!=[] and coordinate!=[]:
                vector=[]
                for coor in center:
                    vector.append([coor[3]/coor[6],coor[4]/coor[6],coor[5]/coor[6]])
                vector=np.matrix(vector)
                distanceMatrix=np.asarray(EuclideanDistances(vector, vector))
                vector=np.asarray(vector)

                xAxis=[]
                yAxis=[]
                for coor in coordinate:
                    xAxis.append(coor[0])
                    yAxis.append(coor[1])
                xAxis=np.array(xAxis)-vector
                yAxis=np.array(yAxis)-vector
                xAxis=xAxis/np.asarray(np.asmatrix(np.sqrt(np.sum(xAxis**2,axis=1))).transpose())
                yAxis=yAxis-xAxis*np.asarray(np.asmatrix(np.sum(np.multiply(xAxis,yAxis),axis=1)).transpose())
                yAxis=yAxis/np.asarray(np.asmatrix(np.sqrt(np.sum(yAxis**2,axis=1))).transpose())
                zAxis=np.cross(xAxis,yAxis)

            N=len(center)
            for i1 in range(N):
                for i2 in range(N):
                    if i1!=i2 and distanceMatrix[i1][i2]<cutoff and \
                       (center[i1][1]!=center[i2][1] or center[i1][2]!=center[i2][2]):
                        if center[i2][0] not in aminoDict[center[i1][0]]:
                            aminoDict[center[i1][0]][center[i2][0]]=1
                        else:
                            aminoDict[center[i1][0]][center[i2][0]]+=1
                        a=np.inner(vector[i2]-vector[i1],xAxis[i1])
                        b=np.inner(vector[i2]-vector[i1],yAxis[i1])
                        c=np.inner(vector[i2]-vector[i1],zAxis[i1])
                        if center[i2][0] not in cdDict[center[i1][0]]:
                            cdDict[center[i1][0]][center[i2][0]]=[[a,b,c]]
                        else:
                            cdDict[center[i1][0]][center[i2][0]].append([a,b,c])

        except IOError as err:
            continue
        finally:
            if 'fr' in locals():
                fr.close()
        print fileNameList.index(string), string, N, "residues processed"


if __name__ == "__main__":
    args = sys.argv[1:]
    if (len(args) != 3):
      print "usage: python preprocess.py <pdb_folder_path> <pdb_filename.txt> <distance_cutoff>"
      exit()

    pdb_folder_path = args[0]
    pdbNames = args[1]
    cutoff = float( args[2] )

    pdbNameList = [ ]
    f = open(pdbNames,'r')
    for ll in f.readlines():
      pdbNameList.append( ll.strip() )

    processPDB2LocalCoord(pdbNameList,cutoff**2, pdb_folder_path)


    for key in aminoDict:
        for key2 in aDict:
            if key2 not in aminoDict[key]:
                aminoDict[key][key2]=0
    for key in ccDict:
        for key2 in aminoDict:
            if key2 not in ccDict[key]:
                ccDict[key][key2]=0.0

    fw=open("coordinate.txt","w")
    a=cdDict.keys()
    a.sort()
    for key in a:
        fw.write(key)
        fw.write("\n")
        b=cdDict[key].keys()
        b.sort()
        for key2 in b:
            fw.write(key2)
            fw.write("\n")
            for List in cdDict[key][key2]:
                fw.write(str(round(List[0],3)))
                fw.write("   ")
                fw.write(str(round(List[1],3)))
                fw.write("   ")
                fw.write(str(round(List[2],3)))
                fw.write("\n")
            fw.write("\n")
    fw.close()

