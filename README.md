# NEPRE--CodePackage

![](https://img.shields.io/badge/license-GNU-seagreen.svg?style=flat-square)
![](https://img.shields.io/badge/version-V2.0-blue.svg?style=flat-square)
![](https://img.shields.io/badge/language-Python-seagreen.svg?style=flat-square)
![](https://img.shields.io/badge/platform-Linux|Windows-blue.svg?style=flat-square&logo=Linux)

[LiuLab](http://liulab.csrc.ac.cn) | [Article](https://www.biorxiv.org/content/10.1101/463554v1) | [WebServer](http://liulab.csrc.ac.cn:10004/index/) | [Code Package](https://github.com/LiuLab-CSRC/Nepre) 

Introduction
---------------
Nepre is a scoring function for calculating the protein's potential energy. And it will help to predict the 3D structure of a protein.
Nepre is designed to be **efficient**, **flexible** and **protable** with **two** algorithms included:
* Calculate the potential energy using **fixed** cutoff (Nepre-F).
* Calculate the potential energy using the cutoff defined by the radius of amino acid (Nepre-R). 

#### Nepre-F
Users can use the cutoff between 4 angstrom and 10 angstrom to calculate the potential energy.
#### Nepre-R  
30,000 high-resolution protein are used to get statistical result of radius and fitted by the gaussian distribution function. Then Nepre will calculate the potential energy according to the gaussian mean data of radius of each kind of amino acid.

CopyRight
-------------
**Nepre** is created by **LiuLab** of **Beijing Computation Science Research Center(CSRC)**.

Usage
-------------
Nepre is built completely on **python2.7**. Some basic
modules are required for using Nepre.  
* numpy
* tqmd

Contact Us
-------------
**Email**: nepre2018@163.com  
**Address**: Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China.
