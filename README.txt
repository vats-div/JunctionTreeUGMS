We provide Matlab code for the algorithms in the paper:

D. Vats and R. D. Nowak
A Junction Tree Framework for Undirected Graphical Model Selection

USAGE
-----
To use the scripts, you will first need to include all the directories in the Matlab path.
Note that many functions included are borrowed from other sources (mainly http://code.google.com/p/pmtk3/)

If you are on a Windows machine, run configurepathWin.

If you are on a MaC/Linux/Unix machine, run configurepathMac.

See ExampleScript.m for examples of using the various functions in the toolbox.

See the folder RealData/ for the results on real data given in the paper.

See the folder SyntheticExamples for scripts used to generate Tables 2-5 in the paper.  If you would like the data we used for the synthetic graphs, email at vats.div@gmail.com. 

UGMS ALGORITHMS
---------------

In the folder UGMSAlgorithms, we have provided implementations of three UGMS algorithms:

- graphical Lasso (UGMS_GLasso.m)
- neighborhood selection using adaptive Lasso (UGMS_NLasso.m)
- PC-Algorithm (UGMS_PC.m)

See the individual scripts or ExampleScript.m on how to use these functions.  

CONTACT
------

Divyanshu Vats
dvats@rice.edu or vats.div@gmail.com