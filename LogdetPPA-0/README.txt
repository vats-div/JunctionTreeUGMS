
[Introduction]
LogdetPPA version 0: 
A MATLAB software for log-determinant semidefinite programming
Copyright (c) 2009 by
Chengjing Wang, Defeng Sun, and Kim-Chuan Toh
----------------------------------------------------------------------------
This is software package for solving standard log-determinant
primal and dual SDP of the form: 

%% (P)  min  (C1,X1) + ... + (CN,XN) - mu1*logdet(X1) - ... - muN*logdet(XN)
%%      s.t. A1(X1) + ... + AN(XN) = b,
%%           X1,...,XN positive definite
%%
%% (D)  max  (b,y) + mu1*logdet(Z1) + ... + muN*logdet(ZN) + constant
%%      s.t. A1t(y) + Z1 = C1, ... ,ANt(y) + ZN = CN,
%%           Z1,..,ZN positive definite
%%
%%      where Xk, Zk are either symmetric positive definite matrices or positive vectors.
%%      The parameters mu1,...,muN are given positive scalars.
%%      The constant term in (D) depends on the demensions of X1,...,XN and mu1,...,muN.
----------------------------------------------------------------------------
The main algorithm is a semi-smooth Newton CG primal proximal 
point algorithm applied to (P); details can be found in the following reference: 

[Reference]
Chengjing Wang, Defeng Sun, and Kim-Chuan Toh, 
Solving log-determinant optimization problems by a Newton-CG primal 
proximal point algorithm,  preprint, National University of Singapore, 
September 2009.

[Copyright] 
See Copyright.txt
----------------------------------------------------------------------------
[Installation] 
The software LogdetPPA requires a few mex-files that you may need
to compile in MATLAB. 
To so, run MATALB in the directory LogdetPPA, then type: 

>> Installmex 

After that, to see whether you have installed the software
correctly, type the following in MATLAB: 

>> logdetPPAdemo
----------------------------------------------------------------------------