%%******************************************************
%% Run this script in Matlab command window 
%% 
%%******************************************************

   function Installmex 

   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4)); 
%%
   if strcmp(computer_model,'PCWIN') 
      str1 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  ']; 
      str2 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
      libstr = [str1,str2];     
   elseif strcmp(computer_model,'PCWIN64')
      str1 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  ']; 
      str2 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
      libstr = [str1,str2];  
   else
      libstr = '-lmwlapack -lmwblas  '; 
   end
   mexcmd = 'mex -O  -largeArrayDims  -output ';    
%%
   if (matlabversion < 7.4) & (matlabrelease < 2009) 
      error(' LogdetPPA works only for MATLAB version 7.4 and above'); 
   end
   fsp = filesep;

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir); 
%%
%% generate mex files in mexfun
%%
   clear fname

   src = [curdir,fsp,'solver/mexfun']; 
   eval(['cd ','solver/mexfun']); 
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);       
   %%
   fname{1} = 'mexbwsolve';
   fname{2} = 'mexfwsolve';
   fname{3} = 'mexmat';
   fname{4} = 'mexsmat'; 
   fname{5} = 'mexsvec'; 
   fname{6} = 'mextriang'; 
   fname{7} = 'mexScalarxMatrix'; 
 
   ext = mexext;
   for k = 1:length(fname)   
      existmex = exist([fname{k},'.',ext]);
      if (existmex ~= 3)
         cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]); 
      end
   end      
   existmexeig = exist(['mexeig','.',ext]);
   if (existmexeig ~= 3)
      if strcmp(computer_model,'PCWIN') | strcmp(computer_model,'PCWIN64')      
         cmd([mexcmd,'mexeig','  ','mexeigwin','.c  ',libstr]);
      else
         cmd([mexcmd,'mexeig','  ','mexeig','.c  ',libstr]);
      end
   end
   fprintf ('\n Compilation of mexFunctions completed.\n'); 
   cd .. 
   cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
