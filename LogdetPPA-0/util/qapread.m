%%********************************************************
%%
%% [A,B] = qapread('~/ZXY/PAL/QAP/qapdata/nug12.dat'); 
%%
%%********************************************************

  function [A,B] = qapread(fname); 

  fid = fopen(fname,'r');
  if (fid == -1); error('file cannot be opened'); end

  [datavec,count] = fscanf(fid,'%c');
  fclose(fid); 
  linefeeds = findstr(datavec,char(10));
  datavec(linefeeds) = blanks(length(linefeeds)); 
  datavec = sscanf(datavec,'%f'); 

  n = datavec(1); n2 = n*n; 
  aa = datavec(2:n2+1); 
  bb = datavec(n2+2:2*n2+1); 
  
  A = reshape(aa,n,n);
  B = reshape(bb,n,n); 
%%********************************************************
