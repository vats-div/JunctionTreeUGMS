%%*********************************************************************
%% for use in sdpNAL.m
%%
%%*********************************************************************

    function plotfun(blk,par,resnrm); 

    sig = par.sig; 
    dd  = par.dd;
    posidx = par.posidx; 
%%
    subplot(131)
    if ~isfield(par,'phi')
       len = length(posidx{1}); negidx = setdiff([1:blk{1,2}],posidx{1}); 
       if (len > 0)
          semilogy([1:len],sort(abs(dd{1}(posidx{1}))),'ob'); hold on;
       end
       if ~isempty(negidx)
          semilogy(negidx,abs(dd{1}(negidx)),'*r'); 
       end
       title(['eig(X-sig*Z),  sig=',sprintf('%2.1e',sig)]); 
       hold off; 
    else
       semilogy(sort(par.phi{1}),'*'); 
       title(['eig(Phi),  sig=',sprintf('%2.1e',sig)]); 
    end
%%
    subplot(132)
    if isfield(par,'Dsch12')
       r = length(posidx{1}); n = blk{1,2}(1); 
       if (r < n)
          semilogy(sort(abs(par.Dsch12{1}(:))),'*g');             
          title('Dsch12')
          if (r > 0)
   	     axis([1,r*(n-r),1e-8,1])
      	     set(gca,'XTick',[1,r*(n-r)]); 
          end
       end
    elseif isfield(par,'Dsch')
       semilogy(sort(abs(par.Dsch{1}(:))),'*g');             
       title('Dsch')
    end
    subplot(133)
    semilogy(resnrm,'*');
    axis([1 length(resnrm) [0.01, 1.2]*resnrm(1)])
    pause(0.1)
%%*********************************************************************
