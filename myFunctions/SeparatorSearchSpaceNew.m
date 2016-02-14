function Sij = SeparatorSearchSpaceNew(nei,nej,i,j,H,pcType,eta)


if (length(nei) < length(nej))
    ne = [i j mysetdiff(nei,j)];
else
    ne = [i j mysetdiff(nej,i)];
end

if (length(ne)-1 < eta)
    Sij = [];
    Sij2 = [];
    return;
end

Sij = setdiff(ne,[i j]);
% %if pcType == 2 || eta == 1
% if pcType == 2
%     Sij = setdiff(ne,[i j]);
%     %Sij2 = Sij;
%     Sij2 = [];
%     return;
% end
% 
% % if nchoosek(ne,eta) is small enough, don't both computing marginal
% % graphs
% 
% % if (nchoosek(length(ne)-1,eta) < 50)
% %     Sij = setdiff(ne,[i j]);
% %     %Sij2 = Sij;
% %     Sij2 = Sij;
% %     return;
% % end
% 
% %Sij2 = setdiff(ne,[i j]);
% 
% %display('here');
% % Find marginal graph of H_ne
% % Note that the indexing will change since Hm 
% % is a matrix of different size
% 
% Hm = MarginalGraph(H,ne); Hm(1,2) = 0; Hm(2,1) = 0;
% 
% % Sij = all nodes connected to both i and j in Hm
% 
% ni = FindNeighborsUndirected(Hm,1,1);
% nj = FindNeighborsUndirected(Hm,2,1);
% Sij = ne(myintersect(ni,nj));
% Sij2 = [];

end