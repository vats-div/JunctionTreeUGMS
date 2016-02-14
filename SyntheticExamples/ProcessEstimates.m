% Processes files generated using scripts

a = 4;
p1 = perc * p;
for i = 1:NumI
    for indn = 1:length(N)
       NewSave_Ghat_NOJT_Oracle_GL{i,indn} = full(NewSave_Ghat_NOJT_Oracle_GL{i,indn});
       NewSave_Ghat_NOJT_EBIC_GL{i,indn} = full(NewSave_Ghat_NOJT_EBIC_GL{i,indn});
       NewSave_Ghat_NOJT_EBIC_PC{i,indn} = full(NewSave_Ghat_NOJT_EBIC_PC{i,indn});
       NewSave_Ghat_NOJT_Oracle_PC{i,indn} = full(NewSave_Ghat_NOJT_Oracle_PC{i,indn});
        NewSave_Ghat_NOJT_EBIC25_NL{i,indn} = full(Save_Ghat_NOJT_EBIC25_NL{i,indn});
        NewSave_Ghat_NOJT_EBIC50_NL{i,indn} = full(Save_Ghat_NOJT_EBIC50_NL{i,indn});
        NewSave_Ghat_NOJT_EBIC75_NL{i,indn} = full(Save_Ghat_NOJT_EBIC75_NL{i,indn});
        NewSave_Ghat_NOJT_EBIC100_NL{i,indn} = full(Save_Ghat_NOJT_EBIC100_NL{i,indn});
        NewSave_Ghat_JT_EBIC_GL{i,indn} = full(Save_Ghat_JT_EBIC_GL{i,indn});
        NewSave_Ghat_NOJT_EDGES_GL{i,indn} = full(Save_Ghat_NOJT_EDGES_GL{i,indn});
        NewSave_Ghat_JT_EBIC25_NL{i,indn} = full(Save_Ghat_JT_EBIC25_NL{i,indn});
        NewSave_Ghat_JT_EBIC50_NL{i,indn} = full(Save_Ghat_JT_EBIC50_NL{i,indn});
        NewSave_Ghat_JT_EBIC75_NL{i,indn} = full(Save_Ghat_JT_EBIC75_NL{i,indn});
        NewSave_Ghat_JT_EBIC100_NL{i,indn} = full(Save_Ghat_JT_EBIC100_NL{i,indn});
        NewSave_Ghat_JT_EBIC_PC{i,indn} = full(Save_Ghat_JT_EBIC_PC{i,indn});
        NewSave_Ghat_NOJT_EDGES_PC{i,indn} = full(Save_Ghat_NOJT_EDGES_PC{i,indn});
    end
end

Gweak = zeros(p,p); Gweak(1:p1,1:p1) = 1; Gweak = Gweak .* Gtrue;
NumWeakEdges = sum(Gweak(:));
fweakEdges = @(G) (sum(sum(Gweak.*G))/NumWeakEdges);
FPfunc = @(cellclus) FalsePositiveGraph(Gtrue,cellclus,1);
TPfunc = @(cellclus) TruePositiveGraph(Gtrue,cellclus,1);
EDfunc = @(cellclus) CompareGraphs(Gtrue,cellclus);
numEdgesfunc = @(cellclus) (sum(cellclus(:))/2);

NumEdges_NOJT_EBIC_GL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC_GL));
NumEdges_NOJT_EBIC25_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC25_NL));
NumEdges_NOJT_EBIC50_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC50_NL));
NumEdges_NOJT_EBIC75_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC75_NL));
NumEdges_NOJT_EBIC100_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC100_NL));
NumEdges_JT_EBIC25_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC25_NL));
NumEdges_JT_EBIC50_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC50_NL));
NumEdges_JT_EBIC75_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC75_NL));
NumEdges_JT_EBIC100_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC100_NL));
NumEdges_JT_EBIC_PC= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_PC));
NumEdges_NOJT_EDGES_PC= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EDGES_PC));
NumEdges_JT_EBIC_GL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_GL));
NumEdges_NOJT_EDGES_GL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EDGES_GL));



FractionWeakEdges_NOJT_Oracle_GL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_Oracle_GL));
FractionWeakEdges_NOJT_Oracle_PC = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_Oracle_PC));
FractionWeakEdges_JT_EBIC_GL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_GL));
FractionWeakEdges_NOJT_EDGES_GL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EDGES_GL));
FractionWeakEdges_JT_EBIC_PC = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_PC));
FractionWeakEdges_NOJT_EDGES_PC = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EDGES_PC));
FractionWeakEdges_NOJT_EBIC25_NL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC25_NL));
FractionWeakEdges_NOJT_EBIC50_NL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC50_NL));
FractionWeakEdges_NOJT_EBIC75_NL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC75_NL));
FractionWeakEdges_NOJT_EBIC100_NL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC100_NL));
FractionWeakEdges_JT_EBIC25_NL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC25_NL));
FractionWeakEdges_JT_EBIC50_NL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC50_NL));
FractionWeakEdges_JT_EBIC75_NL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC75_NL));
FractionWeakEdges_JT_EBIC100_NL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC100_NL));

TPMatrix_Ghat_NOJT_Oracle_GL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_Oracle_GL));
TPMatrix_Ghat_NOJT_Oracle_PC = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_Oracle_PC));
TPMatrix_Ghat_JT_EBIC_GL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_GL));
TPMatrix_Ghat_NOJT_EDGES_GL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EDGES_GL));
TPMatrix_Ghat_JT_EBIC_PC = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_PC));
TPMatrix_Ghat_NOJT_EDGES_PC = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EDGES_PC));
TPMatrix_Ghat_NOJT_EBIC25_NL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC25_NL));
TPMatrix_Ghat_NOJT_EBIC50_NL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC50_NL));
TPMatrix_Ghat_NOJT_EBIC75_NL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC75_NL));
TPMatrix_Ghat_NOJT_EBIC100_NL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC100_NL));
TPMatrix_Ghat_JT_EBIC25_NL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC25_NL));
TPMatrix_Ghat_JT_EBIC50_NL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC50_NL));
TPMatrix_Ghat_JT_EBIC75_NL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC75_NL));
TPMatrix_Ghat_JT_EBIC100_NL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC100_NL));

FPMatrix_Ghat_NOJT_Oracle_GL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_Oracle_GL));
FPMatrix_Ghat_NOJT_Oracle_PC = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_Oracle_PC));
FPMatrix_Ghat_JT_EBIC_GL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_GL));
FPMatrix_Ghat_NOJT_EDGES_GL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EDGES_GL));
FPMatrix_Ghat_JT_EBIC_PC = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_PC));
FPMatrix_Ghat_NOJT_EDGES_PC = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EDGES_PC));
FPMatrix_Ghat_NOJT_EBIC25_NL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC25_NL));
FPMatrix_Ghat_NOJT_EBIC50_NL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC50_NL));
FPMatrix_Ghat_NOJT_EBIC75_NL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC75_NL));
FPMatrix_Ghat_NOJT_EBIC100_NL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC100_NL));
FPMatrix_Ghat_JT_EBIC25_NL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC25_NL));
FPMatrix_Ghat_JT_EBIC50_NL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC50_NL));
FPMatrix_Ghat_JT_EBIC75_NL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC75_NL));
FPMatrix_Ghat_JT_EBIC100_NL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC100_NL));

EDMatrix_Ghat_NOJT_Oracle_GL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_Oracle_GL));
EDMatrix_Ghat_NOJT_Oracle_PC = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_Oracle_PC));
EDMatrix_Ghat_JT_EBIC_GL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_GL));
EDMatrix_Ghat_NOJT_EDGES_GL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EDGES_GL));
EDMatrix_Ghat_JT_EBIC_PC = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_PC));
EDMatrix_Ghat_NOJT_EDGES_PC = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EDGES_PC));
EDMatrix_Ghat_NOJT_EBIC25_NL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC25_NL));
EDMatrix_Ghat_NOJT_EBIC50_NL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC50_NL));
EDMatrix_Ghat_NOJT_EBIC75_NL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC75_NL));
EDMatrix_Ghat_NOJT_EBIC100_NL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC100_NL));
EDMatrix_Ghat_JT_EBIC25_NL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC25_NL));
EDMatrix_Ghat_JT_EBIC50_NL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC50_NL));
EDMatrix_Ghat_JT_EBIC75_NL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC75_NL));
EDMatrix_Ghat_JT_EBIC100_NL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC100_NL));

stdNumEdges_JT_EBIC_GL = std(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_GL));
stdNumEdges_NOJT_EBIC_GL= std(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC_GL));
stdNumEdges_NOJT_EBIC25_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC25_NL));
stdNumEdges_NOJT_EBIC50_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC50_NL));
stdNumEdges_NOJT_EBIC75_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC75_NL));
stdNumEdges_NOJT_EBIC100_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC100_NL));
stdNumEdges_JT_EBIC25_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC25_NL));
stdNumEdges_JT_EBIC50_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC50_NL));
stdNumEdges_JT_EBIC75_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC75_NL));
stdNumEdges_JT_EBIC100_NL= std(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC100_NL));

stdFractionWeakEdges_NOJT_Oracle_GL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_Oracle_GL));
stdFractionWeakEdges_NOJT_Oracle_PC = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_Oracle_PC));
stdFractionWeakEdges_JT_EBIC_GL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_GL));
stdFractionWeakEdges_NOJT_EDGES_GL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EDGES_GL));
stdFractionWeakEdges_JT_EBIC_PC = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_PC));
stdFractionWeakEdges_NOJT_EDGES_PC = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EDGES_PC));

stdFractionWeakEdges_NOJT_EBIC25_NL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC25_NL));
stdFractionWeakEdges_NOJT_EBIC50_NL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC50_NL));
stdFractionWeakEdges_NOJT_EBIC75_NL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC75_NL));
stdFractionWeakEdges_NOJT_EBIC100_NL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC100_NL));
stdFractionWeakEdges_JT_EBIC25_NL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC25_NL));
stdFractionWeakEdges_JT_EBIC50_NL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC50_NL));
stdFractionWeakEdges_JT_EBIC75_NL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC75_NL));
stdFractionWeakEdges_JT_EBIC100_NL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC100_NL));

stdFPMatrix_Ghat_NOJT_Oracle_GL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_Oracle_GL));
stdFPMatrix_Ghat_NOJT_Oracle_PC = std(cellfun(FPfunc,NewSave_Ghat_NOJT_Oracle_PC));
stdFPMatrix_Ghat_JT_EBIC_GL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_GL));
stdFPMatrix_Ghat_NOJT_EDGES_GL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EDGES_GL));
stdFPMatrix_Ghat_JT_EBIC_PC = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_PC));
stdFPMatrix_Ghat_NOJT_EDGES_PC = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EDGES_PC));
stdFPMatrix_Ghat_NOJT_EBIC25_NL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC25_NL));
stdFPMatrix_Ghat_NOJT_EBIC50_NL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC50_NL));
stdFPMatrix_Ghat_NOJT_EBIC75_NL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC75_NL));
stdFPMatrix_Ghat_NOJT_EBIC100_NL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC100_NL));
stdFPMatrix_Ghat_JT_EBIC25_NL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC25_NL));
stdFPMatrix_Ghat_JT_EBIC50_NL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC50_NL));
stdFPMatrix_Ghat_JT_EBIC75_NL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC75_NL));
stdFPMatrix_Ghat_JT_EBIC100_NL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC100_NL));

stdTPMatrix_Ghat_NOJT_Oracle_GL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_Oracle_GL));
stdTPMatrix_Ghat_NOJT_Oracle_PC = std(cellfun(TPfunc,NewSave_Ghat_NOJT_Oracle_PC));
stdTPMatrix_Ghat_JT_EBIC_GL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_GL));
stdTPMatrix_Ghat_NOJT_EDGES_GL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EDGES_GL));
stdTPMatrix_Ghat_JT_EBIC_PC = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_PC));
stdTPMatrix_Ghat_NOJT_EDGES_PC = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EDGES_PC));
stdTPMatrix_Ghat_NOJT_EBIC25_NL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC25_NL));
stdTPMatrix_Ghat_NOJT_EBIC50_NL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC50_NL));
stdTPMatrix_Ghat_NOJT_EBIC75_NL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC75_NL));
stdTPMatrix_Ghat_NOJT_EBIC100_NL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC100_NL));
stdTPMatrix_Ghat_JT_EBIC25_NL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC25_NL));
stdTPMatrix_Ghat_JT_EBIC50_NL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC50_NL));
stdTPMatrix_Ghat_JT_EBIC75_NL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC75_NL));
stdTPMatrix_Ghat_JT_EBIC100_NL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC100_NL));

stdEDMatrix_Ghat_NOJT_Oracle_GL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_Oracle_GL));
stdEDMatrix_Ghat_NOJT_Oracle_PC = std(cellfun(EDfunc,NewSave_Ghat_NOJT_Oracle_PC));
stdEDMatrix_Ghat_JT_EBIC_GL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_GL));
stdEDMatrix_Ghat_NOJT_EDGES_GL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EDGES_GL));
stdEDMatrix_Ghat_JT_EBIC_PC = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_PC));
stdEDMatrix_Ghat_NOJT_EDGES_PC = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EDGES_PC));
stdEDMatrix_Ghat_NOJT_EBIC25_NL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC25_NL));
stdEDMatrix_Ghat_NOJT_EBIC50_NL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC50_NL));
stdEDMatrix_Ghat_NOJT_EBIC75_NL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC75_NL));
stdEDMatrix_Ghat_NOJT_EBIC100_NL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC100_NL));
stdEDMatrix_Ghat_JT_EBIC25_NL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC25_NL));
stdEDMatrix_Ghat_JT_EBIC50_NL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC50_NL));
stdEDMatrix_Ghat_JT_EBIC75_NL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC75_NL));
stdEDMatrix_Ghat_JT_EBIC100_NL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC100_NL));
