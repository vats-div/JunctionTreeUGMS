% Select Best NLasso Method for each case

a = 3;
for i = 1:NumI
    for indn = 1:length(N)
        
        % find number of edges in Ghat_JT_EBIC_GL
        Ghat_JT_EBIC_GL = NewSave_Ghat_JT_EBIC_GL{i,indn};
        numEdges = sum(Ghat_JT_EBIC_GL(:));


        nn = zeros(1,4);
        g{1} = NewSave_Ghat_JT_EBIC25_NL{i,indn};
        g{2} = NewSave_Ghat_JT_EBIC50_NL{i,indn};
        g{3} = NewSave_Ghat_JT_EBIC75_NL{i,indn};
        g{4} = NewSave_Ghat_JT_EBIC100_NL{i,indn};
        nn(1) = sum(g{1}(:));
        nn(2) = sum(g{2}(:));
        nn(3) = sum(g{3}(:));
        nn(4) = sum(g{4}(:));
        diffE = abs(nn - numEdges);
        [diffEMin, indMin] = min(diffE);
        NewSave_Ghat_JT_EBIC_NL{i,indn} = g{indMin}; %#ok<*SAGROW>
        
        numEdges = nn(indMin);
        nn = zeros(1,4);
        g{1} = NewSave_Ghat_NOJT_EBIC25_NL{i,indn};
        g{2} = NewSave_Ghat_NOJT_EBIC50_NL{i,indn};
        g{3} = NewSave_Ghat_NOJT_EBIC75_NL{i,indn};
        g{4} = NewSave_Ghat_NOJT_EBIC100_NL{i,indn};
        nn(1) = sum(g{1}(:));
        nn(2) = sum(g{2}(:));
        nn(3) = sum(g{3}(:));
        nn(4) = sum(g{4}(:));
        diffE = -(nn - numEdges); % choose a positive one
        if sum(diffE > 0) == 4 % if all positive
            [diffEMin, indMin] = min(diffE);
            NewSave_Ghat_NOJT_EBIC_NL{i,indn} = g{indMin};
        else
            diffE(diffE > 0) = 100000;
            [diffEMin, indMin] = min(abs(diffE));
            NewSave_Ghat_NOJT_EBIC_NL{i,indn} = g{indMin};
        end
    end
end

NumEdges_JT_EBIC_PC= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_PC));
NumEdges_JT_EBIC_GL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_GL));
NumEdges_JT_EBIC_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_JT_EBIC_NL));
NumEdges_NOJT_EDGES_GL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EDGES_GL));
NumEdges_NOJT_EDGES_PC= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EDGES_PC));
NumEdges_NOJT_EBIC_NL= mean(cellfun(numEdgesfunc,NewSave_Ghat_NOJT_EBIC_NL));


EDMatrix_Ghat_JT_EBIC_NL = mean(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_NL));
EDMatrix_Ghat_NOJT_EBIC_NL = mean(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC_NL));
TPMatrix_Ghat_JT_EBIC_NL = mean(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_NL));
TPMatrix_Ghat_NOJT_EBIC_NL = mean(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC_NL));
FPMatrix_Ghat_JT_EBIC_NL = mean(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_NL));
FPMatrix_Ghat_NOJT_EBIC_NL = mean(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC_NL));
FractionWeakEdges_JT_EBIC_NL = mean(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_NL));
FractionWeakEdges_NOJT_EBIC_NL = mean(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC_NL));
stdEDMatrix_Ghat_JT_EBIC_NL = std(cellfun(EDfunc,NewSave_Ghat_JT_EBIC_NL));
stdEDMatrix_Ghat_NOJT_EBIC_NL = std(cellfun(EDfunc,NewSave_Ghat_NOJT_EBIC_NL));
stdTPMatrix_Ghat_JT_EBIC_NL = std(cellfun(TPfunc,NewSave_Ghat_JT_EBIC_NL));
stdTPMatrix_Ghat_NOJT_EBIC_NL = std(cellfun(TPfunc,NewSave_Ghat_NOJT_EBIC_NL));
stdFPMatrix_Ghat_JT_EBIC_NL = std(cellfun(FPfunc,NewSave_Ghat_JT_EBIC_NL));
stdFPMatrix_Ghat_NOJT_EBIC_NL = std(cellfun(FPfunc,NewSave_Ghat_NOJT_EBIC_NL));
stdFractionWeakEdges_JT_EBIC_NL = std(cellfun(fweakEdges,NewSave_Ghat_JT_EBIC_NL));
stdFractionWeakEdges_NOJT_EBIC_NL = std(cellfun(fweakEdges,NewSave_Ghat_NOJT_EBIC_NL));

% Make Table

% for n = 

JgL = ['JgL '];
gL = ['gL '];
JPC = ['JPC '];
PC = ['PC '];
JNL = ['JNL '];
NL = ['NL '];

JgL = ['$\CHo$ & \multirow{6}{*}{} $300$ & $\JgL$ '];
gL = ['{$p = 100$} & & $\gL$ '];
JPC = ['\multirow{4}{*}{} && $\JfPC$ '];
PC = ['&& $\fPC$ '];
JNL = ['&&$\JnL$ '];
NL = ['&&$\nL$'];
OrG = ['&&$Oracle$'];
OrP = ['&&$Oracle$'];

for i = [1]
    JgL = [JgL ' & $' num2str(FractionWeakEdges_JT_EBIC_GL(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_GL(i),a) '$ ' ];
    gL = [gL ' & $' num2str(FractionWeakEdges_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EDGES_GL(i),a) '$ ' ];
    JPC = [JPC ' & $' num2str(FractionWeakEdges_JT_EBIC_PC(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_PC(i),a) '$ ' ];
    PC = [PC ' & $' num2str(FractionWeakEdges_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EDGES_PC(i),a) '$ ' ];
    JNL = [JNL ' & $' num2str(FractionWeakEdges_JT_EBIC_NL(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_NL(i),a) '$ ' ];
    NL = [NL ' & $' num2str(FractionWeakEdges_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EBIC_NL(i),a) '$ ' ];
end

JgL = [JgL '\\'];
gL = [gL '\\'];
JPC = [JPC '\\'];
PC = [PC '\\'];
JNL = [JNL '\\'];
NL = [NL '\\'];
OrG = [OrG '\\'];
OrP = [OrP '\\'];


disp(JgL);
disp(gL);
disp(JPC);
disp(PC);
disp(JNL);
disp(NL);

JgL = ['$\CHo$ & \multirow{6}{*}{} $300$ & $\JgL$ '];
gL = ['{$p = 100$} & & $\gL$ '];
JPC = ['\multirow{4}{*}{} && $\JfPC$ '];
PC = ['&& $\fPC$ '];
JNL = ['&&$\JnL$ '];
NL = ['&&$\nL$'];
OrG = ['&&$Oracle$'];
OrP = ['&&$Oracle$'];

for i = [2]
    JgL = [JgL ' & $' num2str(FractionWeakEdges_JT_EBIC_GL(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_GL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_GL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_GL(i),a) '$ ' ];
    gL = [gL ' & $' num2str(FractionWeakEdges_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EDGES_GL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EDGES_GL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EDGES_GL(i),a) '$ ' ];
    JPC = [JPC ' & $' num2str(FractionWeakEdges_JT_EBIC_PC(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_PC(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_PC(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_PC(i),a) '$ ' ];
    PC = [PC ' & $' num2str(FractionWeakEdges_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EDGES_PC(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EDGES_PC(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EDGES_PC(i),a) '$ ' ];
    JNL = [JNL ' & $' num2str(FractionWeakEdges_JT_EBIC_NL(i),a) '$ ($' num2str(stdFractionWeakEdges_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_JT_EBIC_NL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_JT_EBIC_NL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_JT_EBIC_NL(i),a) '$ ' ];
    NL = [NL ' & $' num2str(FractionWeakEdges_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdFractionWeakEdges_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(FPMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdFPMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(TPMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdTPMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' ...
        num2str(EDMatrix_Ghat_NOJT_EBIC_NL(i),a) '$ ($' num2str(stdEDMatrix_Ghat_NOJT_EBIC_NL(i)/sqrt(n),a) '$) & $' num2str(NumEdges_NOJT_EBIC_NL(i),a) '$ ' ];
end


JgL = [JgL '\\'];
gL = [gL '\\'];
JPC = [JPC '\\'];
PC = [PC '\\'];
JNL = [JNL '\\'];
NL = [NL '\\'];
OrG = [OrG '\\'];
OrP = [OrP '\\'];


disp(JgL);
disp(gL);
disp(JPC);
disp(PC);
disp(JNL);
disp(NL);

