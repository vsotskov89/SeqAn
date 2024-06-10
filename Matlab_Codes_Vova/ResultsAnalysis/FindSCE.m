function SCEs=FindSCE(Neurons,Experiment,Params,results)

%% Find Correlated pairs of Neurons
disp('    Finding correlated pairs : in progress')
    Neurons=FindCorrPairs(Neurons,Experiment,results);
fprintf([repmat('\b',1,12) '%s\n'],'Done')
%% Get SCEs Timings (starts/ends)
disp('    Finding SCEs timings : in progress')
    SCEs=[];
    SCEs=GetSCEsTimings(SCEs,Neurons,Experiment,Params,'Dir1');
    SCEs=GetSCEsTimings(SCEs,Neurons,Experiment,Params,'Dir2');
fprintf([repmat('\b',1,12) '%s\n'],'Done')
%% Plot SCEs along animal position
disp('    Plotting SCEs along the position : in progress')
    plotSCEsPosition(SCEs,Experiment,Params)
fprintf([repmat('\b',1,12) '%s\n'],'Done')
%% Plot SCEs Summaries
if Params.SCEs.PrintPDFs
    disp('    Filtering Fluo Traces : in progress')
        sampleAct=[(1:Experiment.SizeData)' results.C_raw'];                                                        % Data to FMA Sample
        FiltAct = FilterLFP(sampleAct,'passband',12.5*Params.Method.filtRange); FiltAct=FiltAct(:,2:end)';          % Filter then Sample to Data
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    disp('    Printing summary PDFs for Dir. 1 : SCE 001')
        plotSCEsSummaries(SCEs,Neurons,results,Experiment,Params,FiltAct,'Dir1')
    fprintf([repmat('\b',1,11) '%s\n'],'Done')
    disp('    Printing summary PDFs for Dir. 2 : SCE 001')    
        plotSCEsSummaries(SCEs,Neurons,results,Experiment,FiltAct,'Dir2')
    fprintf([repmat('\b',1,11) '%s\n'],'Done')   
end


end
