%% Master function that call every other functions
Params=GetParams;

for exp=1:numel(Params.Paths)
    if ~exist('currentResultsFile','var'); currentResultsFile=[]; results=[] ; end      % Initialize a variable used later to avoid the (long) reloading of results.mat if it does not exist
    
    disp(['Starting analysis : experiment ' strrep(Params.Paths{exp},[Params.DataFolder '/'],'')])
        Experiment=[];
        Experiment.path=[Params.Paths{exp} '/'];                                            % Select one experiment
    
    % Call ComputePlaceActivity.m 
    disp('Start Computing Place-related Activities')
        [Params,Neurons,Experiment,results,currentResultsFile]=ComputePlaceActivity(Params,Experiment,currentResultsFile,results);
    
    % Call FindPlaceCells.m
    disp('Start Finding Place Cells')
        [Experiment,Neurons,Params]=FindPlaceCells(Experiment,Neurons,Params,results);
     
    % If SCE Test is wanted, do the SCE detection
    if Params.findSCE
    disp('Start Finding SCEs')
        SCEs=FindSCE(Neurons,Experiment,Params,results);
    end
   
    % If Theta Test is wanted, do the theta detection
    if Params.findTheta
    disp('Start Looking for theta')
        Neurons=FindTheta(Neurons,Experiment,results,Params);
    end
    
    % If Event statistics are wanted (to make the .xls for example)
%     if Params.getEvtsStats || Params.MakeXLS
%         disp('Recovering events statistics')
%         [Neurons,Stats]=GetEventsStats(Neurons,Experiment,results,Params);
%         save([Experiment.path 'Stats'], 'Stats')
%     end
    
    % If a .xls summary is wanted
    if Params.MakeXLS
        disp('Writing summary .xls')
        Params=makeXLS(Neurons,Experiment,Params,exp);
    end
    
    % If saving Matrices is wanted
    if Params.SaveMats
        disp('Saving Matrices')
        save([Experiment.path 'Neurons'], 'Neurons')
        save([Experiment.path 'Params'], 'Params')
        save([Experiment.path 'Experiment'], 'Experiment')
    end
    disp(['Done analyzing experiment ' strrep(Params.Paths{exp},[Params.DataFolder '/'],'')])
    close all
end
disp('Analysis Finished')
