function [Params,Neurons,Experiment,results,currentResultsFile] = ComputePlaceActivity(Params,Experiment,currentResultsFile,results)
% Compute Place-related activity from analysis results for each neurons %

%% Load Data
    [Experiment,results,currentResultsFile] = LoadData(Experiment,currentResultsFile,results);
    Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),'Ends',cell(Experiment.nNeurons,1),...
        'Duration',cell(Experiment.nNeurons,1),'HasTheta',cell(Experiment.nNeurons,1),'EvtsCat',cell(Experiment.nNeurons,1),'inPF',cell(Experiment.nNeurons,1)),...
        'PlaceActivity',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'SummedPlaceActivity',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'PFs',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'Conc',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'Stab',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'Peak',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)),...
        'isPC',struct('Dir1',cell(Experiment.nNeurons,1),'Dir2',cell(Experiment.nNeurons,1)));
    Params.Obsd.Time = size(results.C,2); % Number of time frames. 

        
%%  Get Position-related Data
    % Get instant position and speed of the animal
    disp('    Calculating position and speed : in progress')
        [Experiment,results,Params] = GetPositionSpeed(Params,Experiment,results);
    fprintf([repmat('\b',1,12) '%s\n'],'Done')

    % Get Direction Segments
    disp('    Finding Direction Segments : in progress')
        Experiment = GetDirSegments(Params,Experiment,'Dir1');                  % Segments where the animal run from left to right (Dir. 1)
        Experiment = GetDirSegments(Params,Experiment,'Dir2');                  % Segments where the animal run from right to left (Dir. 2)
        Experiment = GetDirSegments(Params,Experiment,'Still');                 % Segments where the animal do not run (Stillness)
    fprintf([repmat('\b',1,12) '%s\n'],'Done')

    disp('    Merging Direction Segments : in progress')
    Experiment=MergeSegments(Experiment);
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
    
    % Check data
    figure; 
    ax1 = subplot(2,1,1);
    toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
    toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
    plot(Experiment.positionXcmSmooth); hold on; plot(toPlot1); plot(toPlot2)         
    ax2 = subplot(2,1,2);
    plot(Experiment.mouseSpeed); hold on; yline(2); yline(-2);
    linkaxes([ax1,ax2],'x');

    % Get time spent in each spatial bin (for normalization)
    disp('    Calculating time spent in each position : in progress')
        Experiment = GetTimeSpent(Experiment,'Dir1');
        Experiment = GetTimeSpent(Experiment,'Dir2');
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
    % Get speed Map
    disp('    Calculating speed matrix : in progress')
        Experiment = getSpeedMap(Experiment,'Dir1');
        Experiment = getSpeedMap(Experiment,'Dir2');
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
%% Find events in each neuron's activity
    disp('    Finding neuronal events : in progress')
        Neurons = GetRises(results,Params,Experiment,Neurons);   
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
%% Compute place-related activity
    disp('    Computing place-related activity : in progress')
        [Neurons,Experiment] = GetActivity(Neurons,Experiment,'Dir1',Params, results);
        [Neurons,Experiment] = GetActivity(Neurons,Experiment,'Dir2',Params, results);
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
%% "Classical" format data
    AllPADir1 = {Neurons.PlaceActivity.Dir1}; AllPADir1=permute(cat(3,AllPADir1{:}),[3,1,2]);
    AllPADir2 = {Neurons.PlaceActivity.Dir2}; AllPADir2=permute(cat(3,AllPADir2{:}),[3,1,2]);
    AllMeanPADir1 = cell2mat({Neurons.SummedPlaceActivity.Dir1}');
    AllMeanPADir2 = cell2mat({Neurons.SummedPlaceActivity.Dir2}');

%% Plot
% figure('Position',[400 100 700 600])
% % Plot sorted Mean Activities 
%     subplot(10,1,1:7)
%         [~,MaxDir1]  = max(AllMeanPADir1,[],2); [~,sortedMaxDir1]=sort(MaxDir1);
%         imagesc(AllMeanPADir1(sortedMaxDir1,:)./max(AllMeanPADir1(sortedMaxDir1,:),[],2))
%         xlabel('Position (cm)'); ylabel('Neuron #')
%         title('Sorted Normalized Place-related Activities')
%  
% % Check Dir1 and Dir2 segments 
%     subplot(10,1,9:10)
%         toPlot1=Experiment.positionXcmSmooth; toPlot1(~Experiment.Dir1.Segments)=nan;
%         toPlot2=Experiment.positionXcmSmooth; toPlot2(~Experiment.Dir2.Segments)=nan;
%         plot([Experiment.positionXcmSmooth,toPlot1,toPlot2])
%         xlabel('Time (frame #)'); ylabel('Position (cm)')
%         title('Segments of Dir1 and Dir2 along position');legend({'Position','Dir 1','Dir 2'},'Location','eastoutside')
end















