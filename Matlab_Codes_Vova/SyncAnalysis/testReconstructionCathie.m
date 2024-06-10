clear;

% test from scratch including detection of events (to tune parameters).
% on behavior data
% we suppose than a first analysis of the experiment has been done using
% program Master.m  : the structures Neurons and Experiments have been
% generated and can be loaded from Params.DataFolder folder. 

Params.DataFolder = 'G:\Souris154148\22-06-30\130014_sCMOS_154148-awake_FromAllExp\';
Params.Method.ActivityType = 'Spikes';                                                                      % false Params Struct.,
Params.Method.SmallEvents = 'RisesDerivative';                                                                  % Choose the method to extract the rises
Params.RedetectEvents = 1;
nBins=200;                                                                                              % Number of spatial bins (200 --> 1bin = 1cm)
wndw=0.3; %0.02;    

cd(Params.DataFolder)
load('Neurons.mat');
load('Experiment.mat');

% optional redetection of events
if Params.RedetectEvents == 1
    load('Results.mat');
    NeuronsB = GetRises(results,Params,Experiment,Neurons);
else
    NeuronsB = Neurons;
end


Dir1=[Experiment.Dir1.Starts Experiment.Dir1.Ends]; Dir2=[Experiment.Dir2.Starts Experiment.Dir2.Ends]; % Get Mouvement Periods
Dir1 = Dir1(1:34, :); Dir2 = Dir2(1:33, :);                                                                % Det rid of the end of the behavior 
spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat. spikes = list of (t, neuron) couples
EvtMat=vertcat(NeuronsB.Rises.Matrix);EvtMat=full(EvtMat)>0;                                         % matrix of events. need to convert this a sample data; 
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMat(n,:)); 
        spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    
    [spikesMvt,~,~] = InIntervals(spikes(:,1),sort([Dir1 ; Dir2]));                                         % Get Spikes occuring during mouvement periods
    spikes=spikes(spikesMvt,:);                                                                             % Restrict Spikes to mouvement periods

    positions=[(1:Experiment.SizeData-4)'...                                                                  % Normalize Position between [0 0.5]
        0.5*fillmissing((Experiment.positionXcmSmooth-min(Experiment.positionXcmSmooth))/(max(Experiment.positionXcmSmooth)-min(Experiment.positionXcmSmooth)),'previous')]; %min(Experiment.positionXcmSmooth)
    [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);                                                         % Get positions during Dir2
    positions(Mvtd2,2)=1-positions(Mvtd2,2);                                                                % Linearisation (Dir1 = 0 --> 0.5 ; Dir2 = 0.5 --> 1)
    [Mvt,~,~] = InIntervals(positions(:,1),sort([Dir1 ; Dir2]));                                            % Get positions occuring during mouvement periods
    positions=positions(Mvt,:);                                                                             % Restrict positions to mouvement periods
    %figure; plot(positions(:,1), positions(:,2),'.');
    
    spikes(:,1)=spikes(:,1)/100;                                                                            % Timings (in ms) from indices
    positions(:,1)=positions(:,1)/100;                                                                      % Timings (in ms) from indices
    
   
%    Test the reconstruction on Awake itself (just to evaluate the performance)
    [statsTest,~,~] = ReconstructPosition(positions,spikes,'mode','both','window',wndw,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
    % figure; ax1 = subplot(2,1,1); plot(statsTest.positions(:,1), double(toKeep));  ax2 = subplot(2,1,2); plot(statsTest.positions(:,1), statsTest.positions(:,2), '+'); hold on; plot(positions(:,1), positions(:,2)*200, '.'); linkaxes([ax1,ax2],'x');
    %  estim = estim/nBins;
    %  figure; ax1 = subplot(2,1,1);  plot(estim, '.'); ax2 = subplot(2,1,2);  plot(statsTest.positions(:,2)/nBins, '.'); linkaxes([ax1,ax2],'x');
     figure; plot(estim, '.'); hold on; plot(pos, '.');
     ylabel('Mouse Position (cm)');
     xlabel('Time Points (wnd = 20ms)')
     set(gca, 'FontSize',14);

     % Pour comprendre à quoi correspondent les points temporels de
     % statsTest
%      pos2=statsTest.positions(toKeep,:); %save('pos2.mat','pos2');
%      figure; plot(positions(:,1), positions(:,2), '.'); hold on; plot(pos2(:,1), pos2(:,2)/200,'.'); plot(statsTest.positions(:,1),statsTest.positions(:,2)/200,'.')

     
  % quantification de l'erreur
    error = min(abs(estim-pos),nBins-abs(estim-pos));
    errorCm=error*200/nBins;
    mean(errorCm)
    median(errorCm)
%     mean(errorCm(1:floor(end/2)))
%     median(errorCm(1:floor(end/2)))
    figure; histogram(errorCm)
    figure; histogram(errorCm(1:floor(end/2)))

%%
nBinsA=200;%[25 50 100 200]; %nb bins spatiaux
wndwA=0.02;%[0.01 0.02 0.03 0.05 0.1 0.2 0.5];  % fenetre de temps pour la reconstruction 
medMat=zeros(numel(nBinsA),numel(wndwA));
meanMat=zeros(numel(nBinsA),numel(wndwA));
for nb=1:numel(nBinsA)
    nBins=nBinsA(nb);
    for w =1:numel(wndwA)
        wndw=wndwA(w);
        [stats,lambda,Px] = ReconstructPosition(positions,spikes,'mode','train','window',wndw,'nBins',[nBins 1]); % train ; lambda : champs d'activite de tous les neurones, px : répartition des position (histogramme des positions)
        [~,estim]=max(stats.estimations);
        toKeep=sum(~ismember(stats.estimations,ones(nBins,1)/nBins))>0;
        estim=estim(toKeep)';
        pos=stats.positions(toKeep,2);
        error = min(abs(estim-pos),nBins-abs(estim-pos));
        errorCm=error*200/nBins;
        medMat(nb,w)=mean(errorCm);
        meanMat(nb,w)=median(errorCm);
    end
end
%%
figure
    subplot(1,2,1)
        imagesc(meanMat)
        title('Mean error(cm)'),xlabel('Window Size (s)'),ylabel('Spatial Bin Size (cm)'); xticks(1:numel(wndwA)); yticks(1:numel(nBinsA));
        xticklabels(cellfun(@num2str,num2cell(wndwA),'UniformOutput',false)); yticklabels(cellfun(@num2str,num2cell(200*ones(1,numel(nBinsA))./nBinsA),'UniformOutput',false))
    subplot(1,2,2)
        imagesc(medMat)
        title('Median error(cm)'),xlabel('Window Size (s)'),ylabel('Spatial Bin Size (cm)'); xticks(1:numel(wndwA)); yticks(1:numel(nBinsA));
        xticklabels(cellfun(@num2str,num2cell(wndwA),'UniformOutput',false)); yticklabels(cellfun(@num2str,num2cell(200*ones(1,numel(nBinsA))./nBinsA),'UniformOutput',false))
    
%%
[Mx,MxIdx]=max(squeeze(lambda));
PFs=squeeze(lambda);
[~,sMxIdx]=sort(MxIdx);
figure
imagesc((PFs(:,sMxIdx)./Mx(sMxIdx))')
%%
figure
imagesc(stats.estimations)
hold on
plot(stats.positions(:,2),'r','LineWidth',2)
[~,estim]=max(stats.estimations);
scatter(1:numel(stats.positions(:,2)),estim,'y','filled')

%%
figure
E=stats.estimations;
Emax=max(E);
imagesc(E./Emax)









