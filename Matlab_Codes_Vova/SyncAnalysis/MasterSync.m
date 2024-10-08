%% Path of the Exp. Folder + Channels for TTL and HPC ePhy (Exp. Dep.) + Nb of sub Exps ([SleepPRE, Awake, SleepPOST])
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris154147'; ePhysChan=8; TTLChan=18;
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris154148'; ePhysChan=8; TTLChan=18;nbExp=[1,1,1];
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris133656'; ePhysChan=8; TTLChan=19;nbExp=[2,1,1];
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris154148_22-08-04'; ePhysChan=8; TTLChan=18;nbExp=[0,1,1];
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris154148_22-07-07'; ePhysChan=6; TTLChan=18;nbExp=[1,2,1];
% MainPath='/mnt/soma-home/theo/tmpTheo/Data/Souris133656_21-03-02'; ePhysChan=9; TTLChan=19;nbExp=[1,1,2];

%MainPath='D:\Analysis\Souris154148_22-08-04'; ePhysChan=8; TTLChan=18;nbExp=[0,1,1];
%MainPath='D:\Analysis\Souris133656_21-03-08'; ePhysChan=9; TTLChan=19;nbExp=[2,1,1];
%MainPath='D:\Analysis\Souris154148_22-06-30\'; ePhysChan=8; TTLChan=18;nbExp=[1,1,1];
%MainPath='D:\Analysis\Souris154148_22-06-30\'; ePhysChan=8; TTLChan=18;nbExp=[1,1,0];
%MainPath='D:\Analysis\Souris154147_22-06-21\'; ePhysChan=9; TTLChan=18;nbExp=[0,1,1];
%MainPath='D:\Analysis\Souris133656_21-03-02'; ePhysChan=9; TTLChan=19;nbExp=[1,1,2];
%MainPath='G:\Souris133655\21-02-23\'; ePhysChan=3; TTLChan=19;nbExp=[1,1,0];
MainPath='G:\Souris154148\22-06-30'; ePhysChan=8; TTLChan=18;nbExp=[1,1,1];
%MainPath='G:\Souris154148\22-07-07'; ePhysChan=6; TTLChan=18;nbExp=[1,2,1];
%MainPath='G:\Souris154148\22-08-04'; ePhysChan=6; TTLChan=18;nbExp=[0,1,1];
%MainPath='F:\Souris154147\22-06-29'; ePhysChan=10; TTLChan=18;nbExp=[0,1,1];


DeepCAD='';                                                                 % If you want to use DeepCAD or not. In the concatenated result folder, the matrices must
% DeepCAD='_DeepCAD';                                                         % be named Results.mat (if no DC) and Results_DeepCAD.mat (if DC)
nDownsample = 1;                                                               % N_downsample>1 To analyse data at lower resolution. Requires a file results_XXHz.mat computed using CNMFE
                                                                                 % N_downsample = 1 to analyse data at 100Hz   
%% Find path to files
expPaths=GetExpsPaths('NbExp',nbExp,'Path',MainPath);                                                       % Load the paths to the differents file (or ask the user for them if it wasn't done before)

%% Create low rate tag and timing files at low frame rate if necessary

% create the tag to read results and scmos timing data
% the results file will be ['/Results' DeepCAD tagLowRate '.mat']
% the timing files will be ['/sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']
if nDownsample > 1
    tagLowRate = sprintf( '_ds%d',nDownsample); 
elseif nDownsample == 1
    tagLowRate = '';
else
    return
end

if nDownsample > 1
    createLowRateTimingFiles(expPaths,nDownsample,tagLowRate);                     % create Low Rate timing files for the sCMOS camera in the same folder as original timing files, with the name ['/sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']
end

%% Load or extract the calcium results file
ExtractedResults=LoadExtractedCalciumData(MainPath,expPaths,DeepCAD,nDownsample,tagLowRate);                                       % Load Calcium data (extract it from the concatenated file if it wasn't done before)

%% Load ePhysData (HPC channel and TTL
[file,path] = uigetfile([MainPath '/*.dat'],'Select ePhys data to load'); ePhysFile = fullfile(path,file);  % Select the continuous.dat of the experiment    
ePhysData=loadChanLFP(ePhysChan,[ePhysFile(1:end-3) 'lfp'],1);                                              % Load the LFP (create it from the .dat if it doesn't exist)
TTLData=loadChanDat(TTLChan,ePhysFile); TTLData=single(TTLData);                                            % Load the TTL data

%% Extract TTL Timings
TTLstartsTimes=GetTTLtimes(TTLData,ExtractedResults);                                                       % Get the TTL times (which are also the time of Calcium data points)

%% Get calcium events starts
SizeData = size(ExtractedResults.C, 2); Experiment.SizeData = SizeData;                                     % Build false Experiment Struct.,                
nNeurons = size(ExtractedResults.C, 1); Experiment.nNeurons = nNeurons;

Params.Method.ActivityType = 'Spikes';                                                                      % false Params Struct.,
Params.Method.SmallEvents = 'RisesDerivative';                                                                  % Choose the method to extract the rises
%Params.Method.SmallEvents ='RisesSE';
%Params.Method.SmallEvents = 'RisesSEbsd';
%Params.Method.SmallEvents = 'RisesDecaysSEbsd';

Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),...   % and false Neurons Struct. that are required for GetRises
    'Ends',cell(Experiment.nNeurons,1),'Duration',cell(Experiment.nNeurons,1)));
    
Neurons=GetRises(ExtractedResults,Params,Experiment,Neurons);                                               % Get rises
EvtMat=vertcat(Neurons.Rises.Matrix);EvtMat=full(EvtMat)>0;                                                 % Make a binary matrix (ID,timings)
%nNeurons = 200; 
EvtMat = EvtMat(1:nNeurons, :);                                                                             % keep only the first neurons
EvtStarts=[zeros(nNeurons,1) diff(EvtMat,[],2)]>0;    % change compared to theo (before : Experiment.nNeurons)  % Take only the rises starts

%% Get Number of starts and neurons IDs
wndw=10;                                                                                                    % Size of the moving window (in this version the window is [0 wndw] and not [-wdnw/2 wdnw/2])
[nbStarts,NeuronIDs]=GetNbStarts(EvtStarts,wndw);                                                           % Get the number of starts in each window (here : from t=1 to t=SizeData-wndw)

%% Get chance level of number of starts and set threshold for SCEs
nbPerm=1000;                                                                                                % Number of random distribution of the number of start to compute (100 should be enough but 1000 or 10000 is good)
AllNbStarts=zeros(nbPerm,SizeData-wndw);                                                                    % Initialize the nbStarts matrix (Perm #, timings)
parfor p=1:nbPerm                                                                                              % Loop
    shifts=randi(SizeData,1,nNeurons) ;                                                                         % Get a random time shift for each neuron
    permEvt=cell2mat(arrayfun(@(x) circshift(EvtStarts(x,:),[1 shifts(x)]),(1:numel(shifts))','un',0));         % Apply a circular shift of said shifts to each neuron's event starts
    [AllNbStarts(p,:),~]=GetNbStarts(permEvt,wndw);                                                             % Get nbStarts for this time-shifted matrix of Evt Starts
end
nbStartsThr=prctile(AllNbStarts(:),99.9); %99                                                                    % Set the threshold for a SCE to the 99th percentile of these nbStarts

%% Get SCEs
[putSCE,putSCEidx,putSCEtimes]=GetSCEs(nbStarts,nbStartsThr,EvtStarts,wndw,TTLstartsTimes);                 % Find SCEs (moments where NbStarts exceed the threshold)
sceNeuronIDs=NeuronIDs(putSCEidx);                                                                          % Get the neurons involved in each SCE

%% Compute Spectrogram and get power in ripple band
rangeSpect = [0 500];                                                                                       % Range of the Spectrogram in Hz
windowSpect = 0.05;                                                                                         % Window for computing the spectrogram in s
[spectrogram,t,f] = MTSpectrogram(ePhysData,'range',[0 500],'window',0.05);                                 % Compute the spectrogram of the LFP
bands = SpectrogramBands(spectrogram,f);                                                                    % Define physiological rhythms bands (default for ripples = [100 250] Hz)
ripplesPower=[t bands.ripples];                                                                             % Get the power in the ripples band in each time window

%% PETD of ripple-band power around SCEs
durations=[-1 1];                                                                                           % Window around events times to compute the PETD
GetRipScePETD(ripplesPower,putSCEtimes,durations)                                                           % Get the Peri-Events(SCEs)-Time Distribution of the ripples power

%% Find Ripples
RipBand=[100 250];                                                                                          % Define Ripples Band
ePhysDataFilt=FilterLFP(ePhysData,'passband',RipBand);                                                      % Filter the LFP in the Ripples Band

RipThr = [2 3]; % [2 3]; % [1 2]; % Default values in FindRipples function  [2 5]                                                                                          % Set Threshold for finding the ripples 
[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',RipThr);                                             % Find Ripples
%[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',[0.2 0.4],'durations',[30 15 150]);                 % Weird things happen when finding ripples on G6 data...Thresholds are messed up, probably because of some artifacts so You'll need to get very low thr.
%[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',[0.7 1.4],'durations',[30 20 100]);                 % Weird things happen when finding ripples on G6 data...Thresholds are messed up, probably because of some artifacts so You'll need to get very low thr.

% plot ripples on lfp trace to check that the thesholds are OK
% figure; plot(ePhysData(:,1),ePhysData(:,2));
%  hold on;
%  for kripple = 1 : size(ripples, 1)
%      %kripple = 1;
%      xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
%      yBox = [-20000, 20000, 20000, -20000, -20000];
%      patch(xBox, yBox, 'green', 'FaceAlpha', 0.1, 'FaceColor', 'green', 'FaceAlpha', 0.1);
%  end
% 
%    %plot ripples position on EvtMat matrix
% figure; [yPoints,xPoints] = find(EvtStarts==1); plot(TTLstartsTimes(xPoints),yPoints,'.k');
% hold on;
% for kripple = 1 : size(ripples, 1)
%     %kripple = 1;
%     xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
%     yBox = [0, nNeurons, nNeurons, 0, 0];
%     patch(xBox, yBox, 'green', 'FaceAlpha', 0.1, 'FaceColor', 'green', 'FaceAlpha', 0.1);
% end

   % histogram of ripple duration
% figure; histogram(ripples(:,3)-ripples(:,1)); 

%% CCG of SCEs and SWRs
binSize=0.01;  %0.01; 0.2                                                                                              % Bin size for Cross-CorreloGram
durations=[-1 1];  %[-1 1]; [-20 20];                                                                                          % Window around events times to compute the CCG
[RipSceCCG,tCCG]=GetRipSceCCG(ripples,putSCEtimes,durations,binSize);                                       % Get the CCG

%% Percentage of synchronous SCEs and SWRs
MaxSep=0.1; %0.05; 0.1                                                                                                % Max time distance between the beginnings of two events to consider them synchronous 
[syncSCEs,nbSyncSCEs]=GetNbSync(putSCEtimes,ripples(:,1),MaxSep);                                           % Get the number of SCEs, synced with SWRs
[syncSWRs,nbSyncSWRs]=GetNbSync(ripples(:,1),putSCEtimes,MaxSep);                                           % Get the number of SWRs, synced with SCEs

%% Chance level of synchrony
nbPerm=1000;                                                                                                % Number of shifted SCEs distributions to compute
nbPermSyncSCEs=zeros(1,nbPerm);                                                                             % Initialization of matrices
nbPermSyncSWRs=zeros(1,nbPerm);
permCCG=zeros(nbPerm,size(RipSceCCG,1));
for p=1:nbPerm                                                                                              % Loop
    permSCEidx=rem(putSCEidx+randi(numel(putSCE)),numel(putSCE))+1;                                             % Get a circularly shifted distribution of SCEs
    permSCEtimes=TTLstartsTimes(permSCEidx);                                                                                 % Get the corresponding timings
    [~,nbPermSyncSCEs(p)]=GetNbSync(permSCEtimes,ripples(:,1),MaxSep);                                          % Get the number of SCEs, synced with SWRs
    [~,nbPermSyncSWRs(p)]=GetNbSync(ripples(:,1),permSCEtimes,MaxSep);                                          % Get the number of SWRs, synced with SCEs
    [pccg,~]=GetRipSceCCG(ripples,permSCEtimes,durations,binSize);                                              % Compute CCG
    permCCG(p,:)=pccg(:,1,2);                                                                                   % Store the values
end

thrSCEs=prctile(nbPermSyncSCEs,95);                                                                         % Define thresholds for the chance nb of synchronies as the 95th percentile of the synchronies obtained whith shifted data
thrSWRs=prctile(nbPermSyncSWRs,95);

%% Plot various figs
PlotMisc(nbStarts,sceNeuronIDs,putSCEidx,syncSCEs,syncSWRs,ripples,nNeurons,tCCG,RipSceCCG,nbSyncSCEs,nbPermSyncSCEs,thrSCEs,nbSyncSWRs,nbPermSyncSWRs,thrSWRs,permCCG)

%% Study of the Ca transients amplitude and the correlation with the presence of ripples
InfluenceCaTransientsAmpl(ExtractedResults,TTLstartsTimes,EvtStarts,ripples, nNeurons )

%% Plot of the time courses of neuronal activity in the SCEs with neurons ordered as their PFs. 
AnalyseSequences(MainPath,expPaths,DeepCAD)

%% Get coactivation matrix
coAct=GetCoactivationMat(putSCEtimes,sceNeuronIDs,nNeurons);

%% Get rank correlation between SCEs (work in progress)
minNbNeuronsSCE=10;                                                                                         % Can be use to restrict the number of SCEs to analyze : big ones involving more than 'minNbNeuronsSCE' neurons (set as nbStartsThr for all SCEs)
minCommNeurons=7;                                                                                           % Same for restricting the analysis to pair of SCEs that share more than 'minCommNeurons' neurons (set as 0 for all pairs)

bigSCEs=cellfun(@numel,sceNeuronIDs)>=minNbNeuronsSCE;                                                      % Restrict the SCEs to the 'big' ones (over 'minNbNeuronsSCE' neurons)
bigSCEidx=putSCEidx(bigSCEs); bigSCEnIDs=sceNeuronIDs(bigSCEs); bigSCEsync=syncSCEs(bigSCEs);               % Get corresponding data

MatCorCoef=zeros(numel(bigSCEnIDs));                                                                        % Initialization of matrices
MatCorP=ones(numel(bigSCEnIDs));
MatCorN=cell(numel(bigSCEnIDs));
for sce1=1:numel(bigSCEnIDs)-1                                                                              % Loop on SCEs
    for sce2=sce1+1:numel(bigSCEnIDs)                                                                           % Second loop on SCEs
        [MatCorN{sce1,sce2},rank1,rank2]=intersect(bigSCEnIDs{sce1}, bigSCEnIDs{sce2},'stable');                    % Get common neurons and their rank
        if numel(rank1)>=minCommNeurons                                                                             % If enough in common
            [MatCorCoef(sce1,sce2),MatCorP(sce1,sce2)]=corr(rank1,rank2,'Type','Spearman');                             % Get the rank correlation and the p-value
        end
    end
end
MatCorCoef(MatCorP>0.05)=0;                                                                                 % Discard p>0.05 correlations
    
    % Need to commpare with the number/fraction of correlated pairs when
    % the neurons IDs are shuffled or when their ranks are shuffled 
    
%% Analyse of the structure of sequences (based on bayesian reconstruction of position)

% Parameters
    nBins=200;                                                                                              % Number of spatial bins (200 --> 1bin = 1cm)
    wndw=0.02;                                                                                              % Time window used for decoding
    redetectEvents = 1; % redo the detection of events for behavior data
    
% Train Phase (on Behavior)
    % loade behavior data - we suppose that the analysis of activity has
    % already been done and we load the matrices Experiment and Neuron
    ExperimentB = load([cell2mat(expPaths.Behavior.Position) '\Experiment.mat']);    %ExperimentB = load('Experiment.mat');
    ExperimentB = ExperimentB.Experiment; 
    NeuronsB = load([cell2mat(expPaths.Behavior.Position) '\Neurons.mat']);    %NeuronsB = load('Neurons.mat');
    NeuronsB = NeuronsB.Neurons; 
    ParamsB = Params; % use the same parameters for behavior and sleep
    
    % optional redetection of events
    if redetectEvents == 1
        resultsB = load([cell2mat(expPaths.Behavior.Position) '\Results.mat']);
        resultsB = resultsB.results; 
        NeuronsB = GetRises(resultsB,ParamsB,ExperimentB,NeuronsB);
    else
        NeuronsB = Neurons;
    end

    % 1 ----  First trial : we use Dir1 and Dir2.  ---
    
    Dir1=[ExperimentB.Dir1.Starts ExperimentB.Dir1.Ends]; Dir2=[ExperimentB.Dir2.Starts ExperimentB.Dir2.Ends]; % Get Mouvement Periods
    Dir1 = Dir1(1:34, :); Dir2 = Dir2(1:33, :);                                                                % Get rid of the end of the behavior 
    spikesB=[];                                                                                              % Initialize 'spikes' (starts) sample mat. spikes = list of (t, neuron) couples
    EvtMatB=vertcat(NeuronsB.Rises.Matrix);EvtMatB=full(EvtMatB)>0;                                         % matrix of events. need to convert this a sample data; 
    for n=1:ExperimentB.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMatB(n,:)); 
        spikesB=[spikesB ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    
    [spikesMvt,~,~] = InIntervals(spikesB(:,1),sort([Dir1 ; Dir2]));                                         % Get Spikes occuring during mouvement periods
    spikesB=spikesB(spikesMvt,:);                                                                             % Restrict Spikes to mouvement periods

    positions=[(1:ExperimentB.SizeData-4)'...                                                                  % Normalize Position between [0 0.5]
        0.5*fillmissing((ExperimentB.positionXcmSmooth-min(ExperimentB.positionXcmSmooth))/(max(ExperimentB.positionXcmSmooth)-min(ExperimentB.positionXcmSmooth)),'previous')];
    [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);                                                         % Get positions during Dir2
    positions(Mvtd2,2)=1-positions(Mvtd2,2);                                                                % Linearisation (Dir1 = 0 --> 0.5 ; Dir2 = 0.5 --> 1)
    [Mvt,~,~] = InIntervals(positions(:,1),sort([Dir1 ; Dir2]));                                            % Get positions occuring during mouvement periods
    positions=positions(Mvt,:);                                                                             % Restrict positions to mouvement periods
    %figure; plot(positions(:,1), positions(:,2),'.');
    
    spikesB(:,1)=spikesB(:,1)/100;                                                                            % Timings (in ms) from indices
    positions(:,1)=positions(:,1)/100;                                                                      % Timings (in ms) from indices
    
    
%    Test the reconstruction on Awake itself (just to evaluate the performance)
    [statsTest,~,~] = ReconstructPosition(positions,spikesB,'mode','both','window',wndw,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
  %  estim = estim/nBins;
  %  figure; ax1 = subplot(2,1,1);  plot(estim, '.'); ax2 = subplot(2,1,2);  plot(statsTest.positions(:,2)/nBins, '.'); linkaxes([ax1,ax2],'x');
     figure; plot(estim, '.'); hold on; plot(pos, '.');
     ylabel('Mouse Position (cm)');
     xlabel('Time Points (wnd = 20ms)')
     set(gca, 'FontSize',14);
    error = min(abs(estim-pos),nBins-abs(estim-pos));
    errorCm=error*200/nBins;
    mean(errorCm)
    median(errorCm)
%     mean(errorCm(1:floor(end/2)))
%     median(errorCm(1:floor(end/2)))
    
   % train
    [~,lambda,Px] = ReconstructPosition(positions,spikesB,'mode','train','window',wndw,'nBins',[nBins 1]); % Train Reconstruction
 
    % test on same data
    [stats,~,~] = ReconstructPosition([],spikesB,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
    figure; imagesc(stats.estimations)
    xlabel('Time Points (window = 20ms)');
    ylabel('Recontructed position'); 
    set(gca, 'FontSize', 14);
    [~,estim]=max(stats.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    figure; 
    plot(stats.windows(:,1)+0.01, estim, '.'); hold on; plot(positions(:,1), positions(:,2)*200, '.'); 
    
    % Test Phase (on sleep)
 %   load('Px.mat'); load('Lambda.mat');                                                                     % Load Data from Train Phase                                                                      

    % test; On raccourcit la matrice des evts car �a prend beaucoup de
    % temps si je fais la d�tection des �vnts sans la condition sur la
    % d�riv�e. 
    EvtMat2 = EvtMat(:,1:60000);
    
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMat2(n,:)); 
        spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    
    spikes(:,1)=TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
    [stats,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction

    figure; imagesc(stats.estimations)
    xlabel('Time Points (window = 20ms)');
    ylabel('Recontructed position'); 
    set(gca, 'FontSize', 14);
   
    
    
    
       % 1 ----  Second trial : we separate directions  ---
    
    nBins = 100;
    
    Dir1=[ExperimentB.Dir1.Starts ExperimentB.Dir1.Ends]; Dir2=[ExperimentB.Dir2.Starts ExperimentB.Dir2.Ends]; % Get Mouvement Periods
    Dir1 = Dir1(1:34, :); Dir2 = Dir2(1:33, :);                                                                % Get rid of the end of the behavior 
    spikesB=[];                                                                                              % Initialize 'spikes' (starts) sample mat. spikes = list of (t, neuron) couples
    EvtMatB=vertcat(NeuronsB.Rises.Matrix);EvtMatB=full(EvtMatB)>0;                                         % matrix of events. need to convert this a sample data; 
    for n=1:ExperimentB.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMatB(n,:)); 
        spikesB=[spikesB ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    
    [spikesMvtDir1,~,~] = InIntervals(spikesB(:,1),Dir1);                                         % Get Spikes occuring during mouvement periods
    [spikesMvtDir2,~,~] = InIntervals(spikesB(:,1),Dir2);                                         % Get Spikes occuring during mouvement periods
 
    spikesBDir1=spikesB(spikesMvtDir1,:);                                                                             % Restrict Spikes to mouvement periods
    spikesBDir2=spikesB(spikesMvtDir2,:);                                                                             % Restrict Spikes to mouvement periods

    positions=[(1:ExperimentB.SizeData-4)'...                                                                  % Normalize Position between [0 0.5]
        fillmissing((ExperimentB.positionXcmSmooth-min(ExperimentB.positionXcmSmooth))/(max(ExperimentB.positionXcmSmooth)-min(ExperimentB.positionXcmSmooth)),'previous')];
     [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);                                                         % Get positions during Dir2
     [Mvtd1,~,~] = InIntervals(positions(:,1),Dir1);                                                         % Get positions during Dir2
     positionsDir1=positions(Mvtd1,:);                                                                             % Restrict positions to mouvement periods
     positionsDir2=positions(Mvtd2,:);                                                                             % Restrict positions to mouvement periods
%     figure; subplot(2,1,1); plot(positionsDir1(:,1), positionsDir1(:,2),'.'); subplot(2,1,2); plot(positionsDir2(:,1), positionsDir2(:,2),'.');
     
    spikesBDir1(:,1)=spikesBDir1(:,1)/100;                                                                            % Timings (in ms) from indices
    spikesBDir2(:,1)=spikesBDir2(:,1)/100;                                                                            % Timings (in ms) from indices
    positionsDir1(:,1)=positionsDir1(:,1)/100;                                                                      % Timings (in ms) from indices
    positionsDir2(:,1)=positionsDir2(:,1)/100;                                                                      % Timings (in ms) from indices
    
    
%    Test the reconstruction on Awake itself (just to evaluate the performance)

    %dir1
    [statsTest,~,~] = ReconstructPosition(positionsDir1,spikesBDir1,'mode','both','window',wndw,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
     figure; plot(estim, '.'); hold on; plot(pos, '.');
     ylabel('Mouse Position (cm)');
     xlabel('Time Points (wnd = 20ms)')
     set(gca, 'FontSize',14);
    error = abs(estim-pos);
    errorCm=error*100/nBins;
    mean(errorCm)
    median(errorCm)
   
    %dir2
    [statsTest,~,~] = ReconstructPosition(positionsDir2,spikesBDir2,'mode','both','window',wndw,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
     figure; plot(estim, '.'); hold on; plot(pos, '.');
     ylabel('Mouse Position (cm)');
     xlabel('Time Points (wnd = 20ms)')
     set(gca, 'FontSize',14);
    error = abs(estim-pos);
    errorCm=error*100/nBins;
    mean(errorCm)
    median(errorCm)
    
    
    % dir 1, train on full behavior data puis test on full behavior data
    
   % train
    [~,lambda,Px] = ReconstructPosition(positionsDir1,spikesBDir1,'mode','train','window',wndw,'nBins',[nBins 1]); % Train Reconstruction
    
    % order cells using lambda
     lambda1 = squeeze(lambda)';
     [Maxes1,Maxes1Idx]  = max(lambda1,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); lambda1=lambda1./Maxes1;
     figure;
     imagesc(lambda1(sMAxes1Idx,:));
    
    % test on behavior data
    [stats,~,~] = ReconstructPosition([],spikesBDir1,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
    figure; imagesc(stats.estimations)
    xlabel('Time Points (window = 20ms)');
    ylabel('Recontructed position'); 
    set(gca, 'FontSize', 14);
    [~,estim]=max(stats.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    figure; 
    plot(stats.windows(:,1)+0.01, estim, '.'); hold on; plot(positionsDir1(:,1), positionsDir1(:,2)*100, '.');     
      
    figure; % choix d'un systeme de temps; on prend celui donne pour l'instant � ReconstructPosition. Pour le sommeil on verra si on veut prendre le vrai temps pour plotter aussi l'ephy
    for kDir1 = 1 : size(Dir1,1)
        subplot(2,1,1)
        %plot(stats.windows(:,1)+0.01, estim, '.'); hold on; 
        plot(positions(Dir1(kDir1,1):Dir1(kDir1,2),1)/100, positions(Dir1(kDir1,1):Dir1(kDir1,2),2), '.');xlim([Dir1(kDir1,1)/100 Dir1(kDir1,2)/100]) ;
        subplot(2,1,2)
        imagesc(EvtMatB(sMAxes1Idx,Dir1(kDir1,1):Dir1(kDir1,2) ))
        pause;
    end
    
 
    
    % Test on sleep
    EvtMat2 = EvtMat;%(:,1:20000);                                                                            % On raccourcit la matrice des evts car �a prend beaucoup de temps si je fais la d�tection des �vnts sans la condition sur la d�riv�e. 
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMat2(n,:)); 
        spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    spikes(:,1)=TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
    [stats,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
    figure; imagesc(stats.estimations)
    xlabel('Time Points (window = 20ms)');
    ylabel('Recontructed position'); 
    set(gca, 'FontSize', 14);
    save('statsReconstructPos.mat','stats');
    
   
    % Then the estimation would be used to look at eventual trajectory replay during SCEs.

        % order cells using lambda
     %    
%     lambda1 = lambda(:,1:floor(end/2)); %Dir1
%     [Maxes1,Maxes1Idx]  = max(lambda1,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); lambda1=lambda1./Maxes1;
%     figure;
%     imagesc(lambda1(sMAxes1Idx,:));
%     
%     lambda2 = lambda(:, ceil(end/2):end); %Dir2
%     [Maxes2,Maxes2Idx]  = max(lambda2,[],2); [~,sMAxes2Idx]=sort(Maxes2Idx); lambda2=lambda2./Maxes2;
%     figure;
%     imagesc(lambda2(sMAxes2Idx,:));
    
    % find SCEs timings 
    nSyncSCEs = sum(syncSCEs); %nSyncSCEs=floor(nSyncSCEs/4); % we performed testing only on the first part of the trace to gain time
    SyncSCEsIdx = find(syncSCEs);
    nwinplot = 10;
    
    SumEvtStartsSmooth = sum(EvtStarts,1);
    SumEvtStartsSmooth2 = movmean(SumEvtStartsSmooth, 10) *10;
 %   figure; plot(SumEvtStartsSmooth); hold on; plot(SumEvtStartsSmooth2)

%     figure; 
%     for kSCE = 1 : size(putSCEtimes,1) %  nSyncSCEs %
%         kSCE
%             
%         %timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
%         timeSCE = putSCEtimes((kSCE)); %timing dans la ref ephy
%         [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
%      
%         [~,iSCE_wndw] = min(abs(stats.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
%          
%             %to check that the timing is OK I recalculate the estimated
%         %position on the considered time. 
%         
%            EvtMat3 = EvtMat(:,iSCE - nwinplot*(wndw/0.01)+2  : iSCE + nwinplot*(wndw/0.01)-1);                                                                            % On raccourcit la matrice des evts car �a prend beaucoup de temps si je fais la d�tection des �vnts sans la condition sur la d�riv�e. 
%             spikes3=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
%             for n=1:Experiment.nNeurons                                                                             % Loop on neurons
%                 spike3TimesN = find(EvtMat3(n,:)); 
%                 spikes3=[spikes3 ;  spike3TimesN' n*ones(numel(spike3TimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
%             end
%             spikes3(:,1)=TTLstartsTimes(spikes3(:,1));                                                                % Transform the starts idx in timings
%             [stats3,~,~] = ReconstructPosition([],spikes3,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
%         
%         
%         subplot(4,1,1);       
%         imagesc(stats.estimations(:, iSCE_wndw - nwinplot+1  : iSCE_wndw + nwinplot-1),[0 0.3])
%         xlabel('Time Points (window = 20ms)');
%         ylabel('Recontructed position'); 
%         set(gca, 'FontSize', 12);
%         title(['kSCE = ' num2str(kSCE,'%d')]); 
%           
%         subplot(4,1,3)
%  %       imagesc(EvtStarts(sMAxes1Idx,iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1))
%          imagesc(stats3.estimations,[0 0.3])
%         xlabel('Time Points (window = 20ms)');
%         ylabel('Recontructed position'); 
%         set(gca, 'FontSize', 12);
% 
%         subplot(4,1,4)
%         imagesc(EvtMat(sMAxes1Idx,iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1))
%         xlabel('Time Points (at 100Hz)');
%         ylabel('# neuron (ordered)'); 
%         set(gca, 'FontSize', 12);
%         
%         subplot(4,1,2)
%         plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1)); 
%         yline(nbStartsThr,'g');
%         xlabel('Time Points (at 100Hz)');
%         ylabel('SCE detection'); 
%         set(gca, 'FontSize', 12);
%        
%         pause; 
%     end

   figure; 
    for kSCE = 1 : size(putSCEtimes,1) %  nSyncSCEs %
        kSCE
        
        %timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
        timeSCE = putSCEtimes((kSCE)); %timing dans la ref ephy
        [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
        
        [~,iSCE_wndw] = min(abs(stats.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
              
        subplot(4,1,1);       
        imagesc(stats.estimations(:, iSCE_wndw - nwinplot+1  : iSCE_wndw + nwinplot-1),[0 0.3])
        xlabel('Time Points (window = 20ms)');
        ylabel('Recontructed position'); 
        set(gca, 'FontSize', 12);
        title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(4,1,3)
        imagesc(EvtStarts(sMAxes1Idx,iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1))
        xlabel('Time Points (at 100Hz)');
        ylabel('# neuron (ordered)'); 
        set(gca, 'FontSize', 12);

        subplot(4,1,4)
        imagesc(EvtMat(sMAxes1Idx,iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1))
        xlabel('Time Points (at 100Hz)');
        ylabel('# neuron (ordered)'); 
        set(gca, 'FontSize', 12);
        
        subplot(4,1,2)
        plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1)); 
        yline(nbStartsThr,'g');
        xlabel('Time Points (at 100Hz)');
        ylabel('SCE detection'); 
        set(gca, 'FontSize', 12);
       
        pause; 
    end  
    
    % now we want to use findReplayScore pour analyser les s�quences
    
    FindReplayScore()
    
    
    kSCE = 1;
    timeSCE = putSCEtimes((kSCE)); %timing dans la ref ephy
    [~,iSCE_wndw] = min(abs(stats.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
    figure; plot(stats.estimations(:,  iSCE_wndw ));
    
    
    
    %% Reconstruction (work in progress ; not working in that state (need to load awake then sleep...))
    % I keep here what was done by theo as a reference

if 0
% Parameters
    nBins=200;                                                                                              % Number of spatial bins (200 --> 1bin = 1cm)
    wndw=0.02;                                                                                              % Time window used for decoding
    
% Train Phase (on Behavior)
    Dir1=[Experiment.Dir1.Starts Experiment.Dir1.Ends]; Dir2=[Experiment.Dir2.Starts Experiment.Dir2.Ends]; % Get Mouvement Periods
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat. 
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
       spikes=[spikes ; Neurons.Rises(n).Starts' n*ones(numel(Neurons.Rises(n).Starts),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    
    [spikesMvt,~,~] = InIntervals(spikes(:,1),sort([Dir1 ; Dir2]));                                         % Get Spikes occuring during mouvement periods
    spikes=spikes(spikesMvt,:);                                                                             % Restrict Spikes to mouvement periods

    positions=[(1:Experiment.SizeData)'...                                                                  % Normalize Position between [0 0.5]
        0.5*fillmissing((Experiment.positionXcmSmooth-min(Experiment.positionXcmSmooth))/max(Experiment.positionXcmSmooth),'previous')];
    [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);                                                         % Get positions during Dir2
    positions(Mvtd2,2)=1-positions(Mvtd2,2);                                                                % Linearisation (Dir1 = 0 --> 0.5 ; Dir2 = 0.5 --> 1)

    [Mvt,~,~] = InIntervals(positions(:,1),sort([Dir1 ; Dir2]));                                            % Get positions occuring during mouvement periods
    positions=positions(Mvt,:);                                                                             % Restrict positions to mouvement periods
    
    spikes(:,1)=spikes(:,1)/100;                                                                            % Timings (in ms) from indices
    positions(:,1)=positions(:,1)/100;                                                                      % Timings (in ms) from indices
    
    [stats,lambda,Px] = ReconstructPosition(positions,spikes,'mode','train','window',wndw,'nBins',[nBins 1]); % Train Reconstruction
    
% %     If you want to test the reconstruction on Awake itself (just to evaluate the max performance)
%     [statsTest,~,~] = ReconstructPosition(positions,spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Test Reconstruction on itself
%     [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
%     toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
%     pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
%     error = min(abs(estim-pos),nBins-abs(estim-pos)); errorCm=error*200/nBins;                              % Define error as the difference between decoded and real position
% 
% %     figure
% %         imagesc(statsTest.estimations)
% %         hold on
% %         plot(stats.positions(:,2),'r','LineWidth',2)
% %         scatter(1:numel(stats.positions(:,2)),estim,'y','filled') 
    
% Test Phase (on sleep)
    load('Px.mat'); load('Lambda.mat');                                                                     % Load Data from Train Phase                                                                      
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
       spikes=[spikes ;Neurons.Rises(n).Starts' n*ones(numel(Neurons.Rises(n).Starts),1)];                      % Add neuron ID and starts idx to the spikes sample
    end
    spikes(:,1)=TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
    [stats,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
end

    % Then the estimation would be used to look at eventual trajectory replay during SCEs.

    

%%

% compare the different channels to find the best LFP when the ePhy
% recordings show SNR problems

vectSyncSWRs = zeros(1, 17); 
vectSyncSCEs = zeros(1, 17); 
vectnRipples = zeros(1, 17); 

for ePhysChan=1 : 17
    
    ePhysData=loadChanLFP(ePhysChan,[ePhysFile(1:end-3) 'lfp'],1);                                              % Load the LFP (create it from the .dat if it doesn't exist)

    % Compute Spectrogram and get power in ripple band
    rangeSpect = [0 500];                                                                                       % Range of the Spectrogram in Hz
    windowSpect = 0.05;                                                                                         % Window for computing the spectrogram in s
    [spectrogram,t,f] = MTSpectrogram(ePhysData,'range',[0 500],'window',0.05);                                 % Compute the spectrogram of the LFP
    bands = SpectrogramBands(spectrogram,f);                                                                    % Define physiological rhythms bands (default for ripples = [100 250] Hz)
    ripplesPower=[t bands.ripples];                                                                             % Get the power in the ripples band in each time window

    % PETD of ripple-band power around SCEs
    durations=[-1 1];                                                                                           % Window around events times to compute the PETD
    GetRipScePETD(ripplesPower,putSCEtimes,durations)                                                           % Get the Peri-Events(SCEs)-Time Distribution of the ripples power

    % Find Ripples
    RipBand=[100 250];                                                                                          % Define Ripples Band
    ePhysDataFilt=FilterLFP(ePhysData,'passband',RipBand);                                                      % Filter the LFP in the Ripples Band

    RipThr = [2 3]; % [2 3]; % [1 2]; % Default values in FindRipples function  [2 5]                                                                                          % Set Threshold for finding the ripples 
    [ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',RipThr);                                             % Find Ripples
    %[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',[0.2 0.4],'durations',[30 15 150]);                 % Weird things happen when finding ripples on G6 data...Thresholds are messed up, probably because of some artifacts so You'll need to get very low thr.
    %[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',[0.7 1.4],'durations',[30 20 100]);                 % Weird things happen when finding ripples on G6 data...Thresholds are messed up, probably because of some artifacts so You'll need to get very low thr.

    % plot ripples on lfp trace to check that the thesholds are OK
    figure; plot(ePhysData(:,1),ePhysData(:,2));
     hold on;
     for kripple = 1 : size(ripples, 1)
         %kripple = 1;
         xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
         yBox = [-20000, 20000, 20000, -20000, -20000];
         patch(xBox, yBox, 'green', 'FaceAlpha', 0.1, 'FaceColor', 'green', 'FaceAlpha', 0.1);
     end

       %plot ripples position on EvtMat matrix
    figure; [yPoints,xPoints] = find(EvtStarts==1); plot(TTLstartsTimes(xPoints),yPoints,'.k');
    hold on;
    for kripple = 1 : size(ripples, 1)
        %kripple = 1;
        xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
        yBox = [0, nNeurons, nNeurons, 0, 0];
        patch(xBox, yBox, 'green', 'FaceAlpha', 0.1, 'FaceColor', 'green', 'FaceAlpha', 0.1);
    end

       % histogram of ripple duration
    % figure; histogram(ripples(:,3)-ripples(:,1)); 

    % CCG of SCEs and SWRs
    binSize=0.01;  %0.01; 0.2                                                                                              % Bin size for Cross-CorreloGram
    durations=[-1 1];  %[-1 1]; [-20 20];                                                                                          % Window around events times to compute the CCG
    [RipSceCCG,tCCG]=GetRipSceCCG(ripples,putSCEtimes,durations,binSize);                                       % Get the CCG

    % Percentage of synchronous SCEs and SWRs
    MaxSep=0.1; %0.05; 0.1                                                                                                % Max time distance between the beginnings of two events to consider them synchronous 
    [syncSCEs,nbSyncSCEs]=GetNbSync(putSCEtimes,ripples(:,1),MaxSep);                                           % Get the number of SCEs, synced with SWRs
    [syncSWRs,nbSyncSWRs]=GetNbSync(ripples(:,1),putSCEtimes,MaxSep);                                           % Get the number of SWRs, synced with SCEs

    % Chance level of synchrony
    nbPerm=1000;                                                                                                % Number of shifted SCEs distributions to compute
    nbPermSyncSCEs=zeros(1,nbPerm);                                                                             % Initialization of matrices
    nbPermSyncSWRs=zeros(1,nbPerm);
    permCCG=zeros(nbPerm,size(RipSceCCG,1));
    for p=1:nbPerm                                                                                              % Loop
        permSCEidx=rem(putSCEidx+randi(numel(putSCE)),numel(putSCE))+1;                                             % Get a circularly shifted distribution of SCEs
        permSCEtimes=TTLstartsTimes(permSCEidx);                                                                                 % Get the corresponding timings
        [~,nbPermSyncSCEs(p)]=GetNbSync(permSCEtimes,ripples(:,1),MaxSep);                                          % Get the number of SCEs, synced with SWRs
        [~,nbPermSyncSWRs(p)]=GetNbSync(ripples(:,1),permSCEtimes,MaxSep);                                          % Get the number of SWRs, synced with SCEs
        [pccg,~]=GetRipSceCCG(ripples,permSCEtimes,durations,binSize);                                              % Compute CCG
        permCCG(p,:)=pccg(:,1,2);                                                                                   % Store the values
    end

    thrSCEs=prctile(nbPermSyncSCEs,95);                                                                         % Define thresholds for the chance nb of synchronies as the 95th percentile of the synchronies obtained whith shifted data
    thrSWRs=prctile(nbPermSyncSWRs,95);

    % Plot various figs
    PlotMisc(nbStarts,sceNeuronIDs,putSCEidx,syncSCEs,syncSWRs,ripples,nNeurons,tCCG,RipSceCCG,nbSyncSCEs,nbPermSyncSCEs,thrSCEs,nbSyncSWRs,nbPermSyncSWRs,thrSWRs,permCCG)

    vectSyncSWRs(ePhysChan) = nbSyncSWRs;
    vectSyncSCEs(ePhysChan) =  nbSyncSCEs;
    vectnRipples(ePhysChan) = size(ripples, 1);
        
    pause
    close all;
end





