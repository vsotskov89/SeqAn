%% Path of the Exp. Folder + Channels for TTL and HPC ePhy (Exp. Dep.) + Nb of sub Exps ([SleepPRE, Awake, SleepPOST])
root = '/media/sotskov/Vladimir/Analysis/Souris154148/22-07-07/';
ephys_path = '154148-sleeppost_2022-07-07_13-42-44/Record Node 104/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
timing_path = '134246_sCMOS_154148-sleeppost/sCMOS_AbsoluteTimingsSec.csv';
ePhysChan=6; TTLChan=18;

%% Load results and take only the last (i.e., sleeppost) session
load([root, 'AllExp/4-Results/Results.mat'])
timings = csvread([root, timing_path]);
n_conc = size(results.C, 2);
slen = length(timings) - 4; %since n_Downsample = 1, see ExtractExpData
toExtract = n_conc-slen+1:n_conc;

results.C_raw = results.C_raw(:,toExtract);
results.C = results.C(:,toExtract); 
results.S = results.S(:,toExtract);
results.AvgRoi = results.AvgRoi(:,toExtract); 
results.AvgRoiNoBaseline = results.AvgRoiNoBaseline(:,toExtract);

%% Load ePhysData (HPC channel and TTL)
ePhysData=loadChanLFP(ePhysChan,[root, ephys_path, 'continuous.lfp'],1);                                              % Load the LFP (create it from the .dat if it doesn't exist)
TTLData=loadChanDat(TTLChan,[root, ephys_path, 'continuous.dat']); TTLData=single(TTLData);                  % Load the TTL data
TTLstartsTimes=GetTTLtimes(TTLData,results);                                                       % Get the TTL times (which are also the time of Calcium data points)

%% Find ripples
rangeSpect = [0 500];                                                                                       % Range of the Spectrogram in Hz
windowSpect = 0.05;                                                                                         % Window for computing the spectrogram in s
[spectrogram,t,f] = MTSpectrogram(ePhysData,'range',[0 500],'window',0.05);                                 % Compute the spectrogram of the LFP
bands = SpectrogramBands(spectrogram,f);                                                                    % Define physiological rhythms bands (default for ripples = [100 250] Hz)
ripplesPower=[t bands.ripples];                                                                             % Get the power in the ripples band in each time window
% 
% % PETD of ripple-band power around SCEs
% durations=[-1 1];                                                                                           % Window around events times to compute the PETD
% GetRipScePETD(ripplesPower,putSCEtimes,durations)                                                           % Get the Peri-Events(SCEs)-Time Distribution of the ripples power

% Find Ripples
RipBand=[100 250];                                                                                          % Define Ripples Band
ePhysDataFilt=FilterLFP(ePhysData,'passband',RipBand);                                                      % Filter the LFP in the Ripples Band

RipThr = [2 3]; % [2 3]; % [1 2]; % Default values in FindRipples function  [2 5]                                                                                          % Set Threshold for finding the ripples 
[ripples,~,~] = FindRipples(ePhysDataFilt,'thresholds',RipThr);                                             % Find Ripples
%% Save results in readable form
csvwrite(['/export/home1/Sequences/220707/', 'SL_traces.csv'], [TTLstartsTimes, results.C_raw'])


%% main loop

files = dir([root, 'SL_events_sigma*.mat']);
for f = 1:length(files)
    load([root, files(f).name])
    %% Get calcium events starts
    % SizeData = size(results.C, 2); Experiment.SizeData = SizeData;                                              % Build false Experiment Struct.,                
    % nNeurons = size(results.C, 1); Experiment.nNeurons = nNeurons;
    % 
    % Params.Method.ActivityType =  'Binary';  %'Spikes';                                                                     % false Params Struct.,
    % Params.Method.SmallEvents = 'RisesDerivative';  % 'RisesSE'    'RisesSEbsd'  'RisesDecaysSEbsd'             % Choose the method to extract the rises
    % 
    % Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),...   % and false Neurons Struct. that are required for GetRises
    %     'Ends',cell(Experiment.nNeurons,1),'Duration',cell(Experiment.nNeurons,1)));
    %     
    % Neurons=GetRises(results,Params,Experiment,Neurons);                                                        % Get rises
    % EvtMat=vertcat(Neurons.Rises.Matrix);EvtMat=full(EvtMat)>0;                                                 % Make a binary matrix (ID,timings)
    % EvtStartsD=vertcat(Neurons.Rises.EvtPos);
    %% Get starts from already calculated events
    nNeurons = size(Evts, 2);
    nFrames = size(traces,1);
    EvtStarts = zeros(nNeurons, nFrames, 'logical');
    for n = 1:nNeurons
        nEvts = size(Evts{n},2);
        for ev = 1:nEvts
            st = find(TTLstartsTimes > Evts{n}{ev}.t0, 1);
            EvtStarts(n,st) = 1;
        end
    end
    
    
    %% Get Number of starts and neurons IDs
    wndw=10;                                                                                                    % Size of the moving window (in this version the window is [0 wndw] and not [-wdnw/2 wdnw/2])
    [nbStarts,NeuronIDs]=GetNbStarts(EvtStarts,wndw);                                                           % Get the number of starts in each window (here : from t=1 to t=SizeData-wndw)
    
    %% Get chance level of number of starts and set threshold for SCEs
    nbPerm=1000;                                                                                                % Number of random distribution of the number of start to compute (100 should be enough but 1000 or 10000 is good)
    AllNbStarts=zeros(nbPerm, nFrames -wndw);                                                                    % Initialize the nbStarts matrix (Perm #, timings)
    parfor p=1:nbPerm                                                                                              % Loop
        shifts=randi(nFrames,1,nNeurons) ;                                                                         % Get a random time shift for each neuron
        permEvt=cell2mat(arrayfun(@(x) circshift(EvtStarts(x,:),[1 shifts(x)]),(1:numel(shifts))','un',0));         % Apply a circular shift of said shifts to each neuron's event starts
        [AllNbStarts(p,:),~]=GetNbStarts(permEvt,wndw);                                                             % Get nbStarts for this time-shifted matrix of Evt Starts
    end
    nbStartsThr=prctile(AllNbStarts(:),99.9); %99  99.9                                                                  % Set the threshold for a SCE to the 99th percentile of these nbStarts
    
    %% Get SCEs
    [putSCE,putSCEidx,putSCEtimes]=GetSCEs(nbStarts,nbStartsThr,EvtStarts,wndw,TTLstartsTimes);                 % Find SCEs (moments where NbStarts exceed the threshold)
    sceNeuronIDs=NeuronIDs(putSCEidx);                                                                          % Get the neurons involved in each SCE
        
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
    
    thrSCEs=prctile(nbPermSyncSCEs,99);      %95                                                                   % Define thresholds for the chance nb of synchronies as the 95th percentile of the synchronies obtained whith shifted data
    thrSWRs=prctile(nbPermSyncSWRs,99);      %95 
    %% Plot CCGs
    figure; 
    bar(tCCG,RipSceCCG(:,1,2)); 
    xlabel('Delay (s)'); ylabel('Occurence (#)');title(['CCG of Ripples x SCEs , sigma = ', string(f+1)]) 
    hold on
    plot(tCCG,prctile(permCCG,95),'color','r','LineWidth',2)
    %% Plot various figs
    % PlotMisc(nbStarts,sceNeuronIDs,putSCEidx,syncSCEs,syncSWRs,ripples,nNeurons,tCCG,RipSceCCG,nbSyncSCEs,nbPermSyncSCEs,thrSCEs,nbSyncSWRs,nbPermSyncSWRs,thrSWRs,permCCG)
    %
    % %% Study of the Ca transients amplitude and the correlation with the presence of ripples
    % InfluenceCaTransientsAmpl(ExtractedResults,TTLstartsTimes,EvtStarts,ripples, nNeurons )
    % 
    % %% Plot of the time courses of neuronal activity in the SCEs with neurons ordered as their PFs. 
    % AnalyseSequences(MainPath,expPaths,DeepCAD)
    
    %% Get coactivation matrix
    %coAct=GetCoactivationMat(putSCEtimes,sceNeuronIDs,nNeurons);
    disp(['sigma = ', string(f+1)])
    disp(nnz(EvtStarts))
    disp(length(syncSCEs)) %total amount of SCEs
    disp(nbSyncSCEs)       %amount of syncchronized SCEs 
end    

