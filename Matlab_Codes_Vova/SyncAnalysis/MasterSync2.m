%% Path of the Exp. Folder + Channels for TTL and HPC ePhy (Exp. Dep.) + Nb of sub Exps ([SleepPRE, Awake, SleepPOST])

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
% this section is used if we want to analyse data at a slower frame rate
% create the tag to read results and scmos timing data
% the results file will be ['/Results' DeepCAD tagLowRate '.mat']
% the timing files will be ['/sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']
if nDownsample > 1
    tagLowRate = sprintf( '_ds%d',nDownsample); 
elseif nDownsample == 1 , tagLowRate = '';
else , return
end
if nDownsample > 1
    createLowRateTimingFiles(expPaths,nDownsample,tagLowRate);                     % create Low Rate timing files for the sCMOS camera in the same folder as original timing files, with the name ['/sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']
end

%% Load or extract the calcium results file
[ExtractedResults,  nameExtract]=LoadExtractedCalciumData(MainPath,expPaths,DeepCAD,nDownsample,tagLowRate);                                       % Load Calcium data (extract it from the concatenated file if it wasn't done before)

%% Load ePhysData (HPC channel and TTL)
[file,path] = uigetfile([MainPath '/*.dat'],'Select ePhys data to load'); ePhysFile = fullfile(path,file);  % Select the continuous.dat of the experiment    
ePhysData=loadChanLFP(ePhysChan,[ePhysFile(1:end-3) 'lfp'],1);                                              % Load the LFP (create it from the .dat if it doesn't exist)
TTLData=loadChanDat(TTLChan,ePhysFile); TTLData=single(TTLData);                                            % Load the TTL data

%% Extract TTL Timings
TTLstartsTimes=GetTTLtimes(TTLData,ExtractedResults);                                                       % Get the TTL times (which are also the time of Calcium data points)

%% Get calcium events starts
SizeData = size(ExtractedResults.C, 2); Experiment.SizeData = SizeData;                                     % Build false Experiment Struct.,                
nNeurons = size(ExtractedResults.C, 1); Experiment.nNeurons = nNeurons;

Params.Method.ActivityType =  'Binary';  %'Spikes';                                                                     % false Params Struct.,
Params.Method.SmallEvents = 'RisesDerivative';  % 'RisesSE'    'RisesSEbsd'  'RisesDecaysSEbsd'             % Choose the method to extract the rises

Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),...   % and false Neurons Struct. that are required for GetRises
    'Ends',cell(Experiment.nNeurons,1),'Duration',cell(Experiment.nNeurons,1)));
    
Neurons=GetRises(ExtractedResults,Params,Experiment,Neurons);                                               % Get rises
EvtMat=vertcat(Neurons.Rises.Matrix);EvtMat=full(EvtMat)>0;                                                 % Make a binary matrix (ID,timings)
nNeurons = 80;                                                                                            % if we want to restrict the number of neurons for this analysis
EvtMat = EvtMat(1:nNeurons, :);                                                                             % keep only the first neurons
neuronsPCDir2 = [2 5 11 14 19 20 22 23 25 26 30 33 35 37 38 40 42 43 54 56 57 61 66 69 71 72 74 78 87 89 90 94 101 111]; nNeurons = numel(neuronsPCDir2);  
%EvtStarts=[zeros(nNeurons,1) diff(EvtMat,[],2)]>0;    % change compared to theo (before : Experiment.nNeurons)  % Take only the rises starts
EvtStarts=vertcat(Neurons.Rises.EvtPos);
EvtStarts = EvtStarts(neuronsPCDir2, :);   %EvtStarts = EvtStarts(1:nNeurons, :);   
EvtLoc1 = vertcat(Neurons.Rises.EvtPos1);
EvtLoc1 = EvtLoc1(neuronsPCDir2, :); %EvtLoc1 = EvtLoc1(1:nNeurons, :);

% for behavior experiment, we need to select time windows when the mouse is still
% use nameExtract
%behav = ~isempty(cat(strfind(cell2mat(nameExtract),'Behav'), strfind(cell2mat(nameExtract),'awak'))) ;
% behav = ~isempty(strfind(cell2mat(nameExtract),'Behav')) ; % incomplet, a changer
% if behav 
%     load(fullfile(cell2mat(expPaths.Behavior.TimingsCMOS),'Experiment.mat'));
% %      figure; 
% %         ax1 = subplot(3,1,1);
% %         toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
% %         toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
% %         plot(Experiment.positionXcmSmooth); hold on; plot(toPlot1); plot(toPlot2)         
% %         ax2 = subplot(3,1,2);
% %         plot(Experiment.mouseSpeed); hold on; yline(2); yline(-2);
% %         ax3 = subplot(3,1,3);
% %         plot(Experiment.Still.Segments); ylim([-0.2 1.2]);
% %         linkaxes([ax1,ax2,ax3],'x');
%     % still segments : Experiment.Still.Segments
%     Experiment.Still.Segments = [false; false; false; false; Experiment.Still.Segments];
%     EvtStarts(:,not(Experiment.Still.Segments))=[];
%     TTLstartsTimes(not(Experiment.Still.Segments))=[];
%     SizeData = numel(TTLstartsTimes);
% end
% verifier que si on prend la solution load extracted data on a bien le bon
% nom dans nameextract. 


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
nbStartsThr=prctile(AllNbStarts(:),99.9); %99  99.9                                                                  % Set the threshold for a SCE to the 99th percentile of these nbStarts

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

% eliminer les ripples pendant la locomotion de la souris 
% if behav 
%     
% end
%plot ripples on lfp trace to check that the thesholds are OK
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

thrSCEs=prctile(nbPermSyncSCEs,99);      %95                                                                   % Define thresholds for the chance nb of synchronies as the 95th percentile of the synchronies obtained whith shifted data
thrSWRs=prctile(nbPermSyncSWRs,99);      %95 

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
    
    
 %% Analyse sequences like malvache et al

%parametres  
t0 = 20; t1 = 50; % en nb de frames 
 
% calcul d'un template comme celui utilise dans bsd. 
dt = 0.01; % frame rate en s
tauRise = 0.005; % rise time en s
tauDecay = 0.1; % decay time en s
delta = exp(-dt/tauRise - dt/tauDecay);
gamma = exp(-dt/tauRise)+exp(-dt/tauDecay)-exp(-dt/tauRise - dt/tauDecay) ;
eta = (tauRise/tauDecay)^(tauDecay/(tauDecay - tauRise))*( tauDecay/tauRise - 1)./(exp(-dt/tauDecay)-exp(-dt/tauRise));
s = zeros(t0+t1+1, 1);
s(t0+1) = 1;
c=filter(1,[eta, -eta * (gamma+delta)  ,  eta * delta  ],s);
%figure; plot(s); hold on; plot(c);

% calcul des positions des place fields. 
NeuronsB = load([cell2mat(expPaths.Behavior.Position) '\Neurons.mat']);   NeuronsB = NeuronsB.Neurons;

% il faut garder seulement certains neurones. 
neuronsPCDir2 = [2 5 11 14 19 20 22 23 25 26 30 33 35 37 38 40 42 43 54 56 57 61 66 69 71 72 74 78 87 89 90 94 101 111];

% SPAdir1=cell2mat({NeuronsB.SummedPlaceActivity.Dir1}');
% for neuron = 1:size(SPAdir1, 1), SPAdir1(neuron, :) = movmean(SPAdir1(neuron, :), 9); end
% [Maxes1,Maxes1Idx]  = max(SPAdir1,[],2); [~,sMaxes1Idx]=sort(Maxes1Idx); SPAdir1=SPAdir1./Maxes1;
% %figure; imagesc(SPAdir1(sMAxes1Idx,:));

SPAdir2=cell2mat({NeuronsB.SummedPlaceActivity.Dir2}');
for neuron = 1:size(SPAdir2, 1), SPAdir2(neuron, :) = movmean(SPAdir2(neuron, :), 9); end
[Maxes2,Maxes2Idx]  = max(SPAdir2,[],2); [~,sMaxes2Idx]=sort(Maxes2Idx); SPAdir2=SPAdir2./Maxes2;
%figure; imagesc(SPAdir2(sMaxes2Idx,:));

plotFig = true; %false; % 
nFig = 1; 
if plotFig == true
    for kFig = 1 : nFig
        figure(kFig);
    end
end
OnlySyncSCEs = false; %false : all SCEs 
if OnlySyncSCEs == true
    nSCEs = nbSyncSCEs;
else
    nSCEs = numel(putSCEidx);
end
slopes  = NaN(nSCEs,2);
vectRsq  = NaN(nSCEs,2); 
PearsonCorr = NaN(nSCEs,2);
PearsonPVal = NaN(nSCEs,2);
SpearmanCorr = NaN(nSCEs,2);
SpearmanPVal = NaN(nSCEs,2);


% kSCE = 2; %On va considerer toutes les SCEs synchrones dans un premier temps. Il y en a nbSyncSCEs. kSCE est un nombre entre 1 et nbSyncSCEs donnant le numero de la SCEs synchrone consideree.

%for kSCE = 1 : nSCEs 
for k = 1 : numel(replay4)   
kSCE = replay4(k);% %kSCE = replay2(k);
    
    % syncSCEs : vecteur dont la taille est le nb de SCEs et dont les éléments sont égaus 1 si la SCE correspondante est synchrone et 0 sinon. 
    % nbStarts : vecteur de longueur le nb d'images acquise - wndw. Les elements correspondent au nombre de départ d'evenements dans une fenetre wndw (100ms en general)
    % NeuronIDs : cell array de la meme taille que nbStarts et qui donne les IDs des neurones qui ont un EvtStart dans une fenetre de 100ms. 
    % putSCEtimes : vecteur dont le nombre d'elements est le nb de SCES et qui donne timing dans la ref ephy de chaque SCE
    % putSCEidx :  vecteur dont le nombre d'elements est le nb de SCES et qui donne le #frame dans l'acquisition optique de chaque SCE
    % putSCE :  vecteur de longueur le nb d'images acquises dans l'expériences et dont les éléments égaux à 1 correspondent aux positions des SCEs, les autres éléments sont égaux à 0. 
    % TTLStartTimes: vecteur de longueur le nb d'images acquises dans l'expériences et qui donne pour chaque image le timing dans le referentiel d'ephy
    % sceNeuronIDs : IDs of neurons involved in SCEs

    if OnlySyncSCEs == true
         SyncSCEsIdx = find(syncSCEs); % vecteur donnant les # des SCEs synchrones
         k2SCE = SyncSCEsIdx(kSCE);  % numero de la sce parmi toute les SCEs.
    else
         k2SCE = kSCE;
    end

    timeSCEfr = putSCEidx(k2SCE); %timing de la SCEs en # de frames optique
    neuronsSCE = sceNeuronIDs{k2SCE}; %neurones qui participent à la SCE. 
    tmin =  timeSCEfr - t0; tmax =  timeSCEfr + t1; % bornes d'extraction des traces calciques

    % methode 1 : covariance avec un template
    
    % choix d'un template
%    medianTrace = median( ExtractedResults.C_raw(neuronsSCE, tmin : tmax), 1 ); % calcul de la trace calcique mediane du SCE
%    template = medianTrace;
%    template = c;

    %plot des traces calciques de tous les neurones et du template
%     figure(1);
%     plot(template);
%     hold on; 
%     for kneuron = 1 : numel(neuronsSCE)
%         plot(kneuron * 20 + ExtractedResults.C_raw(neuronsSCE(kneuron), tmin : tmax));
%     end
%     xline(t0+1);
%  %   pause; 
%     hold off; 
     

%     % computation of the normalised covariance between each transient and
%     % the median trace and determination of the position of the max. we
%     % keep this as the position of the event if the max value of the
%     % covariance is larger than 0.6. Later we will find the position of the
%     % max using a parabolic fit. 
%     covmat = zeros(numel(neuronsSCE), 2*(t0+t1)+1);
%     evtPosMat = false(nNeurons, 2*(t0+t1)+1); % matrice binaire representant les départs d'activite de tous les neurones pendant la fenetre concernee.  de taille tous les neurones et 
%     evtPosVect = nan(numel(neuronsSCE),1);
%     for kneuron = 1 : numel(neuronsSCE)
%         covmat(kneuron, :) =  xcov(ExtractedResults.C_raw(neuronsSCE(kneuron), tmin : tmax)-ExtractedResults.C_raw(neuronsSCE(kneuron), tmin ),template,'normalized');
%         [maxCov,xmax] = max(covmat(kneuron, :));
%         if maxCov >0.6
%             evtPosMat(neuronsSCE(kneuron),xmax) = true; 
%             evtPosVect(kneuron) = xmax;
%         end
%     end
%     % figure; imagesc(evtPos);
%     
%     %figure; % plot des images
% %     subplot(1,2,1); imagesc(evtPosMat(sMaxes1Idx',60:80)); hold off; 
% %     subplot(1,2,2); imagesc(evtPosMat(sMaxes2Idx',60:80)); hold off; pause; 
%     
%      %figure; % plot des nuages de points;
%     figure(2);
%     subplot(2,2,1); scatter(Maxes1Idx(neuronsSCE),evtPosVect); hold off; 
%     subplot(2,2,3); scatter(Maxes2Idx(neuronsSCE),evtPosVect); hold off;
%     subplot(2,2,2:2:4); 
%         for kneuron = 1 : numel(neuronsSCE)
%            if  ~isnan(evtPosVect(kneuron)) 
%             plot(covmat(kneuron, :)); hold on;
%            end
%         end
%         hold off;
%     pause;
%      
    
  % methode 2 : calcul de la position de l'évènement comme le fit du max de la derivee de la trace calcique 
    evtPosVect = zeros(numel(neuronsSCE),1);
    evtPosVectTot = NaN(nNeurons,1);

    for kneuron = 1 : numel(neuronsSCE)
           dataC = ExtractedResults.C_raw(neuronsSCE(kneuron), tmin : tmax);
           gdatac = gradient(dataC);   %gdatac(2:end) = gdatac(1:end-1); 
           %[~,xmax2] = max(gdatac(20:31)); xmax2 = xmax2+19;
           %evtPosVect(kneuron) = xmax2;
            xmin = 20; xmax = 30; 
            [~,xmax2] = max(gdatac(xmin:xmax)); 
            % fit de la derivee par une gaussienne autour du max
            xminFit = xmax2+xmin-4; xmaxFit =xmax2+xmin+3;
            x = xminFit:1:xmaxFit;
            p =gaussfitn(x',gdatac(xminFit:xmaxFit)', [], {-Inf, 0, xminFit},{+Inf, +Inf, xmaxFit}  ); %, {-100, 0, xminFit},{+100, 200, xmaxFit} ); 
            [~,~,xmax2,~]=p{:}; 
            evtPosVect(kneuron) = xmax2;
            evtPosVectTot(neuronsSCE(kneuron)) = xmax2;
    end
       
   % figure(2); % plot des nuages de points;
   % plot pour la dir 2. On va comparer les résultats si on prend les plus
   % belles cellules ou si on prend toutes les cellules. 
   
        % we first look at only the nicest PCs in dir 2
        neuronsPCsDir2SCE = intersect(neuronsSCE,neuronsPCDir2);  % l'intersection entre les vecteurs neuronsSCE et neuronsPCDir2. 
        if numel(neuronsPCsDir2SCE)>5
           x = Maxes2Idx(neuronsPCsDir2SCE);
           y = evtPosVectTot(neuronsPCsDir2SCE);
           P1 = polyfit(x,y,1); % [P1 ErrorStruct]
           yfit = P1(1)*x+P1(2); %y_fit = polyval(P1,x); %[y_fit,delta] = polyval(P1,x,ErrorStruct);
           slopes(kSCE,1) = P1(1); 
           SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
           SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
           vectRsq(kSCE,1) = 1-SSres/SStot;                                    % R^2
           [SpearmanCorr(kSCE,1),SpearmanPVal(kSCE,1)] = corr(x,y,'Type','Spearman');
           [PearsonCorr(kSCE,1),PearsonPVal(kSCE,1)] = corr(x,y,'Type','Pearson');
           PearsonCorr(kSCE,1)
           PearsonPVal(kSCE,1)
           SpearmanCorr(kSCE,1)
           SpearmanPVal(kSCE,1)
           if plotFig == true
               subplot(2,1,1); 
               scatter(x,y); hold on;
               plot(x,yfit,'r-');
               %plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--');
               %legend('Data','Linear Fit');
               legend('Data',sprintf("Linear fit, R2 = %f",vectRsq(kSCE,1)));
               %text(70, 30, sprintf("R2 = %f",vectRsq(kSCE,1)))
               hold off;
           end
             
        end
     % pour le 2ème subplot on prend toutes les cellules comme avant pour comparer. On reste sur la dir2. 
        x = Maxes2Idx(neuronsSCE);
        y = evtPosVect;
        P1 = polyfit(x,y,1); % [P1 ErrorStruct]
        yfit = P1(1)*x+P1(2); %y_fit = polyval(P1,x); [y_fit,delta] = polyval(P1,x,ErrorStruct);
        slopes(kSCE,2) = P1(1);
        SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
        SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
        vectRsq(kSCE,2) = 1-SSres/SStot;                        % R^2
        [SpearmanCorr(kSCE,2),SpearmanPVal(kSCE,2)] = corr(x,y,'Type','Spearman');
        [PearsonCorr(kSCE,2),PearsonPVal(kSCE,2)] = corr(x,y,'Type','Pearson');
        if plotFig == true
            subplot(2,1,2);
            scatter(x,y);  hold on;
            plot(x,yfit,'r-'); 
            %legend('Data','Linear Fit');
            legend('Data',sprintf("Linear fit, R2 = %f",vectRsq(kSCE,2)));
            hold off; pause;
            %plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
        end
end

figure; 
subplot(2,1,1); histogram(slopes(:,1), 100); 
subplot(2,1,2); histogram(slopes(:,2), 100); 
replay1 = find(slopes(:,1)<-0.05);
replay2 = find(slopes(:,1)>+0.05);

figure; 
subplot(2,1,1); histogram(vectRsq(:,1), 100); 
subplot(2,1,2); histogram(vectRsq(:,2), 100); 
replay3 = find(vectRsq(:,1)>0.5);

figure; 
subplot(2,1,1); histogram(PearsonPVal(:,1), 100); 
subplot(2,1,2); histogram(PearsonPVal(:,2), 100); 

figure; 
subplot(2,1,1); histogram(SpearmanPVal(:,1), 100); 
subplot(2,1,2); histogram(SpearmanPVal(:,2), 100); 

replay4 = find(PearsonPVal(:,1)<+0.05);
replay5 = find(SpearmanPVal(:,1)<+0.05);

figure; scatter(PearsonCorr(:,1),PearsonPVal(:,1));

figure; histogram(PearsonCorr(replay4,1), 100); 
sum(PearsonCorr(replay4,1)<0)
sum(PearsonCorr(replay4,1)>0)


% syncSCEs : vecteur dont la taille est le nb de SCEs et dont les éléments sont égaux 1 si la SCE correspondante est synchrone et 0 sinon. 
syncSCEsList = find(syncSCEs>0);
replay4SyncSWRs = intersect(replay4,syncSCEsList); % SCEs correlees avec les SWRs et qui ont un pval <0.05 pour la correlation de pearson
ConsideredSCEs = find(~isnan(PearsonPVal(:,1)));
ConsideredSCEsSyncSWRs =  intersect(ConsideredSCEs,syncSCEsList);
numel(replay4SyncSWRs)
numel(ConsideredSCEsSyncSWRs)



% aides diverses aux calculs faits juste au dessus

figure;
for kneuron = 1 : numel(neuronsSCE)
    kneuron
%    [~,xmax] = max(covmat(kneuron, :));
    dataC = ExtractedResults.C_raw(neuronsSCE(kneuron), tmin : tmax);
    gdatac = gradient(dataC);   %gdatac(2:end) = gdatac(1:end-1); 
    xmin = 20; xmax = 30; 
    [~,xmax2] = max(gdatac(xmin:xmax)); 
    % fit de la derivee par une courbe parabolique autour du max
    xminFit = xmax2+xmin-3; xmaxFit =xmax2+xmin+2;
    x = xminFit:1:xmaxFit;

    % first attempt of fitting : function fit with 'gauss1'. but gauss1
    % does not incude an offset...
%     options = fitoptions('gauss1');
%     options.Lower = [0 -Inf 0];
%     [f, gof] = fit(x',gdatac(xminFit:xmaxFit)','gauss1', options);
%     gof
%     xplot = xminFit:0.1:xmaxFit;
%     fittedCurve =  (f.a1)*exp(-((xplot-(f.b1))/(f.c1)).^2); 
%     xmax2 = f.b1;
    
    % second attempt : function gaussfitn
    [p,resnorm, residual] =gaussfitn(x',gdatac(xminFit:xmaxFit)',  [], {-Inf, 0, xminFit},{+Inf, +Inf, xmaxFit}  ); %, {-100, 0, xminFit},{+100, 200, xmaxFit} ); 
   % resnorm 
    [v0,amp,mu,var]=p{:}; sig=sqrt(var);
    xplot = xminFit:0.1:xmaxFit;
    fittedCurve =  amp*exp(-(((xplot-mu).^2)/(2*sig.^2)))+ v0;
    xmax2 = mu;
    
    template = c;
%   subplot(2,1,1); 
    title(sprintf('%d',kneuron)); plot(template*20); hold on; plot(20 + dataC); plot(gdatac*3); plot(xplot, fittedCurve*3);  xline(xmax2, 'g'); %xline(xmax-50);
    hold off; 
%    subplot(2,1,2); plot(covmat(kneuron, :)); xline(xmax);
    pause; 
end

figure; plot( xcov( c,c,'normalized'));

% tentatives de calcul de trace de reference 
% utiliser les spikes trouvés par BSD ? On prend les spikes isolés et on
% mesure un spike triger average de la trace calcique. 
kneuron = 4; 
figure; plot(results.C_raw(kneuron,:)); hold on; plot(results.C(kneuron,:)); plot(results.S(kneuron,:));

    %fit to find the position and value of the maximum 
    
mu=5;
sig=2.3;
amp=4;
v0=1.7;
x=linspace(0,10,100);
y=amp*exp(-(((x-mu).^2)/(2*sig.^2)))+ v0;
 y=y+randn(size(y))*0.1;
 
p=gaussfitn(x(:),y(:));  
    
    
%end

% kneuron = 2;
% figure; plot(covmat(kneuron, :))
% x = linspace(-2*t0,1,2*t1);
% [maxCov,xmax] = max(covmat(kneuron, :));
% 
% p = polyfit(x,y,4);

%     % plot pour toutes les cellules
%     % figure(2); % plot des nuages de points;
%    % subplot(2,1,1); 
%         x = Maxes1Idx(neuronsSCE);
%         y = evtPosVect;
%    %     scatter(x,y);
%         P1 = polyfit(x,y,1); % [P1 ErrorStruct]
%    %     yfit = P1(1)*x+P1(2); %y_fit = polyval(P1,x);
%         %[y_fit,delta] = polyval(P1,x,ErrorStruct);
%         slopes(kSCE,1) = P1(1); 
%    %     hold on;
%    %     plot(x,yfit,'r-')
%         %plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
%    %     SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
%    %     SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
%    %     Rsq = 1-SSres/SStot                                    % R^2
%    %     legend('Data','Linear Fit');
%    %     hold off; 
%    % subplot(2,1,2); 
%         x = Maxes2Idx(neuronsSCE);
%         y = evtPosVect;
%    %     scatter(x,y);
%         P1 = polyfit(x,y,1); % [P1 ErrorStruct]
%    %     yfit = P1(1)*x+P1(2); %y_fit = polyval(P1,x);
%         %[y_fit,delta] = polyval(P1,x,ErrorStruct);
%         slopes(kSCE,2) = P1(1);
% %         hold on;
% %         plot(x,yfit,'r-')
% %         %plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
% %         SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
% %         SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
% %         Rsq = 1-SSres/SStot                                    % R^2
% %         legend('Data','Linear Fit');
% %         hold off; 
% %     pause;
%     


%% Analyse of the structure of sequences (based on bayesian reconstruction of position)

    
% Load/compute behavior data
    redetectEvents = 0; % redo the detection of events for behavior data
    ExperimentB = load([cell2mat(expPaths.Behavior.Position) '\Experiment.mat']); ExperimentB = ExperimentB.Experiment;   % load behavior data - we suppose that the analysis of activity has already been done and we load the matrices Experiment and Neuron
    NeuronsB = load([cell2mat(expPaths.Behavior.Position) '\Neurons.mat']);   NeuronsB = NeuronsB.Neurons; 
    %ParamsB = Params;                                                                                           % use the same parameters for behavior and sleep
    ParamsB.Method.ActivityType =  'Binary'; %  'Binary'; %                                                                    % false Params Struct.,
    ParamsB.Method.SmallEvents = 'RisesDerivative';% 'Sbsd';   % 'RisesSE'    'RisesSEbsd'  'RisesDecaysSEbsd'             % Choose the method to extract the rises
  
    if redetectEvents == 1                                                                                     % optional redetection of events during behavior
        resultsB = load([cell2mat(expPaths.Behavior.Position) '\Results.mat']); resultsB = resultsB.results;    % load results
        NeuronsB = GetRises(resultsB,ParamsB,ExperimentB,NeuronsB);                                             % compute activity using GetRises
    %else , NeuronsB = Neurons; % je ne comprends pas pourquoi j'ai mis cette ligne
    end


     
    % 1 ----  First trial : we use Dir1 and Dir2.  ---

%     nBins=200;                                                                                              % Number of spatial bins (200 --> 1bin = 1cm)
%     wndw=0.02;                                                                                              % Time window used for decoding
%     Dir1=[ExperimentB.Dir1.Starts ExperimentB.Dir1.Ends]; Dir2=[ExperimentB.Dir2.Starts ExperimentB.Dir2.Ends]; % Get Mouvement Periods
%     Dir1 = Dir1(1:34, :); Dir2 = Dir2(1:33, :);                                                                % Get rid of the end of the behavior (on mouse 154148, session 30/06/22 the end of behavior is bad quality)
%     spikesB=[];                                                                                              % Initialize 'spikes' (starts) sample mat. spikes = list of (t, neuron) couples
%     EvtMatB=vertcat(NeuronsB.Rises.Matrix);EvtMatB=full(EvtMatB)>0;                                         % matrix of events. need to convert this a sample data; 
%     for n=1:ExperimentB.nNeurons                                                                             % Loop on neurons
%         spikeTimesN = find(EvtMatB(n,:)); 
%         spikesB=[spikesB ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
%     end
%     
%     [spikesMvt,~,~] = InIntervals(spikesB(:,1),sort([Dir1 ; Dir2]));                                         % Get Spikes occuring during mouvement periods
%     spikesB=spikesB(spikesMvt,:);                                                                             % Restrict Spikes to mouvement periods
% 
%     positions=[(1:ExperimentB.SizeData-4)'...                                                                  % Normalize Position between [0 0.5]
%         0.5*fillmissing((ExperimentB.positionXcmSmooth-min(ExperimentB.positionXcmSmooth))/(max(ExperimentB.positionXcmSmooth)-min(ExperimentB.positionXcmSmooth)),'previous')];
%     [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);                                                         % Get positions during Dir2
%     positions(Mvtd2,2)=1-positions(Mvtd2,2);                                                                % Linearisation (Dir1 = 0 --> 0.5 ; Dir2 = 0.5 --> 1)
%     [Mvt,~,~] = InIntervals(positions(:,1),sort([Dir1 ; Dir2]));                                            % Get positions occuring during mouvement periods
%     positions=positions(Mvt,:);                                                                             % Restrict positions to mouvement periods
%     %figure; plot(positions(:,1), positions(:,2),'.');
%     
%     spikesB(:,1)=spikesB(:,1)/100;                                                                            % Timings (in ms) from indices
%     positions(:,1)=positions(:,1)/100;                                                                      % Timings (in ms) from indices
%         
% %    Test the reconstruction on Awake itself (just to evaluate the performance)
%     [statsTest,~,~] = ReconstructPosition(positions,spikesB,'mode','both','window',wndw,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
%     [~,estim]=max(statsTest.estimations);                                                                   % Define Estimated Position as the bin with the max probability
%     toKeep=sum(~ismember(statsTest.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
%     pos=statsTest.positions(toKeep,2); estim=estim(toKeep)';                                                % Restrict positions and estimations to mouvement periods
%   %  figure; ax1 = subplot(2,1,1);  plot(estim, '.'); ax2 = subplot(2,1,2);  plot(statsTest.positions(:,2)/nBins, '.'); linkaxes([ax1,ax2],'x');
%     figure; plot(estim, '.'); hold on; plot(pos, '.'); ylabel('Mouse Position (cm)'); xlabel('Time Points (wnd = 20ms)');  set(gca, 'FontSize',14);
%     error = min(abs(estim-pos),nBins-abs(estim-pos));  errorCm=error*200/nBins;
%     mean(errorCm) ,    median(errorCm)
%     
%    % train on all behavior data
%     [~,lambda,Px] = ReconstructPosition(positions,spikesB,'mode','train','window',wndw,'nBins',[nBins 1]); % Train Reconstruction
%  
%     % test on same data
%     [stats,~,~] = ReconstructPosition([],spikesB,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
%     figure; imagesc(stats.estimations); xlabel('Time Points (window = 20ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
%     [~,estim]=max(stats.estimations);                                                                   % Define Estimated Position as the bin with the max probability
%     figure; plot(stats.windows(:,1)+0.01, estim, '.'); hold on; plot(positions(:,1), positions(:,2)*200, '.'); 
%     
%     % Test on sleep
%  %   load('Px.mat'); load('Lambda.mat');                                                                     % Load Data from Train Phase                                                                      
%     EvtMat2 = EvtMat(:,1:60000);                                                                            % I make the matrix shorter for the program to run faster
%     spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
%     for n=1:Experiment.nNeurons                                                                             % Loop on neurons
%         spikeTimesN = find(EvtMat2(n,:)); 
%         spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
%     end
%     spikes(:,1)=TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
%     [stats,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px); % Try Reconstruction
%     figure; imagesc(stats.estimations); xlabel('Time Points (window = 20ms)'); ylabel('Recontructed position');  set(gca, 'FontSize', 14);
%            
%     % order cells using lambda    
%     lambda1 = lambda(:,1:floor(end/2)); %Dir1
%     [Maxes1,Maxes1Idx]  = max(lambda1,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); lambda1=lambda1./Maxes1;
%     figure; imagesc(lambda1(sMAxes1Idx,:));
%     lambda2 = lambda(:, ceil(end/2):end); %Dir2
%     [Maxes2,Maxes2Idx]  = max(lambda2,[],2); [~,sMAxes2Idx]=sort(Maxes2Idx); lambda2=lambda2./Maxes2;
%     figure; imagesc(lambda2(sMAxes2Idx,:));
%     
       % 1 ----  Second trial : we separate directions  ---
  
    % finish computing behavior data in order to train the algorithm    
    Dir1=[ExperimentB.Dir1.Starts ExperimentB.Dir1.Ends]; Dir2=[ExperimentB.Dir2.Starts ExperimentB.Dir2.Ends]; % Get Mouvement Periods
    Dir1 = Dir1(1:34, :); Dir2 = Dir2(1:33, :);                                                                % Get rid of the end of the behavior 
    spikesB=[];                                                                                              % Initialize 'spikes' (starts) sample mat. spikes = list of (t, neuron) couples
    EvtMatB=vertcat(NeuronsB.Rises.Matrix);EvtMatB=full(EvtMatB)>0;                                         % matrix of events. need to convert this a sample data; 
    for n=1:ExperimentB.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMatB(n,:)); 
        spikesB = [spikesB ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];                     % Add neuron ID and starts idx to the spikes sample
    end
    [spikesMvtDir1,~,~] = InIntervals(spikesB(:,1),Dir1); [spikesMvtDir2,~,~] = InIntervals(spikesB(:,1),Dir2);     % Get Spikes occuring during mouvement periods
    spikesBDir1 = spikesB(spikesMvtDir1,:); spikesBDir2=spikesB(spikesMvtDir2,:);                                     % Restrict Spikes to mouvement periods
    spikesBDir1 = sortrows(spikesBDir1,1); spikesBDir2 = sortrows(spikesBDir2,1);
    
     positions=[(1:ExperimentB.SizeData-4)'...                                                                  % Normalize Position between [0 1]
     fillmissing((ExperimentB.positionXcmSmooth-min(ExperimentB.positionXcmSmooth))/(max(ExperimentB.positionXcmSmooth)-min(ExperimentB.positionXcmSmooth)),'previous')];
     [Mvtd1,~,~] = InIntervals(positions(:,1),Dir1);  [Mvtd2,~,~] = InIntervals(positions(:,1),Dir2);          % Get positions during Dir1 and Dir2
     positionsDir1=positions(Mvtd1,:); positionsDir2=positions(Mvtd2,:);                                       % Restrict positions to mouvement periods
%     figure; subplot(2,1,1); plot(positionsDir1(:,1), positionsDir1(:,2),'.'); subplot(2,1,2); plot(positionsDir2(:,1), positionsDir2(:,2),'.');

    spikesBDir1(:,1)=spikesBDir1(:,1)/100; spikesBDir2(:,1)=spikesBDir2(:,1)/100;                               % Timings (in ms) from indices
    positionsDir1(:,1)=positionsDir1(:,1)/100; positionsDir2(:,1)=positionsDir2(:,1)/100;                       % Timings (in ms) from indices
    
    
    % analyse EvtMatB
    EvtMatBDir1 = EvtMatB(:,logical(ExperimentB.Dir1.Segments));
    %figure; imagesc(EvtMatBDir1);
     SumEvtMatBDir1 = sum(EvtMatBDir1,1); mean(SumEvtMatBDir1) 
     SumEvtMatB = sum(EvtMatB,1); mean(SumEvtMatB) ;
    

%    First test : test the reconstruction on awake itself (just to evaluate the performance) Train on 60% of the data and test on 40%
    wndw2=0.01;                                                                                              % Time window used for decoding
    nBins=100;                                                                                              % Number of spatial bins (200 --> 1bin = 1cm)

    %dir1
    [statsTest1,~,~] = ReconstructPosition(positionsDir1,spikesBDir1,'mode','both','window',wndw2,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estimTest1]=max(statsTest1.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest1.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest1.positions(toKeep,2); estimTest1=estimTest1(toKeep)';                                                % Restrict positions and estimations to mouvement periods
    figure; plot(estimTest1, '.'); hold on; plot(pos, '.'); ylabel('Mouse Position (cm)'); xlabel('Time Points (wnd = 10ms)'); set(gca, 'FontSize',14);
    error = abs(estimTest1-pos); errorCm=error*100/nBins;  mean(errorCm), median(errorCm)
   
    [~,estimTest1]=max(statsTest1.estimations); 
    figure; 
    ax1 = subplot(2,1,1); plot(statsTest1.windows(:,1), estimTest1, '.'); hold on; plot(statsTest1.windows(:,1), statsTest1.positions(:,2), '.'); plot(positionsDir1(:,1), positionsDir1(:,2)*100,'.'); ylabel('Mouse Position (cm)'); xlabel('Time Points (wnd = 10ms)'); set(gca, 'FontSize',14);
    ax2 = subplot(2,1,2); plot(positions(:,1)/100,SumEvtMatB(1,5:end),'.');
    linkaxes([ax1,ax2],'x');
    
    %dir2
    [statsTest2,~,~] = ReconstructPosition(positionsDir2,spikesBDir2,'mode','both','window',wndw2,'nBins',[nBins 1],'training', 0.6);  % Test Reconstruction on itself
    [~,estimTest2]=max(statsTest2.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    toKeep=sum(~ismember(statsTest2.estimations,ones(nBins,1)/nBins))>0;                                     % Find Mouvement period indices
    pos=statsTest2.positions(toKeep,2); estimTest2=estimTest2(toKeep)';                                                % Restrict positions and estimations to mouvement periods
    figure; plot(estimTest2, '.'); hold on; plot(pos, '.'); ylabel('Mouse Position (cm)'); xlabel('Time Points (wnd = 10ms)'); set(gca, 'FontSize',14);
    error = abs(estimTest2-pos); errorCm=error*100/nBins; mean(errorCm), median(errorCm)
    
    
%  now train on 100% of behavior data and test on sleep
    
    %dir 1
    
    [~,lambda1,Px1] = ReconstructPosition(positionsDir1,spikesBDir1,'mode','train','window',wndw2,'nBins',[nBins 1]); % Train Reconstruction
    save('lambda_10ms.mat', 'lambda1'); save('Px_10ms.mat', 'Px1'); 
    
    lambda1s = squeeze(lambda1)'; [Maxes1,Maxes1Idx]  = max(lambda1s,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); lambda1s=lambda1s./Maxes1; % order cells using lambda
    figure; imagesc(lambda1s(sMAxes1Idx,:));
    
    [stats1b,~,~] = ReconstructPosition([],spikesBDir1,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1,'Px',Px1); % Try Reconstruction on behvior data first
    figure; imagesc(stats1b.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    [~,estim1]=max(stats1b.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    figure; plot(stats1b.windows(:,1)+0.01, estim1, '.'); hold on; plot(positionsDir1(:,1), positionsDir1(:,2)*100, '.');     
%    figure; % choix d'un systeme de temps; on prend celui donne pour l'instant à ReconstructPosition. Pour le sommeil on verra si on veut prendre le vrai temps pour plotter aussi l'ephy
%     for kDir1 = 1 : size(Dir1,1)
%         subplot(2,1,1);  %plot(stats.windows(:,1)+0.01, estim, '.'); hold on; 
%         plot(positions(Dir1(kDir1,1):Dir1(kDir1,2),1)/100, positions(Dir1(kDir1,1):Dir1(kDir1,2),2), '.');xlim([Dir1(kDir1,1)/100 Dir1(kDir1,2)/100]) ;
%         subplot(2,1,2) ; imagesc(EvtMatB(sMAxes1Idx,Dir1(kDir1,1):Dir1(kDir1,2) ));
%         pause;
%     end
    
    % Test on sleep
    
    %First make a matrix corresponding to SCEs only. This will speed up calculation
%     nwin = 20; % half window size
%     SCEsIntervals = [putSCEidx' - nwin + 1 putSCEidx' + nwin ];
%     SCEsIntervalsSegments = InIntervals(1 : size(EvtMat,2), SCEsIntervals);
%     EvtMat2 = EvtMat(:,SCEsIntervalsSegments);
     
    EvtMat2 = EvtMat(:,1:5000);                                                                          % On raccourcit la matrice des evts car ça prend beaucoup de temps si je fais la détection des évnts sans la condition sur la dérivée. 
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMat2(n,:)); spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];     % Add neuron ID and starts idx to the spikes sample
    end
    spikes(:,1) = TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
    spikes=sortrows(spikes, 1); 
    [stats1,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1,'Px',Px1); % Try Reconstruction
    figure; imagesc(stats1.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    save('statsReconstructPos_10ms.mat','stats1');
 %   load('statsReconstructPos_10ms.mat'); 
 
%     % test on sleep with EvtStarts
    spikesStart=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    EvtStarts2 = EvtStarts(:,1:5000);
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtStarts2(n,:)) ; %find(2*movmean(EvtStarts(n,:), 2)); 
        spikesStart=[spikesStart ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];     % Add neuron ID and starts idx to the spikes sample
    end
    spikesStart(:,1) = TTLstartsTimes(spikesStart(:,1));                                                                % Transform the starts idx in timings
    spikesStart = sortrows(spikesStart,1);
    [stats1Start,~,~] = ReconstructPosition([],spikesStart,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1,'Px',Px1); % Try Reconstruction
    figure; imagesc(stats1Start.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    save('statsReconstructPosDir1Start_10ms.mat','stats1Start');
%    

    
    % dir2 
    
    [~,lambda2,Px2] = ReconstructPosition(positionsDir2,spikesBDir2,'mode','train','window',wndw2,'nBins',[nBins 1]); % Train Reconstruction
    save('lambda2_10ms.mat', 'lambda2'); save('Px2_10ms.mat', 'Px2'); 
    
    lambda2s = squeeze(lambda2)'; [Maxes2,Maxes2Idx]  = max(lambda2s,[],2); [~,sMAxes2Idx]=sort(Maxes2Idx); lambda2s=lambda2s./Maxes2; % order cells using lambda
    figure; imagesc(lambda2s(sMAxes2Idx,:));
    
    [stats2,~,~] = ReconstructPosition([],spikesBDir2,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda2,'Px',Px2); % Try Reconstruction on behvior data first
    figure; imagesc(stats2.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    [~,estim2]=max(stats2.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    figure; plot(stats2.windows(:,1)+0.01, estim2, '.'); hold on; plot(positionsDir2(:,1), positionsDir2(:,2)*100, '.');     
    
    % Test on sleep
    [stats2,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda2,'Px',Px2); % Try Reconstruction
    figure; imagesc(stats2.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    save('statsReconstructPosDir2_10ms.mat','stats2');
    
    
    % Analyse trajectory
 
    % find SCEs timings 
%     nSyncSCEs = sum(syncSCEs); %nSyncSCEs=floor(nSyncSCEs/4); % we performed testing only on the first part of the trace to gain time
%     SyncSCEsIdx = find(syncSCEs);
    nwinplot = 20;
    
    SumEvtStartsSmooth = sum(EvtStarts,1); SumEvtStartsSmooth1 = movmean(SumEvtStartsSmooth, 5) *5; %SumEvtStartsSmooth2 = movmean(SumEvtStartsSmooth, 10) *10;
 %   figure; plot(SumEvtStartsSmooth); hold on; plot(SumEvtStartsSmooth1)
    
    figure; 
    for kSCE = 2 : size(putSCEtimes,1) %  nSyncSCEs %
        %timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
        timeSCE = putSCEtimes(kSCE); %timing dans la ref ephy
        [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_1] = min(abs(stats1.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_2] = min(abs(stats2.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
              
        subplot(4,2,1); imagesc(stats1.estimations(:, iSCE_wndw2_1 - nwinplot+1  : iSCE_wndw2_1 + nwinplot-1),[0 0.3])
        xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(4,2,3); plot(SumEvtStartsSmooth1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1));
        % hold on; plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); hold off;
        yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
        
        subplot(4,2,5); imagesc(EvtStarts(sMAxes1Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);

        subplot(4,2,7); imagesc(EvtMat(sMAxes1Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
              
        subplot(4,2,2); imagesc(stats2.estimations(:, iSCE_wndw2_2 - nwinplot+1  : iSCE_wndw2_2 + nwinplot-1),[0 0.3])
        xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(4,2,4); plot(SumEvtStartsSmooth1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1));
        % hold on; plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); hold off;
        yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
        
        subplot(4,2,6); imagesc(EvtStarts(sMAxes2Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);

        subplot(4,2,8); imagesc(EvtMat(sMAxes2Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
        
        pause; 
    end  
    
    % now we want to use findReplayScore pour analyser les séquences
    
%     FindReplayScore()
%     
%     
%     kSCE = 1;
%     timeSCE = putSCEtimes((kSCE)); %timing dans la ref ephy
%     [~,iSCE_wndw2] = min(abs(stats.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
%     figure; plot(stats.estimations(:,  iSCE_wndw2 ));
% 

%%

% debugging code    
    
% test on sleep with EvtStarts but with basic average of maximum probability of active neurons
    spikesStart2=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    EvtStarts2 = EvtStarts; %(:,1:5000);
    positionEstimationDir1 = zeros(1,size(EvtStarts2,2)); positionEstimationDir2 = zeros(1,size(EvtStarts2,2));
    for n=1: size(EvtStarts2,2)                                                                             % Loop on neurons
        activeNeurons = find(EvtStarts2(:,n));
        if numel(activeNeurons) == 0
            positionEstimationDir1(n) = NaN; positionEstimationDir2(n) = NaN;
        else
            positionEstimationDir1(n) = mean(Maxes1Idx(activeNeurons)); positionEstimationDir2(n) = mean(Maxes2Idx(activeNeurons));
        end
    end
    nwinplot = 20;  
    SumEvtStartsSmooth = sum(EvtStarts,1); SumEvtStartsSmooth1 = movmean(SumEvtStartsSmooth, 5) *5; %SumEvtStartsSmooth2 = movmean(SumEvtStartsSmooth, 10) *10;    
    figure; 
    for kSCE = 2 : size(putSCEtimes,1) %  nSyncSCEs %
        timeSCE = putSCEtimes(kSCE); %timing dans la ref ephy
        [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_1] = min(abs(stats1.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_2] = min(abs(stats2.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
              
        subplot(3,2,1); plot(positionEstimationDir1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1), '.');
        xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(3,2,3); plot(SumEvtStartsSmooth1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1));
        % hold on; plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); hold off;
        yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
        
        subplot(3,2,5); imagesc(EvtStarts(sMAxes1Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);

%         subplot(4,2,7); imagesc(EvtMat(sMAxes1Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
%         xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
              
        subplot(3,2,2); plot(positionEstimationDir2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1), '.');
        xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(3,2,4); plot(SumEvtStartsSmooth1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1));
        % hold on; plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); hold off;
        yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
        
        subplot(3,2,6); imagesc(EvtStarts(sMAxes2Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);

%         subplot(4,2,8); imagesc(EvtMat(sMAxes2Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
%         xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
        
        pause; 
    end  
    
    
    
% debugging de la reconstruction pour EvtStarts

wndw2=0.01;
%EvtStarts3 = EvtStarts(:,1:10); %        
EvtStarts3 = zeros(111,10);   EvtStarts3(2,:)=2;  EvtStarts3(4,:)=4;     
        spikes3=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
        for n=1:Experiment.nNeurons                                                                             % Loop on neurons
                spike3TimesN = find(EvtStarts3(n,:)); spikes3=[spikes3 ;  spike3TimesN' n*ones(numel(spike3TimesN),1)]; % Add neuron ID and starts idx to the spikes sample
        end
spikes3(:,1)=spikes3(:,1)/100;                                                               % Transform the starts idx in timings
spikes3 = sortrows(spikes3,1);
[stats3,~,~] = ReconstructPositionCathie([],spikes3,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1,'Px',Px1, 'type', 'll'); % Try Reconstruction
figure; imagesc(stats3.estimations);
stats3.windows
    
%  now train on 100% of behavior data and test on sleep
    
    %dir 1
    
    [~,lambda1c,Px1c] = ReconstructPositionCathie(positionsDir1,spikesBDir1,'mode','train','window',wndw2,'nBins',[nBins 1]); % Train Reconstruction
 %   save('lambda_10ms.mat', 'lambda1'); save('Px_10ms.mat', 'Px1'); 
    
    lambda1sc = squeeze(lambda1c)'; [Maxes1c,Maxes1Idxc]  = max(lambda1sc,[],2); [~,sMAxes1Idxc]=sort(Maxes1Idxc); lambda1sc=lambda1sc./Maxes1c; % order cells using lambda
    figure; imagesc(lambda1sc(sMAxes1Idxc,:));
    
    [stats1bc,~,~] = ReconstructPositionCathie([],spikesBDir1,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1c,'Px',Px1c); % Try Reconstruction on behvior data first
    figure; imagesc(stats1bc.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
    [~,estim1]=max(stats1bc.estimations);                                                                   % Define Estimated Position as the bin with the max probability
    figure; plot(stats1bc.windows(:,1)+0.01, estim1, '.'); hold on; plot(positionsDir1(:,1), positionsDir1(:,2)*100, '.');     
%    figure; % choix d'un systeme de temps; on prend celui donne pour l'instant à ReconstructPosition. Pour le sommeil on verra si on veut prendre le vrai temps pour plotter aussi l'ephy
%     for kDir1 = 1 : size(Dir1,1)
%         subplot(2,1,1);  %plot(stats.windows(:,1)+0.01, estim, '.'); hold on; 
%         plot(positions(Dir1(kDir1,1):Dir1(kDir1,2),1)/100, positions(Dir1(kDir1,1):Dir1(kDir1,2),2), '.');xlim([Dir1(kDir1,1)/100 Dir1(kDir1,2)/100]) ;
%         subplot(2,1,2) ; imagesc(EvtMatB(sMAxes1Idx,Dir1(kDir1,1):Dir1(kDir1,2) ));
%         pause;
%     end
    
    % Test on sleep
    
    %First make a matrix corresponding to SCEs only. This will speed up calculation
%     nwin = 20; % half window size
%     SCEsIntervals = [putSCEidx' - nwin + 1 putSCEidx' + nwin ];
%     SCEsIntervalsSegments = InIntervals(1 : size(EvtMat,2), SCEsIntervals);
%     EvtMat2 = EvtMat(:,SCEsIntervalsSegments);
     
    EvtMat2 = EvtMat(:,1:5000);                                                                          % On raccourcit la matrice des evts car ça prend beaucoup de temps si je fais la détection des évnts sans la condition sur la dérivée. 
    spikes=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtMat2(n,:)); spikes=[spikes ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];     % Add neuron ID and starts idx to the spikes sample
    end
    spikes(:,1) = TTLstartsTimes(spikes(:,1));                                                                % Transform the starts idx in timings
    spikes=sortrows(spikes, 1); 
    [stats1c,~,~] = ReconstructPositionCathie([],spikes,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1c,'Px',Px1c); % Try Reconstruction
    figure; imagesc(stats1c.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
 %   save('statsReconstructPos_10ms.mat','stats1');
 %   load('statsReconstructPos_10ms.mat'); 
 
%     % test on sleep with EvtStarts
    spikesStart=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
    EvtStarts2 = EvtStarts(:,1:5000);
    for n=1:Experiment.nNeurons                                                                             % Loop on neurons
        spikeTimesN = find(EvtStarts2(n,:)) ; %find(2*movmean(EvtStarts(n,:), 2)); 
        spikesStart=[spikesStart ;  spikeTimesN' n*ones(numel(spikeTimesN),1)];     % Add neuron ID and starts idx to the spikes sample
    end
    spikesStart(:,1) = TTLstartsTimes(spikesStart(:,1));                                                                % Transform the starts idx in timings
    spikesStart = sortrows(spikesStart,1);
    [stats1Startc,~,~] = ReconstructPositionCathie([],spikesStart,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1c,'Px',Px1c); % Try Reconstruction
    figure; imagesc(stats1Startc.estimations); xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 14);
 %   save('statsReconstructPosDir1Start_10ms.mat','stats1Start');
%    
    nwinplot = 20;
    
    SumEvtStartsSmooth2 = sum(EvtStarts2,1); SumEvtStartsSmooth3 = movmean(SumEvtStartsSmooth2, 5) *5; %SumEvtStartsSmooth2 = movmean(SumEvtStartsSmooth, 10) *10;
 %   figure; plot(SumEvtStartsSmooth); hold on; plot(SumEvtStartsSmooth1)
    
    figure; 
    for kSCE = 2 : size(putSCEtimes,1) %  nSyncSCEs %
        %timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
        timeSCE = putSCEtimes(kSCE); %timing dans la ref ephy
        [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_1] = min(abs(stats1.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
        [~,iSCE_wndw2_2] = min(abs(stats2.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
              
        %subplot(4,1,1); imagesc(stats1c.estimations(:, iSCE_wndw2_1 - nwinplot+1  : iSCE_wndw2_1 + nwinplot-1),[0 0.3]);
        subplot(4,1,1); imagesc(stats1Startc.estimations(:, iSCE_wndw2_1 - nwinplot+1  : iSCE_wndw2_1 + nwinplot-1),[0 0.3])
        xlabel('Time Points (window = 10ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
          
        subplot(4,1,2); plot(SumEvtStartsSmooth3(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1));
        % hold on; plot(SumEvtStartsSmooth2(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); hold off;
        yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
        
        subplot(4,1,3); imagesc(EvtStarts2(sMAxes1Idxc,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);

        subplot(4,1,4); imagesc(EvtMat2(sMAxes1Idxc,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
        xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
        
        pause; 
    end  


%     % verification que les fenetres temporelles sont les bonnes dans le graphe precedent    
%     wndw2 = 0.01;
%     figure; 
%     for kSCE = 2 : size(putSCEtimes,1) %  nSyncSCEs %
%         %timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
%         timeSCE = putSCEtimes(kSCE); %timing dans la ref ephy
%         [~,iSCE] = min(abs(TTLstartsTimes-timeSCE)); %timing en # de frames optiques
%         [~,iSCE_wndw] = min(abs(stats1.windows(:,1)+0.01-timeSCE)); %timing en # de frames optiques
%          
%         %to check that the timing is OK I recalculate the estimated position on the considered time.       
%         EvtMat3 = EvtMat(:,iSCE - nwinplot*(wndw2/0.01)+2  : iSCE + nwinplot*(wndw2/0.01)-1);                                                                            % On raccourcit la matrice des evts car ça prend beaucoup de temps si je fais la détection des évnts sans la condition sur la dérivée. 
%         spikes3=[];                                                                                              % Initialize 'spikes' (starts) sample mat.
%         for n=1:Experiment.nNeurons                                                                             % Loop on neurons
%                 spike3TimesN = find(EvtMat3(n,:)); spikes3=[spikes3 ;  spike3TimesN' n*ones(numel(spike3TimesN),1)]; % Add neuron ID and starts idx to the spikes sample
%         end
%         spikes3(:,1)=TTLstartsTimes(spikes3(:,1));                                                                % Transform the starts idx in timings
%         spikes3 = sortrows(spikes3,1);
%         [stats3,~,~] = ReconstructPosition([],spikes3,'mode','test','window',wndw2,'nBins',[nBins 1],'lambda',lambda1,'Px',Px1); % Try Reconstruction
%        
%         subplot(4,1,1);       
%         imagesc(stats1.estimations(:, iSCE_wndw - nwinplot+1  : iSCE_wndw + nwinplot-1),[0 0.3])
%         xlabel('Time Points (window = 20ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12); title(['kSCE = ' num2str(kSCE,'%d')]); 
%           
%         subplot(4,1,3)
%  %      imagesc(EvtStarts(sMAxes1Idx,iSCE - nwinplot*(wndw/0.01)+1  : iSCE + nwinplot*(wndw/0.01)-1))
%         imagesc(stats3.estimations,[0 0.3])
%         xlabel('Time Points (window = 20ms)'); ylabel('Recontructed position'); set(gca, 'FontSize', 12);
% 
%         subplot(4,1,4)
%         imagesc(EvtMat(sMAxes1Idx,iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1))
%         xlabel('Time Points (at 100Hz)'); ylabel('# neuron (ordered)'); set(gca, 'FontSize', 12);
%         
%         subplot(4,1,2)
%         plot(SumEvtStartsSmooth1(iSCE - nwinplot*(wndw2/0.01)+1  : iSCE + nwinplot*(wndw2/0.01)-1)); 
%         yline(nbStartsThr,'g'); xlabel('Time Points (at 100Hz)'); ylabel('SCE detection'); set(gca, 'FontSize', 12);
%        
%         pause; 
%     end
    
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
    vectSyncSCEs(ePhysChan) = nbSyncSCEs;
    vectnRipples(ePhysChan) = size(ripples, 1);
        
    pause
    close all;
end





