function Params=GetParams
%% tests
Params.test10 = 0;
%% Folders
Params.DataFolder = '/export/home1/RawCalciumData/sleeppost/';                        % Folder where the data folders are located
Params.XLSsFolder = '/export/home1/RawCalciumData/sleeppost/';    % Folder where the XLSs are saved
Params.PDFsFolder = '//export/home1/RawCalciumData/sleeppost/';    % Folder where the PDFs are saved

%% Position and speed
Params.rewardDist = 95;                                 % Distance between rewards in cm (for px to cm conversion)
Params.speedThreshold = 2;                              % if speed > speedThreshold (cm/s) --> mouvement
Params.speedMax = 250;                                  % if speed > speedMax (cm/s) --> artifact
Params.lengthThreshold = 10;                            % if distance travelled < lengthThreshold (cm) --> segment discarded
Params.stopThreshold = 1;                               % if stops duration < stopThreshold (s) --> not considered as a real stop
Params.smoothActivitySpace = 9;             % time window for moving average applied to the activity as a function of space just before finding PF bounderies
Params.smoothActivitySpaceSR = 1;             % time window for moving average applied to the activity as a function of space for single runs, before averaging
Params.normActivitySR = 0;                  % 1 to normalize activity for single runs, 0 to not normalize
Params.maxNumberPFs = 3;                     % maximum number of PFs
Params.sizePFsThr = 0.2;                    % Threshold for defining the PFs limits (in multiple of peak height)


%% Events Detection
Params.Method.SmallEvents = 'RisesDerivative';               % 'RisesSE' ; 'RisesImprovedSE'; 'RisesAndDecaysSE'; 'RisesDerivative' ;
    if strcmp(Params.Method.SmallEvents,'RisesImprovedSE')       % Parameters for improved SE detection
        Params.Method.filtRange=[0 3];                          % Fluo signal will be filtered in this range     
        Params.Method.minDur = 20;                              % Minimal Duration for an events (in frames)
    end

                                                        % 'Calcium' : activity will be the smoothed calcium traces (C_raw)
                                                        % 'Derivative' : smoothed derivative of the smoothed calcium traces (C_raw)
Params.Method.ActivityType =  'Binary'; %'Spikes'; % 'Spikes';                 % 'Spikes' : smoothed spikes found by BSD
                                                        % 'Binary' : 1 during rises, 0  elsewhere
                                                        
%% BSD
Params.Obsd = struct; % Struct of experimental conditions & decoding options.
%Params.Obsd.Time  % Number of time frames. Will be assigned later using the data
Params.Obsd.nNeurons = 1; % Number of neurons.
Params.Obsd.dt = 0.01; % interval duration. (s)
Params.Obsd.adaptive = 0; % 0 : Not adaptive. Will use provided values for parameters, and estimate the unknown ones.
Params.Obsd.iterations = 10; % Maximal number of iterations. Default: 5.

Params.Pbsd = struct; % Struct of generative model properties.
Params.Pbsd.tauRise = 0.005; % 0.07; % Fluorescence rise time (s)
Params.Pbsd.tauDecay = 0.1; % 0.8; % Fluorescence decay time (s)
Params.Pbsd.th = 0;  % Threshold on spike amplitude 
Params.Pbsd.b = 0; % Baseline position                                                       
                                                        
%% PC criteria
Params.fractPassagePF = 0.5;                % Fraction of the PF which must be crossed to validate the passage
Params.critPCs.thrPFs = 1.8;  %1.2;                % PF is defined as the region over thrPF * meanActivity that contains the max%Params.critPCs.thrPeak = 0.3;               % Peak must exceed thrPeak (fraction of the max theroretical value) ...
%Params.critPCs.StabDiscFactor = 0.3;        %   ... minus StabDiscFactor*Stability in its PF (PC : Peak > thrPeak - StabDiscFactor*Stab
Params.critPCs.thrStab = 0.35;               % stability threshold
Params.critPCs.thrAct = 1; %2.4;
Params.PrintPDFs = 1;                       % 1 : prints pdf, 0 does not. /!\ Folder must not contain .pdf ; pdfs will be saved in folder then merged and resulting pdf will be moved in containing directory

%% SCE detection
Params.findSCE = 0;                         % 1 : Find SCEs; 0 : do not find SCEs
    Params.SCEs.Window=10;                      % Size of the window in frames
    Params.SCEs.Thr=4;                          % Number of event's starts
    Params.SCEs.maxInterval=0;                  % Merge SCEs closer than this interval in frames
    Params.SCEs.PrintPDFs = 0;                  % Create summary PDFs
    % Params.SCEs.maxSize=200;
    % Params.SCEs.excludeCorr = 1;   
    
%% Theta detection
Params.findTheta = 0;                       % 1 : Find Theta; 0 : do not find Theta
    Params.Theta.boundsTheta=[5 10];            % Frequency bounds for theta band
    Params.Theta.minOverlap=0.2;                % fraction of events that must overlap with a theta region to decide if event has theta
    Params.Theta.Durations=[350 300 10000];     % [interval, min , max] : minimum inter-ripple interval, and minimum and maximum ripple durations, in ms
    Params.Theta.Thresholds=[2 3];              % [bounds peak] theta regions : Theta power must be over bounds*std and reach at least once peak*std
    Params.Theta.Normalization=0;               % Normalize Theta Power by power in an other frequency band   
    Params.Theta.boundsFreqNorm=[10 17];        % Frequency bounds for the normalization

%% Events Statistics
Params.getEvtsStats=1;                      % 1 : Get Events stats; 0 : do not get them
    Params.EvtsCategories = [2 4];          % [bound1 bound2] : 3 categories of events : < bound1*std ; bound1*std > and < bound2*std ; > bound2*std  

%% XLS summary
Params.MakeXLS=0;                           % 0 : do not ; 1 : do

%% Save Matrices of results (Neurons, Experiment, Params)
Params.SaveMats=1;                           % 0 : do not ; 1 : do

%% Aligning LFP
Params.alignLFP=0;                          % 0 : do not ; 1 : do
    Params.LFP.TTLchan = 19;                % 19 ; should not be changed a priori
    Params.LFP.FreqDat = 20000;             % sampling freq (in Hz) of the .dat file (20kHz)
    Params.LFP.CheckPlot = 0;               % Plot to check if Experiments start/stop are well placed on the TTL recording 
    Params.LFP.chanLFP=4;                   % Choose LFP Channel
%% Experiments to analyze
Params.Paths=uigetdir2(Params.DataFolder,'Select Folder to Analyze');   % Choose the experiment(s) to analyse

