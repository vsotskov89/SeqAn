function Params=GetParams
%% tests
Params.test10 =0; %test dans GetPositionSpeed
%% Folders
Params.DataFolder = 'G:\OldAnalysis_SmallErrorMotionCorr\Souris132830\20-12-15\161508_sCMOS_132830-awake-10ms-220mW\';                        % Folder where the data folders are located
Params.XLSsFolder = '\mnt\cortex-data-97\Calcium\XLSs\';    % Folder where the XLSs are saved
Params.PDFsFolder = 'G:\OldAnalysis_SmallErrorMotionCorr\Souris132830\20-12-15\161508_sCMOS_132830-awake-10ms-220mW\';    % Folder where the PDFs are saved

%% Position and speed
Params.rewardDist = 95;                                 % Distance between rewards in cm (for px to cm conversion)
Params.speedThreshold = 2;                              % if speed > speedThreshold (cm/s) --> mouvement
Params.speedMax = 250;                                   % if speed > speedMax (cm/s) --> artifact

%% Events Detection
Params.Method.SmallEvents = 'improvedSE';                       % 'noSE' ; 'SE' ; 'improvedSE'
    if strcmp(Params.Method.SmallEvents,'improvedSE')   % Parameters for improved SE detection
        Params.Method.filtRange=[0 3];                      % Fluo signal will be filtered in this range     
        Params.Method.minDur = 20;                          % Minimal Duration for an events (in frames)
    end

%% PC criteria
Params.critPCs.thrPFs = 1.2;                % PF is defined as the region over thrPF * meanActivity that contains the max
Params.critPCs.thrPeak = 0.3;               % Peak must exceed thrPeak (fraction of the max theroretical value) ...
Params.critPCs.StabDiscFactor = 0.3;        %   ... minus StabDiscFactor*Stability in its PF (PC : Peak > thrPeak - StabDiscFactor*Stab
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
    Params.Theta.boundsTheta=[4 10];            % Frequency bounds for theta band
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

