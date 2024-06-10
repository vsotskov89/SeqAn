function [Experiment,results,currentResultsFile]=LoadData(Experiment,currentResultsFile,results)

Experiment.session='';                                             % Init session
file =strsplit(Experiment.path,'/'); file = file{end-1};           % Get file name
Experiment.file=file;

% if contains(file,'AllExperiments')                      % Files from Experiments with Sleep have a session prefix
%     Experiment.session='awake_';
% end

% Paths to files                                                       
Experiment.PositionFile=[Experiment.path Experiment.session 'mousePos1.mat'];                                                           % .csv position file
Experiment.BehaviorTimingFile=[Experiment.path Experiment.session 'Basler_AbsoluteTiming.csv'];                                         % .csv file : "Basler_AbsoluteTiming"
Experiment.BehaviorTimingFileSec=[Experiment.path Experiment.session 'Basler_AbsoluteTimingSec.csv'];                                   % .csv file : "Basler_AbsoluteTimingSec"
Experiment.rBehaviorTimingFile=[Experiment.path Experiment.session 'Basler_RelativeTiming.csv'];                                      % .csv file : "Basler_RelativeTiming"
Experiment.FluoTimingFile=[Experiment.path Experiment.session 'sCMOS_AbsoluteTimings.csv'];                                             % .csv file : "sCMOS_AbsoluteTimings"
Experiment.FluoTimingFileSec=[Experiment.path Experiment.session 'sCMOS_AbsoluteTimingsSec.csv'];                                       % .csv file : "sCMOS_AbsoluteTimingsSec"
Experiment.rFluoTimingFile=[Experiment.path Experiment.session 'sCMOS_RelativeTimings.csv'];                                            % .csv file : "sCMOS_RelativeTimings"
Experiment.ResultFile =[Experiment.path Experiment.session 'Results.mat'];                                                              % results file from the Analysis

% Load Data
if isequal(currentResultsFile,file)                                     
    disp('    Results file already loaded')
else
    disp('    Loading new results file : in progress')
    try
        load(Experiment.ResultFile)
        fprintf([repmat('\b',1,12) '%s\n'],'Done')
    catch
        fprintf([repmat('\b',1,12) '%s\n'],'Unable to load')
        disp('        Manual selection requested')
        newResultsFile=uigetfile(Experiment.path,'Select Results file'); newResultsFile=fullfile(Experiment.path,newResultsFile);
        movefile(newResultsFile,Experiment.ResultFile)
        load(Experiment.ResultFile)
        disp('    Results file renamed and loaded')
    end
    
end
if and(strcmp(Experiment.session,'awake_'),exist('subresults','var'))                      % Files from Experiments with Sleep are named subresults
    results=subresults;
end
currentResultsFile=file;                                                        % Keep track of the loaded results file

% Load Timings files
try
%     if and(isfile(Experiment.FluoTimingFile),isfile(Experiment.BehaviorTimingFile))                       % If Timings Files exist in Sec : load them
%         Experiment.timeFluo = csvread(Experiment.FluoTimingFile);
%         Experiment.timeBasler =csvread(Experiment.BehaviorTimingFile);
%     else                                                                            % Else load Timing Files in machine time
%        FluoTimingFile2=Experiment.FluoTimingFile ; FluoTimingFile2(end-6:end-4)=[]; 
%        BehaviorTimingFile2=Experiment.BehaviorTimingFile; BehaviorTimingFile2(end-6:end-4)=[];
%        Experiment.timeFluo = csvread(FluoTimingFile2);
%        Experiment.timeBasler =csvread(BehaviorTimingFile2);
%     end
%     if and(isfile(Experiment.rFluoTimingFile),isfile(Experiment.BehaviorTimingFile))                      % If Relative Timings file exists, use it to correct the Absolute Timings File if it's in Sec
%         Experiment.rtimeFluo=csvread(Experiment.rFluoTimingFile);
%         t=zeros(size(Experiment.timeFluo,1),size(Experiment.timeFluo,2));
%         t=t+Experiment.timeFluo(1)+Experiment.rtimeFluo;
%         Experiment.timeFluo=t;
%     end
    if and(exist(Experiment.FluoTimingFileSec,'file'),exist(Experiment.BehaviorTimingFileSec,'file'))
        Experiment.timeFluo = csvread(Experiment.FluoTimingFileSec);
        Experiment.timeBasler =csvread(Experiment.BehaviorTimingFileSec);
    else
        Experiment.timeFluo = csvread(Experiment.FluoTimingFile);
        Experiment.timeBasler =csvread(Experiment.BehaviorTimingFile);
    end
    if exist(Experiment.rFluoTimingFile,'file')
        rtimeFluo=csvread(Experiment.rFluoTimingFile);
        Experiment.timeFluo=rtimeFluo+Experiment.timeFluo(1);
    end
    if exist(Experiment.rBehaviorTimingFile,'file')
        rtimeBehav=csvread(Experiment.rBehaviorTimingFile);
        Experiment.timeBasler=rtimeBehav+Experiment.timeBasler(1);
    end
catch
    if ~exist(Experiment.FluoTimingFile,'file'); newFluoFile=uigetfile('*.csv','Select absolute Fluo Timings file',Experiment.path);
        if newFluoFile; newFluoFile=fullfile(Experiment.path,newFluoFile); movefile(newFluoFile,Experiment.FluoTimingFile); end 
    end
    if ~exist(Experiment.FluoTimingFileSec,'file'); newFluoFileSec=uigetfile('*.csv','Select absolute Fluo Timings file in Sec',Experiment.path);
        if newFluoFileSec; newFluoFileSec=fullfile(Experiment.path,newFluoFileSec); movefile(newFluoFileSec,Experiment.FluoTimingFileSec); end
    end
    if ~exist(Experiment.rFluoTimingFile,'file'); newFluoFileRel=uigetfile('*.csv','Select relative Fluo Timings file',Experiment.path);
        if newFluoFileRel; newFluoFileRel=fullfile(Experiment.path,newFluoFileRel); movefile(newFluoFileRel,Experiment.rFluoTimingFile); end
    end
    if ~exist(Experiment.BehaviorTimingFile,'file'); newBehavFile=uigetfile('*.csv','Select absolute Behavior Timings file',Experiment.path);
        if newBehavFile; newBehavFile=fullfile(Experiment.path,newBehavFile); movefile(newBehavFile,Experiment.BehaviorTimingFile); end
    end
    if ~exist(Experiment.BehaviorTimingFileSec,'file'); newBehavFileSec=uigetfile('*.csv','Select absolute Behavior Timings file in Sec',Experiment.path);
        if newBehavFileSec; newBehavFileSec=fullfile(Experiment.path,newBehavFileSec); movefile(newBehavFileSec,Experiment.BehaviorTimingFileSec); end
    end
    if ~exist(Experiment.rBehaviorTimingFile,'file'); newBehavFileRel=uigetfile('*.csv','Select relative Behavior Timings file',Experiment.path);
        if newBehavFileRel; newBehavFileRel=fullfile(Experiment.path,newBehavFileRel); movefile(newBehavFileRel,Experiment.rBehaviorTimingFile); end
    end
    
    if and(exist(Experiment.FluoTimingFileSec,'file'),exist(Experiment.BehaviorTimingFileSec,'file'))
        Experiment.timeFluo = csvread(Experiment.FluoTimingFileSec);
        Experiment.timeBasler =csvread(Experiment.BehaviorTimingFileSec);
    else
        Experiment.timeFluo = csvread(Experiment.FluoTimingFile);
        Experiment.timeBasler =csvread(Experiment.BehaviorTimingFile);
    end
    if exist(Experiment.rFluoTimingFile,'file')
        rtimeFluo=csvread(Experiment.rFluoTimingFile);
        Experiment.timeFluo=rtimeFluo+Experiment.timeFluo(1);
    end
    if exist(Experiment.rBehaviorTimingFile,'file')
        rtimeBehav=csvread(Experiment.rBehaviorTimingFile);
        Experiment.timeBasler=rtimeBehav+Experiment.timeBasler(1);
    end
end

Experiment.SizeData = size(Experiment.timeFluo(5 : end), 1);                      % Size of Fluo Data (nb of frames)
Experiment.nNeurons = size(results.C, 1);    
