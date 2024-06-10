function [CaLFP,LFP,SCEs,ripples,Experiment,Params,subresults]=SleepAnalysis()

%% Find sessions
Folder = 'D:\For Analysis\';
Path=uigetdir(Folder,'Select experiment folder'); Path=[Path '\'];

Files = dir([Path '*.csv']); Files=struct2cell(Files); Files=Files(1,:);
Files=cellfun(@(C)strsplit(C,'_'),Files,'UniformOutput',0); Files=cellfun(@(C)cell2mat(C(1)),Files,'UniformOutput',0);

Sessions=unique(Files); Sessions(cellfun(@(C)isequal(C,'awake'),Sessions))=[]; 
Sessions={Sessions{contains(Sessions,'sleeppre')} Sessions{contains(Sessions,'sleeppost')}};
disp('Sessions found : ');disp(Sessions)

%% Choose Session to analyze
screenSize = get(0,'screensize');
fSize=[0.45*screenSize(3) 0.45*screenSize(4) 110 30+(numel(Sessions)+1)*16];
f=uifigure('Position',fSize);
bg = uibuttongroup(f,'Position',[5 5 100 85]); f.CloseRequestFcn=@(f,~)closeFig(f,bg);

BtnNames=cell(1,numel(Sessions));
for i = 1:numel(Sessions)
BtnNames{i} = strcat('btn_',Sessions{i});
Btns.(BtnNames{i}) = uiradiobutton(bg,'Position',[10 (numel(Sessions)-i+1)*16 80 15],'Text',Sessions{i});
end
global sess;

waitfor(f);
disp(['Starting Analyzing session ' sess ' from ' Path])
%% Load Results 
Params.Method.SmallEvents = 'improvedSE'; Params.Method.filtRange=[0 3]; Params.Method.minDur = 20;

Experiment.path=Path;
Experiment.Session=sess;
load([Experiment.path '\' Experiment.Session '_Results'],'subresults')
Experiment.nNeurons = size(subresults.C, 1); 
Experiment.SizeData = size(subresults.C, 2);

%% Find Neuronal Events
Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),'Ends',cell(Experiment.nNeurons,1),...
    'Duration',cell(Experiment.nNeurons,1)));
Neurons = GetRises(subresults,Params,Experiment,Neurons);

%% Find SCEs
Params.SCEs.Window=10;                      % Size of the window in frames
Params.SCEs.Thr=3;                          % Threshold for number of event's starts (>Thr*mean)
Params.SCEs.maxInterval=0;                  % Merge SCEs closer than this interval in frames
Params.SCEs.PrintPDFs = 0;                  % Create summary PDFs   
Params.SCEs.PCinfo=0;
Params.SCEs.CheckSleep=0;

if Params.SCEs.PCinfo
    %%%%% do PC analysis
else
    Experiment.PCs.Sleep=1:Experiment.nNeurons;
    Direction='Sleep';
end
if Params.SCEs.CheckSleep
    %%%%% do speed check to verify if mouse is moving
else
    Experiment.Still.Segments=ones(1,Experiment.SizeData);
end

SCEs=[];
SCEs=GetSCEsTimings(SCEs,Neurons,Experiment,Params,Direction);

%% Get Aligned LFP
Params.LFP.TTLchan = 19;                % 19 ; should not be changed a priori
Params.LFP.FreqDat = 20000;             % sampling freq (in Hz) of the .dat file (20kHz)
Params.LFP.CheckPlot = 0;               % Plot to check if Experiments start/stop are well placed on the TTL recording
Params.LFP.chooseChan=0;
Params.LFP.chanOrder=[9 8 10 7 11 6 12 5 13 4 14 3 15 2 16 1];
    
if Params.LFP.chooseChan
    nbRip=zeros(2,16);
    LFPpath=uigetfile({'*.lfp' 'LFP file';'*.lfp' 'LFP file'},'Select .lfp file',[Experiment.path 'ePhy']);
    LFPpath=[Experiment.path 'ePhy\' LFPpath];
    for chan =1:16
        LFP=loadChanLFP(chan,LFPpath,1);                       % Load LFP
        filtLFP= FilterLFP(LFP,'passband','ripples');
        [ripples,~,~] = FindRipples(filtLFP,'thresholds',[2 3],'show','off');
        nbRip(1,chan)=size(ripples,1); nbRip(2,chan)=median(ripples(:,4));
    end
    [~,Params.LFP.chanLFP]=max(nbRip(1,:));
else
    Params.LFP.chanLFP =6;
end

[CaLFP,LFP]=AlignLFP(subresults,Params,Experiment);
%% Get SPW-Rs
% Params.Ripples.stdThr=[2 3];
% filtLFP= FilterLFP(LFP,'passband','ripples');    
% [ripples,~,~] = FindRipples(filtLFP,'thresholds',Params.Ripples.stdThr,'show','on');
ripples = FindRipples2([6 12],[100 250],LFPpath);

ripples(ripples(:,1)<CaLFP(1,1) | ripples(:,1)>CaLFP(end,1),:)=[];
ripplesStarts=zeros(1,size(ripples,1));
for r=1:size(ripples,1)
    ripplesStarts(r)=find(CaLFP(:,1)>=ripples(r,1),1);
end
%% Annexe function
function closeFig(f,bg); sess=bg.SelectedObject.Text; delete(f); end
end