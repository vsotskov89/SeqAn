function [results,timingsCMOS,nameExtract]=ExtractExpData(expPaths,MainPath,DeepCAD,nDownsample,tagLowRate)
% Does NOT work as a data splitter (which would have been easier maybe),
% but if you can't split the data by experiments, you CAN extract "blocks" 
% of (consecutive or not) experiment data (e.g. : all sleep POST session). 
% If you do want to split the data you just have to run the program multiple
% times 

ExpPhases={'SleepPRE','Behavior', 'SleepPOST'};
expNames={};
for p=1:numel(ExpPhases)
    phase = ExpPhases{p};
    for n=1:numel(expPaths.(phase).TimingsCMOS)
        expNames=[expNames {[phase num2str(n)]}];                       % Cell array representing the list of all experiments names. Ex :  {'SleepPRE1'}    {'Behavior1'}
    end
end

[idxExtract,~] = listdlg('ListString',expNames,'PromptString','Select Experiment to extract :','ListSize',[200 100],'Name','List of Experiments','OKString','Extract'); %choose which experiment to load
nameExtract= expNames(idxExtract); % cell array of the names of the experiments to extract. Ex :  {'SleepPRE1'}. 

timingsCMOS=[]; 
TotalTimePoints=0;
TimeIdx=ones(numel(expNames),2);
for exp=1:numel(expNames) % loop on all the experiments that have been selected in the main program 
    try
        % read sCMOS timing files of the current experiment (exp)    
        timingsCMOS.(expNames{exp})=csvread([expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1} '\sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']);
    catch
        FluoFileSec=expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1};
        [newFluoFileSec,path]=uigetfile(FluoFileSec,'Select absolute Fluo Timings file in Sec');
        newFluoFileSec=fullfile(path,newFluoFileSec); movefile(newFluoFileSec,FluoFileSec);
        timingsCMOS.(expNames{exp})=csvread(FluoFileSec);
    end
    if nDownsample == 1
        ExpTimePoints=numel(timingsCMOS.(expNames{exp}))-4;
    else
        ExpTimePoints=numel(timingsCMOS.(expNames{exp}));
    end
    TimeIdx(exp,:)= TimeIdx(exp,:)+[TotalTimePoints TotalTimePoints+ExpTimePoints-1];
    TotalTimePoints=TotalTimePoints+ExpTimePoints;
end

ConcResults=load([expPaths.Concatenated.Calcium{1} '/Results' DeepCAD tagLowRate '.mat']);ConcResults=ConcResults.results;

if size(ConcResults.C_raw,2) ~= TotalTimePoints
    error('The number of Timings points and the number of frames are different (but why ... ?)')
end

toExtract=[];
for i = 1:numel(idxExtract)
    toExtract=[toExtract TimeIdx(idxExtract,1):TimeIdx(idxExtract,2)];
end

results=ConcResults;results.C_raw=results.C_raw(:,toExtract);
results.C=results.C(:,toExtract); results.S=results.S(:,toExtract);
results.AvgRoi=results.AvgRoi(:,toExtract); results.AvgRoiNoBaseline=results.AvgRoiNoBaseline(:,toExtract);

SaveExtractedData= questdlg('Do you want to save extracted data ?','Extracted data Saving','Yes','No','Yes');
switch SaveExtractedData
    case 'Yes'
       % uisave('ExtractedResults',[MainPath '/ExtractedResults_' cell2mat(nameExtract)])
        uisave('results',[MainPath  'Results.mat'])
end






