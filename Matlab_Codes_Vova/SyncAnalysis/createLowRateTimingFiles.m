function createLowRateTimingFiles(expPaths,nDownsample,tagLowRate)

ExpPhases={'SleepPRE','Behavior', 'SleepPOST'};
expNames={};
for p=1:numel(ExpPhases)
    phase = ExpPhases{p};
    for n=1:numel(expPaths.(phase).TimingsCMOS)
        expNames=[expNames {[phase num2str(n)]}];
    end
end

remainder = 0;
for exp=1:numel(expNames) % loop on all the experiments that have been selected in the main program 
    try
        % read sCMOS timing files of the current experiment (exp)    
        timingsCMOS=csvread([expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1} '\sCMOS_AbsoluteTimingsSec.csv']);
    catch
        FluoFileSec=expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1};
        [newFluoFileSec,path]=uigetfile(FluoFileSec,'Select absolute Fluo Timings file in Sec');
        newFluoFileSec=fullfile(path,newFluoFileSec); movefile(newFluoFileSec,FluoFileSec);
        timingsCMOS=csvread(FluoFileSec);
    end
    ExpTimePoints=numel(timingsCMOS)-4;
    timingsLowRate = timingsCMOS(5+remainder : nDownsample : end);
    %csvwrite([expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1} '\sCMOS_AbsoluteTimingsSec' tagLowRate '.csv'],timingsLowRate) 
    writematrix(timingsLowRate,[expPaths.(expNames{exp}(1:end-1)).TimingsCMOS{str2double(expNames{exp}(end)),1} '\sCMOS_AbsoluteTimingsSec' tagLowRate '.csv']) 
    remainder = mod(ExpTimePoints,nDownsample);
    
end

