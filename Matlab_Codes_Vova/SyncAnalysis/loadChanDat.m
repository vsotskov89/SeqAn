function data=loadChanDat(channels,inputFile)

% [file,path] = uigetfile([path '/*.dat'],'Select ePhys data to load');

% inputFile = fullfile(path,file);
disp(['Loading : ' inputFile])
data= LoadBinary(inputFile,'channels',channels,'nChannels',24);
