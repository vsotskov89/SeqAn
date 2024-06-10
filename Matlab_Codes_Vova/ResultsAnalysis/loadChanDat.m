function data=loadChanDat(channels,path)

[file,path] = uigetfile({'*.dat' 'ePhys Data'},'Select .dat to load channel from',[path 'ePhy']);

inputFile = fullfile(path,file);
disp(['Loading : ' inputFile])
data= LoadBinary(inputFile,'channels',channels,'nChannels',24);
