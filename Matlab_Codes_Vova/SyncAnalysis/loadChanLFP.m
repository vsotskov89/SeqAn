function lfp=loadChanLFP(channels,path,full)

if nargin<3
    [file,path] = uigetfile({'*.dat;*.lfp' 'LFP file or Data to downsample';'*.lfp' 'LFP file'; '*.dat' 'Data to downsample'},'Select .lfp or .dat to downsample to lfp',[path 'ePhy']);
    inputFile = fullfile(path,file);
else
    inputFile = path;
end
if strcmp(inputFile(end-3:end),'.dat')
    dat2lfp(inputFile)
    
    inputFile=[inputFile(1:end-3) 'lfp'];
end
try
    lfp = LoadBinary([inputFile(1:end-3) 'lfp'],'channels',channels,'nChannels',24);
catch
    inputFile=[inputFile(1:end-3) 'dat'];
    dat2lfp(inputFile)
    lfp = LoadBinary([inputFile(1:end-3) 'lfp'],'channels',channels,'nChannels',24);
end

t = 0:1/1250:(0+(length(lfp)-1)/1250);t=t';
lfp = [t lfp];

