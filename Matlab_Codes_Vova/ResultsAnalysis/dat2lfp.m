function dat2lfp(inputFile)
nChannels = 24; up = 1 ; down = 16;

ResampleBinary(inputFile,[inputFile(1:end-3) 'lfp'],nChannels,up,down)
