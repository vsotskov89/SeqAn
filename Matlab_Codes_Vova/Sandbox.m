%addpath(genpath('E:\SCRIPTS\CNMF_E-master'))
%cnmfe_setup;
addpath(genpath('C:\SCRIPTS\NoRMCorre-master'));

root = 'D:\\Work\\FiberData\\';
fnames = dir([root, '*.tif']);


DAT = uint16(read_file(strcat(fnames(1).folder, '\', fnames(1).name)));
for f = 1:3 %length(fnames)
    DAT = cat(3,DAT,uint16(read_file(strcat(fnames(1).folder, '\', fnames(1).name))));
end
%% Playing normalized movie
DAT_im = double(DAT)/max(max(max(double(DAT))));
implay(DAT_im)

%% Calculating 
    