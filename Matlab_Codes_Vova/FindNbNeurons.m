%Report number of neurons
root = 'D:\\Sequences\';
files = dir([root, '*\\*esults*.mat']);
for f = 1:length(files)
    load([files(f).folder, '\\', files(f).name])
    nbNeurons = size(results.C, 1);
    nbFrames = size(results.C, 2);
    disp([files(f).folder, '\\', files(f).name, '\tneurons:', num2str(nbNeurons), '\tframes:', num2str(nbFrames)])
    clear results
end