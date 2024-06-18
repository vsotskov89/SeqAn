root = '/export/home1/RawCalciumData/';
ses_name = 'awake';
ePhysFile = '/export/home1/Sequences/2_Awake/continuous.dat';
 
%gathering files
load([root, ses_name, '_results.mat'])   %file with CNMF results
TTLData= single(LoadBinary(ePhysFile,'channels', 18, 'nChannels', 24));
TTLstartsTimes = GetTTLtimes(TTLData,results);  
%% Calculate events by Derivative method

SizeData = size(results.C, 2); Experiment.SizeData = SizeData;                                     % Build false Experiment Struct.,                
nNeurons = size(results.C, 1); Experiment.nNeurons = nNeurons;

Params.Method.ActivityType =  'Binary';  %'Spikes';                                                                     % false Params Struct.,
Params.Method.SmallEvents = 'RisesDerivative';  % 'RisesSE'    'RisesSEbsd'  'RisesDecaysSEbsd'             % Choose the method to extract the rises

Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),...   % and false Neurons Struct. that are required for GetRises
    'Ends',cell(Experiment.nNeurons,1),'Duration',cell(Experiment.nNeurons,1)));
    
Neurons=GetRises(results,Params,Experiment,Neurons);

%% Calculating some stats
summ = 0;
for i =1:numel(Neurons.Rises(1).Starts)
    summ = summ + Neurons.Rises(1).Ends(i) - Neurons.Rises(1).Starts(i);
end
disp(summ)
disp(sum(Neurons.Rises(1).Duration))
disp(sum(Neurons.Rises(1).Matrix))
%disp(sum(Neurons.Rises(1).EvtPos1))

%% Get some info
nNeurons = size(results.F0,2);
imsize = size(results.Cn);

%% Save results in a readable form
% Save .csv traces with a timestmap
csvwrite([root, ses_name, '_traces.csv'], [TTLstartsTimes, results.C_raw'])
% Save .mat files for further CellReg matching

A = reshape(full(results.A'), nNeurons, imsize(1), imsize(2));
save([root, ses_name, '_spatials.mat'],'A')




%% Draw traces and cell footprints
figure, hold on
for i = 1:nNeurons
    plot(TTLstartsTimes, results.C_raw(i,:)/max(results.C_raw(i,:)) + i-1, 'Color', sd_colornum_metro(i-1), 'LineWidth', 2);
end
savefig([root, ses_name, '_traces.fig'])
%% Draw traces without a timestamp
nNeurons = 20; %size(results.F0,2);
nFrames = size(results.C,2);
TTimes = linspace(1, nFrames, nFrames);
figure, hold on
for i = 1:nNeurons
    %EvtPos = find(Neurons.Rises(i).EvtPos(1:nFrames));
    %EvtPos1 = find(Neurons.Rises(i).EvtPos1(1:nFrames*10))/10;
    plot(TTimes, results.C_raw(i,1:nFrames)/max(results.C_raw(i,1:nFrames)) + i-1, 'Color', sd_colornum_metro(i-1), 'LineWidth', 2);
    %scatter(EvtPos, ones(numel(EvtPos))*(i-1) - 0.2, 20, sd_colornum_metro(i-1), 'filled')
    %scatter(EvtPos1, ones(numel(EvtPos1))*(i-1) - 0.3, 20, sd_colornum_metro(i-1), 'filled')
    scatter(Neurons.Rises(i).Starts, ones(numel(Neurons.Rises(i).Starts))*(i-1) - 0.2, 20, sd_colornum_metro(i-1), 'filled')
    %scatter(Neurons.Rises(i).Ends, ones(numel(Neurons.Rises(i).Ends))*(i-1) - 0.2, 20, sd_colornum_metro(i-1), 'filled')
end
% savefig('/export/home1/RawCalciumData/awake_traces_with_events.fig')

%% The same for several files stored in different folders
root = '/export/home1/RawCalciumData/';
ses_name = 'sleeppost';
files = dir([root, ses_name, '/CNMF*_070_*/Results.mat']);
ePhysFile = '/export/home1/Sequences/3_SleepPOST/continuous.dat';
sigma = 0;

TTLData= single(LoadBinary(ePhysFile,'channels', 18, 'nChannels', 24));


for f = 1:numel(files)
    load([files(f).folder, filesep, files(f).name])
    TTLstartsTimes = GetTTLtimes(TTLData,results); 

    nNeurons = size(results.F0,2);
    imsize = size(results.Cn);
    par_name = strsplit(files(f).folder, 'results');

    csvwrite([root, ses_name, par_name{2}, '_traces.csv'], [TTLstartsTimes, results.C_raw'])
    A = reshape(full(results.A'), nNeurons, imsize(1), imsize(2));
    %Do some smoothing
    if sigma
        A = imgaussfilt(A,sigma);
    end

    save([root, ses_name, par_name{2}, '_spatials.mat'],'A')

    figure, title([root, ses_name, par_name{2}])
    subplot(1,2,1), imagesc(squeeze(sum(A, 1)))
    subplot(1,2,2), hold on
    for i = 1:nNeurons
        plot(TTLstartsTimes, results.C_raw(i,:)/max(results.C_raw(i,:)) + i-1, 'Color', sd_colornum_metro(i-1), 'LineWidth', 2);

    end
end
