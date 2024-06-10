% given events, signal : 
% events=Neurons.Rises(2).Starts'/100; signal=[[1:Experiment.SizeData]'/100 results.C_raw(2,:)'];


durations = [-1 1];
resolution = 1;
dt = mode(diff(signal(:,1)));
indices = round((durations(1):dt:durations(2))/dt);
spectrograms = cell(length(events),1);
for i=1:size(events,1)
    start = find(signal(:,1)>events(i) + durations(1),1); if start+indices(1)<1, continue; end
    [w,wt,f0] = WaveletSpectrogramRaw(signal(start+indices,:),'range',[1 50],'resolution',resolution);
    spectrograms{i,1} = w;
    
end
% spectrograms(cellfun(@isempty,spectrograms))=[];
stackedSpectrograms = cat(3,spectrograms{:});

%%

subplot(1,2,1);
power = mean(stackedSpectrograms,3);
% Get equally spaced frequencies (rather than logscale)
f = (min(f0):min(diff(f0)):max(f0))';
means = interp1(f0,power,f);

wt = linspace(durations(1),durations(2),length(indices));
PlotColorMap(means,'x',wt,'y',f);
title('mean spectrogram');

set(get(colorbar,'YLabel'),'String','power');
xlabel('time (s)');
axis square
ylabel('frequency (Hz)');

subplot(1,2,2);
power = median(stackedSpectrograms,3);
medians = interp1(f0,power,f);
PlotColorMap(medians,'x',wt,'y',f);
title('median spectrogram');
set(get(colorbar,'YLabel'),'String','power');
xlabel('time (s)');
axis square
ylabel('frequency (Hz)');