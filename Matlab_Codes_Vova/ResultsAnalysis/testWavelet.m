% given events, signal

region='Starts';
% region='Ends';
N=find([Neurons.isPC.Dir1] + [Neurons.isPC.Dir2]);

durations = [-1 2];
resolution = 1;
range=[1 20];
dt = 0.01; indices = round((durations(1):dt:durations(2))/dt);

% spectrograms = cell(sum(arrayfun(@(S)numel(S.(region)),Neurons.Rises(N))),1); % All events
% spectrograms = cell(sum(arrayfun(@(S)sum(S.HasTheta),Neurons.Rises(N))),1); % Events with Theta
% spectrograms = cell(sum(arrayfun(@(S)sum(S.inPF),Neurons.Rises(N))),1); % Events in PF
spectrograms = cell(sum(arrayfun(@(S)sum(S.inPF & S.EvtsCat==3),Neurons.Rises(N))),1); % Big Events in PF

count=0;
for neuron = N
disp(neuron); count=count+1;

% events=Neurons.Rises(neuron).(region)'/100;
% events=Neurons.Rises(neuron).(region)(logical(Neurons.Rises(neuron).HasTheta))'/100; % For events with theta theta
% events=Neurons.Rises(neuron).(region)(logical(Neurons.Rises(neuron).inPF))'/100; % For events in PF
events=Neurons.Rises(neuron).(region)(Neurons.Rises(neuron).inPF & Neurons.Rises(neuron).EvtsCat==3)'/100; % For Big events in PF

signal=[[1:Experiment.SizeData]'/100 results.C_raw(neuron,:)']; 

    for i=1:size(events,1)
        start = find(signal(:,1)>events(i) + durations(1),1); if start+indices(1)<1, continue; end
%         [w,~,f0] = WaveletSpectrogramRaw(signal(start+indices,:),'range',range,'resolution',resolution);
        [w,f0]=cwt(signal(start+indices,2),100); freqs=find(f0<25);f0=flipud(f0(freqs)); w=abs(w(freqs,:));
        spectrograms{i+count,1} = w/median(w(:));
        count=count+1;
    end
end
stackedSpectrograms = cat(3,spectrograms{:});

%%
figure
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