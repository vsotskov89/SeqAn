function Neurons=testFindTheta2(Neurons,results,Experiment,Params)

for n=1:Experiment.nNeurons
fRes=[(1:numel(results.C_raw(n,:)))'/10 results.C_raw(n,:)'];                           % Data to FMA samples

fResNorm=FilterLFP(fRes,'passband',12.5*[10 20]); NSSnorm=testGetNSS(fResNorm(:,2))';   % Get NSS in control band 10-20Hz
fRes = FilterLFP(fRes,'passband',12.5*[5 10]);NSS=testGetNSS(fRes(:,2))';               % Get NSS in theta band 5-10Hz

nNSS=NSS-NSSnorm; nNSS=(nNSS-mean(nNSS))/std(nNSS);                                     % Get the difference between them /!\ Not a very good way to make a normalized NSS
durations=[300 500];                                                                    % [minInter minLength] : merge periods separated by less than minInter and discard periods shorter than minLength
thresholds=[1.5 3];                                                                     % [LowThr HighThr] : theta periods = over LowThr with at least one value over HighThr


%%
time=(1:numel(results.C_raw(n,:)))/100;                                                  % Get timestamps
thresholded=nNSS>thresholds(1);                                                         % Get periods over LowThr
start = find(diff(thresholded)>0); stop = find(diff(thresholded)<0);                    % Find their starts and ends

if length(stop) == length(start)-1;	start = start(1:end-1); end                         % Discard last start if period is incomplete
if length(stop)-1 == length(start); stop = stop(2:end); end                             % Discard first stop is first period is incomplete
if start(1) > stop(1); stop(1) = [];start(end) = []; end                                % Discard both if both are incomplete and realign

Candidates = [start;stop]';                                                             % Get candidate periods matrix

mergedCands = Candidates;                                                               % Start the merging of close periods
inter = time(mergedCands(2:end,1)) - time(mergedCands(1:end-1,2));                      % Compute inter-period durations
toMerge = inter<durations(1)/1000 ;                                                     % Find periods to merge
while any(toMerge)                                                                      % Loop until every period period that need merging is merged
    merge1 = strfind([0 toMerge],[0 1])';                                                   % Find peridods that need merging : 1st ones
    merge2 = merge1+1;                                                                      % The second periods are the next ones
    mergedCands(merge1,2) = mergedCands(merge2,2);                                          % Merge pairs
    mergedCands(merge2,:) = [];                                                             % Delete the now redundant 2nd ones
    inter = time(mergedCands(2:end,1)) - time(mergedCands(1:end-1,2));                      % Recompute inter-period durations
    toMerge = inter<durations(1)/1000;                                                      % Find if new periods need to be merged
end

ThetaEvts = [];                                                                         % Initialize Matrix
for i = 1:size(mergedCands,1)                                                           % Loop on the candidate periods
	[maxValue, posMV] = max(nNSS(mergedCands(i,1):mergedCands(i,2)));                       % Find their max NSS value
	if maxValue > thresholds(2)                                                             % If it's over HighThr
        ThetaEvts = ...                                                                     % These periods is defined as a theta period
            [ThetaEvts ; mergedCands(i,1) mergedCands(i,1)+posMV-1 mergedCands(i,2) maxValue];  % ThetaEvts : [start peak stop PeakPower]
	end
end

duration = time(ThetaEvts(:,3))-time(ThetaEvts(:,1));                                               % Compute the duration of the periods   
ThetaEvts(duration<durations(2)/1000,:) = [];                                           % Discard periods that are shorter than minLength

tTheta=zeros(1,numel(fRes(:,2)));
for t=1:size(ThetaEvts,1)
    tTheta(ThetaEvts(t,1):ThetaEvts(t,3))=1;
end
Neurons.Theta(n).Matrix=tTheta;
Neurons.Theta(n).Starts= ThetaEvts(:,1)';
Neurons.Theta(n).Ends= ThetaEvts(:,3)';

end

% Evts=fRes(:,2); Evts(~tTheta)=nan;
% [phase,amplitude,unwrapped] = Phase(fRes);
% instantFreq= diff(unwrapped(:,2))*100/2/pi;
% figure
% plot(results.C_raw(n,:))
% hold on
% plot(fRes(:,2),'LineWidth', 2.5)
% plot(Evts,'LineWidth', 2.5)
