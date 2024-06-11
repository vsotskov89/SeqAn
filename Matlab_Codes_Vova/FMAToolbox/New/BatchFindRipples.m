function ripples = BatchFindRipples(session,channels, noisyChannel)

% this is the algorithm for finding ripple events developed in May 2016 for
% Céline's data. The reasoning is to take the filtered (ripple-band) signal's
% envelope from multiple sources (carefully picked informative channels) together
% to find ripple events.
% The noise channel is used as a control, so that the envelop of this channel
% is subtracted from the mean envelop of the other channels. This way, noisy
% artefects are not detected as ripples. 

%% Initial parameters
here = max(strfind(session,'/'));
try cd([session(1:here) '/events/']); catch cd(session(1:here));end
threshold = 1;
peakThreshold = 3;

minInterRippleInterval = 0.030;
minRippleDuration = 0.020;
maxRippleDuration = 0.110;
aroundArtefact = 0.05; epsilon = 0.15; %for consolidating noisy periods
SetCurrentSession(session,'verbose','off');

try
    lfp = GetLFP(noisyChannel);
    filtered = FilterLFP(Detrend(lfp), 'passband', 'ripples');
    [~,amplitude] = Phase(filtered);
    noisy = amplitude(:,2);
catch
    %if there was not a valid noiseChannel, use the envelops for higher frequencies
    % artefacts tend to bleed to all frequencies, so this way they will be subtracted
%     noise = zeros(length(t),numel(channels));
    for i=1:length(channels),
        lfp = GetLFP(channels(i));
        filtered = FilterLFP(Detrend(lfp), 'passband', [300 500]);
        [~,amplitude] = Phase(filtered);
        noise(:,i) = amplitude(:,2);
    end
    noisy = nanmean(noise,2);
end

t = lfp(:,1);

% Get envelopes

envelope = zeros(length(t),numel(channels));

for i=1:numel(channels),
    lfp = GetLFP(channels(i));
    filtered = FilterLFP(Detrend(lfp), 'passband', 'ripples');
    [~,amplitude] = Phase(filtered);
    envelope(:,i) = amplitude(:,2);
end
%%
me = nanmean(envelope,2); %mean envelope

% Get rid of (saturated) noisy periods

badperiods = [];
if any(abs(lfp(:,2))>=32766),% saturation @ 32766;
    bad = abs(lfp(:,2))>=32766;
    starts = strfind(bad', [0 1])';
    ends = strfind(bad', [1 0])';
    if starts(1)>ends(1); starts = [1; starts]; end
    if starts(end)>ends(end); ends = [ends; length(t)]; end
    badperiods = [t(starts)-aroundArtefact t(ends)+aroundArtefact];
    % Consolidate bad periods that are separated by less than epsilon
    badperiods = ConsolidateIntervals(badperiods, 'epsilon', epsilon);
elseif str2double(session(end-15:end-13)) == 291; % Rat 291;
    lfp = Detrend(GetLFP(channels(1))); 
    bad = abs(lfp(:,2))>=500;
    starts = strfind(bad', [0 1])';
    ends = strfind(bad', [1 0])';
    if starts(1)>ends(1); starts = [1; starts]; end
    if starts(end)>ends(end); ends = [ends; length(t)]; end
    badperiods = [t(starts)-aroundArtefact t(ends)+aroundArtefact];
    % Consolidate bad periods that are separated by less than epsilon
    badperiods = ConsolidateIntervals(badperiods, 'epsilon', epsilon);
end

if ~isempty(badperiods),me(InIntervals(t,badperiods)) = nan;end

% Find ripples

signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(signal(me>noisy));
signal(~(me>noisy)) = 0;

% get events above threshold sd-s
sig = signal>threshold;
starts = strfind([0 sig' 0],[0 1])';
stops = strfind([0 sig' 0],[1 0])' - 1;
r = t([starts stops]);
if starts(1)==1, r(1,:) = []; end
if stops(end)==length(t), r(end,:) = []; end

% get the local minima around the event as the start/stop points
minima = t(FindLocalMinima(signal)); %original signal, no nans
minima = [minima(1:end-1) minima(2:end)];
if r(1)<minima(1),r(1,:) = [];end
if r(end)>minima(end),r(end,:) = [];end

starts = minima(CountInIntervals(r(:,1)+eps,minima)>0,1);
stops = minima(CountInIntervals(r(:,2)-eps,minima)>0,2);
r = [starts stops];

% merge close ones
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;

while sum(toMerge)>0,
    % Get the index of the first ripple in a pair to be merged
    rippleStart = strfind([0 toMerge'],[0 1])';
    rippleEnd = rippleStart+1; % next ripple
    % Adjust end (incorporate second ripple into first)
    r(rippleStart,2) = r(rippleEnd,2);
    % Remove now-redundant second ripple
    r(rippleEnd,:) = [];
    iri = r(2:end,1) - r(1:end-1,2);
    newSize =  r(2:end,2) - r(1:end-1,1);
    toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
end

% peak threshold check
[in,w] = InIntervals(t,r);
[peakV,peak] = accumax(w(in),signal(in));
ok = peakV>peakThreshold;
tIn = t(in);
peak = tIn(peak(ok));
r = [r(ok,1) peak r(ok,end) peakV(ok)];

% duration check
d = r(:,3) - r(:,1);
r(d<minRippleDuration | d>maxRippleDuration,:) = [];

ripples = r;
% Save ripples

% dlmwrite('ripples',r,'precision','%f8');
% SaveRippleEvents('ripples.rip.evt', r, 0,'overwrite', 'on');
