function SCEs=GetSCEsTimings(SCEs,Neurons,Experiment,Params,Direction)

SCEs.(Direction).StillStarts=zeros(numel(Experiment.PCs.(Direction)),Experiment.SizeData);  % SCEs Starts
PCrank=0;                                                                                   % PC's ranks
for PC = Experiment.PCs.(Direction)                                                         % Loop on PCs
   PCrank=PCrank+1;                                                                             % increase rank by one
   for event = 1:numel(Neurons.Rises(PC).Starts)                                                % Loop on PCs events
      SCEs.(Direction).StillStarts(PCrank,Neurons.Rises(PC).Starts(event))=1;                   % Keep event's start timing 
   end
end
SCEs.(Direction).StillStarts(:,~logical(Experiment.Still.Segments))=0;                      % All starts outside stillness are deleted

SCEs.(Direction).NbStarts=movsum(sum(SCEs.(Direction).StillStarts),[0 Params.SCEs.Window]);                     % Moving sum : Number of PCs whith an event start within the given window
SCEs.(Direction).MeanNbStarts=mean(SCEs.(Direction).NbStarts(:,logical(Experiment.Still.Segments)));            % Average start rate

dNbStartsOverThr=[0 diff(SCEs.(Direction).NbStarts>=round(Params.SCEs.Thr*SCEs.(Direction).MeanNbStarts))];     % Get periods where starts rate is over a threshold

HighRatePeriodsStarts=find(dNbStartsOverThr==1);                                                                % Get periods starts
HighRatePeriodsEnds=find(dNbStartsOverThr==-1);                                                                 % Get periods ends
if HighRatePeriodsStarts(1)>HighRatePeriodsEnds(1); HighRatePeriodsStarts=[1 HighRatePeriodsStarts]; end                    % Boundaries checks
if HighRatePeriodsStarts(end)>HighRatePeriodsEnds(end); HighRatePeriodsEnds=[HighRatePeriodsEnds Experiment.SizeData]; end

for start = 1: numel(HighRatePeriodsStarts)                                                                                 % Loop on periods
    roi=sum(SCEs.(Direction).StillStarts(:,HighRatePeriodsStarts(start):HighRatePeriodsEnds(start)+Params.SCEs.Window));        % get the window of putative SCE
    tmpStart=HighRatePeriodsStarts(start)+find(roi,1,'first')-1;                                                                % Find first start 
    HighRatePeriodsEnds(start)=HighRatePeriodsStarts(start)+find(roi,1,'last')-1;                                               % Find last start -> end of the period
    HighRatePeriodsStarts(start)=tmpStart;                                                                                      % First start -> start of the period 
end

HighRatePeriods=[HighRatePeriodsStarts; HighRatePeriodsEnds];                               % Merge starts/ends matrices                                
Interval=HighRatePeriods(1,2:end) - HighRatePeriods(2,1:end-1);                             % Duration of each periods
% sizeMerge=HighRatePeriodsEnds(2:end) - HighRatePeriodsStarts(1:end-1);
toMerge = Interval<=Params.SCEs.maxInterval ;%& sizeMerge<=Params.SCEs.maxSize;             % Find close periods (closer than maxInterval)
while sum(toMerge)>0                                                                        % Loop while there is still close periods to merge
    first = strfind([0 toMerge],[0 1])';                                                        % Find first period index
    second = first+1;                                                                           % Find second period index (the next one)
    HighRatePeriods(2,first) = HighRatePeriods(2,second);                                       % Merge them
    HighRatePeriods(:,second) = [];                                                             % Delete the second one
    Interval=HighRatePeriods(1,2:end) - HighRatePeriods(2,1:end-1);                             % Get new durations
    toMerge = Interval<=Params.SCEs.maxInterval ;%& sizeMerge<=Params.SCEs.maxSize;             % Check for periods to merge
end

SCEs.(Direction).Starts=HighRatePeriods(1,:);                                               % Keep periods timings
SCEs.(Direction).Ends=HighRatePeriods(2,:);
SCEs.(Direction).NbSCE = numel(SCEs.(Direction).Starts);

for start = 1:SCEs.(Direction).NbSCE
    [SCEs.(Direction).PCsRanks{start},ActTime]=find(SCEs.(Direction).StillStarts(:,SCEs.(Direction).Starts(start):SCEs.(Direction).Ends(start)));
    SCEs.(Direction).ActTime{start}=SCEs.(Direction).Starts(start)+ActTime-1;
    SCEs.(Direction).CellsRanks{start}=Experiment.PCs.(Direction)(SCEs.(Direction).PCsRanks{start});
end
