function Experiment=GetTimeSpent(Experiment,Direction)

Experiment.(Direction).TimeSpent = zeros(Experiment.(Direction).NbPassage,ceil(max(Experiment.positionXcmSmooth)));        % Time spend in each spatial bin for each mouvement ; #mouvement, #number of data point in each spatial bin

for nSeg = 1:Experiment.(Direction).NbPassage                                                                                                                                      % For each mouvement period
    for time = Experiment.(Direction).Starts(nSeg):Experiment.(Direction).Ends(nSeg)                                                                                              % Time period is from start to end of the mouvement period
        Experiment.(Direction).TimeSpent(nSeg,ceil(Experiment.positionXcmSmooth(time)))=Experiment.(Direction).TimeSpent(nSeg,ceil(Experiment.positionXcmSmooth(time)))+1;  % Increment
    end
end