function Experiment=getSpeedMap(Experiment,Direction)
                                                % Number of passages
Experiment.(Direction).SpeedMap = zeros(Experiment.(Direction).NbPassage,ceil(max(Experiment.positionXcmSmooth)));          % Initialize Speed Matrix # mouvement periods , Positions

for nSeg = 1:Experiment.(Direction).NbPassage                                          % loop on mouvement periods 
    for time = Experiment.(Direction).Starts(nSeg):Experiment.(Direction).Ends(nSeg)                      % loop on times within mouvement periods
        Experiment.(Direction).SpeedMap(nSeg,ceil(Experiment.positionXcmSmooth(time))) = ...
            Experiment.(Direction).SpeedMap(nSeg,ceil(Experiment.positionXcmSmooth(time)))+Experiment.mouseSpeed(time);  %add speed to the corresponding mouvement period, spatial bin
    end
end    

Experiment.(Direction).SpeedMap=Experiment.(Direction).SpeedMap./squeeze(Experiment.(Direction).TimeSpent);                      % normalize speed by time spend in each in each time bin