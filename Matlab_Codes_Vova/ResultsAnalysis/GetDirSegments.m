function Experiment=GetDirSegments(Params,Experiment,Direction)

if strcmp(Direction,'Dir1')
    segments = single(and(Experiment.mouseSpeed > Params.speedThreshold,Experiment.mouseSpeed<Params.speedMax));   % periods of positive (left to right) mouvement
elseif strcmp(Direction,'Dir2')
    segments = single(and(Experiment.mouseSpeed < -Params.speedThreshold,Experiment.mouseSpeed>-Params.speedMax));
elseif strcmp(Direction,'Still')
    segments=~(Experiment.Dir1.Segments+Experiment.Dir2.Segments);
end

starts =find(diff(segments)>0);                                         % Find mouvement periods' starts 
ends =find(diff(segments)<0);                                           % Find mouvement periods' ends 
if starts(1)>ends(1)                                                    % if first end is before first start it means that the mouse was moving since the beginning
    starts = [1; starts];                                               % so we add time 1 as the first start
end
if ends(end)<starts(end)                                                % if last end is before last start it means that the mouse was still moving at the end of the movie
    ends = [ends; Experiment.SizeData-Params.nStart+1];                                             % so we add the last time as the last end
end

if ~strcmp(Direction,'Still')    
    toDel=[];
    for s =1:length(starts)                                                    % Loop for deleting small mouvements
        if abs(Experiment.positionXcmSmooth(ends(s))-Experiment.positionXcmSmooth(starts(s))) < Params.lengthThreshold     % If the mouse moved less than lengthThreshold cm
            toDel=[toDel s];                                                      % get period # to delete
        end
    end
    starts(toDel)=[];                                                          % Delete small periods
    ends(toDel)=[];
    
    toDel=[];
    for s =1:length(starts)-1                                                    % Loop for forgiving small periods of immobility
        if Experiment.timeFluo(starts(s+1))-Experiment.timeFluo(ends(s)) < Params.stopThreshold     % If the mouse was too slow for less than stopThreshold
            ends(s)=ends(s+1);                                                      % consider it as a mouvement and not a real stop (not the most beautiful way to do it but hey, it works ! 
            toDel=[toDel s+1]; 
        end
    end
    starts(toDel)=[];                                                          % Delete small periods
    ends(toDel)=[];
    
    segments=zeros(size(segments,1),1);                                     
    idxDir=arrayfun(@(B,C)B:C,starts,ends,'UniformOutput',0);               % Get idx of mouvement periods
    for seg = 1 : numel(idxDir)
         segments(idxDir{seg}) = 1;                                            % Rebuild segments Mat with new periods
    end
end

Experiment.(Direction).Segments=segments;
Experiment.(Direction).Starts=starts;
Experiment.(Direction).Ends=ends;
Experiment.(Direction).NbPassage=numel(Experiment.(Direction).Starts);



