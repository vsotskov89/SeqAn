function Experiment=GetDirSegments(Params,Experiment,Direction)

if strcmp(Direction,'Dir1')
    segments = single(and(Experiment.mouseSpeed > Params.speedThreshold,Experiment.mouseSpeed<Params.speedMax));   % periods of positive (left to right) mouvement
elseif strcmp(Direction,'Dir2')
    segments = single(and(Experiment.mouseSpeed < -Params.speedThreshold,Experiment.mouseSpeed>-Params.speedMax));
elseif strcmp(Direction,'Still')
    segments = single((abs(Experiment.mouseSpeed) < Params.speedThreshold));
end

starts =find(diff(segments)>0);                                         % Find mouvement periods' starts 
ends =find(diff(segments)<0);                                           % Find mouvement periods' ends 
if starts(1)>ends(1)                                                    % if first end is before first start it means that the mouse was moving since the beginning
    starts = [1; starts];                                               % so we add time 1 as the first start
end
if ends(end)<starts(end)                                                % if last end is before last start it means that the mouse was still moving at the end of the movie
    ends = [ends; Experiment.SizeData];                                            % so we add the last time as the last end
end

if ~strcmp(Direction,'Still')
    toDel=[];                                                                       
    for s =1:length(starts)                                                    % Loop for deleting small mouvements
        if abs(Experiment.positionXcmSmooth(ends(s))-Experiment.positionXcmSmooth(starts(s))) < 25     % If the mouse mouved less than 40cm % changed from 40 to 25
            toDel=[toDel s];                                                      % get period # to delete
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


