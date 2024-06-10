function [Neurons,Experiment]=GetActivity(Neurons,Experiment,Direction,Params, results)

MvtZone = sum(Experiment.(Direction).TimeSpent>0,'omitnan')>0.4*Experiment.(Direction).NbPassage; %MvtZone represente les portions du track que la souris a visité pendant au moins 40% des passages

for neuron = 1:Experiment.nNeurons                                      % loop on neurons
    Neurons.PlaceActivity(neuron).(Direction) = zeros(Experiment.(Direction).NbPassage,ceil(max(Experiment.positionXcmSmooth)));          % Matrix of Place-related Activities ; dims : # neuron, # mouvement, Activity
    Neurons.SummedPlaceActivity(neuron).(Direction) = zeros(1,ceil(max(Experiment.positionXcmSmooth)));          % Matrix of Average place-related Activities ; dims : # neuron, Activity averaged on every mouvement period
    
    for nSeg = 1:Experiment.(Direction).NbPassage                                                                                      % Loop on mouvement periods
        for time = Experiment.(Direction).Starts(nSeg):Experiment.(Direction).Ends(nSeg)                    % Loop on times within mouvement periods
            Neurons.PlaceActivity(neuron).(Direction)(nSeg,ceil(Experiment.positionXcmSmooth(time))) = ...
                Neurons.PlaceActivity(neuron).(Direction)(nSeg,ceil(Experiment.positionXcmSmooth(time)))+Neurons.Rises(neuron).Matrix(time);       % Add activity to the corresponding neuron, mouvement period, spatial bin
%                Neurons.PlaceActivity(neuron).(Direction)(nSeg,ceil(Experiment.positionXcmSmooth(time)))+results.C_raw(neuron, time);   

        end
    end    
    Neurons.PlaceActivity(neuron).(Direction)=Neurons.PlaceActivity(neuron).(Direction)(:,MvtZone)./Experiment.(Direction).TimeSpent(:,MvtZone);    % normalize activity by time spend in each in each time bin
    for nSeg = 1:Experiment.(Direction).NbPassage 
         Neurons.PlaceActivity(neuron).(Direction)(nSeg,:) = movmean(Neurons.PlaceActivity(neuron).(Direction)(nSeg,:), Params.smoothActivitySpaceSR, 'omitnan'); % smoothing of the activity for each run
    end
    if Params.normActivitySR == 1
        maxes=max(Neurons.PlaceActivity(neuron).(Direction),[],2); maxes(maxes==0)=-1;                                  % Take maxes of each passage
    %     Neurons.PlaceActivity(neuron).(Direction)=Neurons.PlaceActivity(neuron).(Direction)./maxes;                     % Normalize activity by the maxes
        Neurons.PlaceActivity(neuron).(Direction)=Neurons.PlaceActivity(neuron).(Direction)./repmat(maxes,1,size(Neurons.PlaceActivity(neuron).(Direction),2));
    end
    %if 
    Neurons.SummedPlaceActivity(neuron).(Direction)=sum(Neurons.PlaceActivity(neuron).(Direction),'omitnan');       % compute mean
    Neurons.SummedPlaceActivity(neuron).(Direction) = movmean(Neurons.SummedPlaceActivity(neuron).(Direction),Params.smoothActivitySpace);  % smooth mean
    Neurons.SummedPlaceActivity(neuron).(Direction)(isnan(Neurons.SummedPlaceActivity(neuron).(Direction)))=0;        % remove nan values
    Neurons.NormedSummedPlaceActivity(neuron).(Direction)=Neurons.SummedPlaceActivity(neuron).(Direction)/max(Neurons.SummedPlaceActivity(neuron).(Direction));
    Experiment.(Direction).MvtZone =MvtZone;
end

