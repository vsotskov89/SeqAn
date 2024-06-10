function Neurons=GetChars(Experiment,Neurons,Direction,Params)

for neuron = 1:Experiment.nNeurons
%     Neurons.Conc(neuron).(Direction)= ...
%         sum(Neurons.SummedPlaceActivity(neuron).(Direction)(Neurons.PFs(neuron).(Direction)))/sum(Neurons.SummedPlaceActivity(neuron).(Direction));
%     Neurons.Stab(neuron).(Direction)= ...
%         sum(sum(Neurons.PlaceActivity(neuron).(Direction)(:,Neurons.PFs(neuron).(Direction)),2,'omitnan')>0)/Experiment.(Direction).NbPassage;
%     Neurons.Peak(neuron).(Direction)= ...
%         max(Neurons.SummedPlaceActivity(neuron).(Direction)(:,Neurons.PFs(neuron).(Direction)))/Experiment.(Direction).NbPassage;
%     if isempty(Neurons.Peak(neuron).(Direction)); Neurons.Peak(neuron).(Direction)=0; end
    
    summedActivityCorrected = sum(Neurons.PlaceActivity(neuron).(Direction)(Neurons.PassagePFs(neuron).(Direction),:),'omitnan');
    Neurons.Conc(neuron).(Direction)= ...
        sum(summedActivityCorrected(Neurons.PFs(neuron).(Direction)))/sum(summedActivityCorrected);
    Neurons.Stab(neuron).(Direction)= ...
        sum(sum(Neurons.PlaceActivity(neuron).(Direction)(Neurons.PassagePFs(neuron).(Direction),Neurons.PFs(neuron).(Direction)),2,'omitnan')>0)/numel(Neurons.PassagePFs(neuron).(Direction));
    Neurons.Peak(neuron).(Direction)= ...
        max(summedActivityCorrected(:,Neurons.PFs(neuron).(Direction)))/numel(Neurons.PassagePFs(neuron).(Direction));
    if isempty(Neurons.Peak(neuron).(Direction)); Neurons.Peak(neuron).(Direction)=0; end
    filteredActPos = movmean(Neurons.SummedPlaceActivity(neuron).(Direction), Params.smoothActivitySpace);
    boundsPF = Neurons.PFs(neuron).(Direction);
    activityOutsidePF = (sum(filteredActPos) - sum(filteredActPos(boundsPF)))/(size(filteredActPos,2)-size(boundsPF,2));
    Neurons.ActvityInPFvsOutPF(neuron).(Direction) = mean(filteredActPos(boundsPF)) / activityOutsidePF;

end