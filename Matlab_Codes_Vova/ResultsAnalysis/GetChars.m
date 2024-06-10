function Neurons=GetChars(Experiment,Neurons,Direction,Params)

for neuron = 1:Experiment.nNeurons
%     Neurons.Conc(neuron).(Direction)= ...
%         sum(Neurons.SummedPlaceActivity(neuron).(Direction)(Neurons.PFs(neuron).(Direction)))/sum(Neurons.SummedPlaceActivity(neuron).(Direction));
%     Neurons.Stab(neuron).(Direction)= ...
%         sum(sum(Neurons.PlaceActivity(neuron).(Direction)(:,Neurons.PFs(neuron).(Direction)),2,'omitnan')>0)/Experiment.(Direction).NbPassage;
%     Neurons.Peak(neuron).(Direction)= ...
%         max(Neurons.SummedPlaceActivity(neuron).(Direction)(:,Neurons.PFs(neuron).(Direction)))/Experiment.(Direction).NbPassage;
%     if isempty(Neurons.Peak(neuron).(Direction)); Neurons.Peak(neuron).(Direction)=0; end
    
    nPF = size(Neurons.PFs(neuron).(Direction),1);
    Neurons.Conc(neuron).(Direction) = cell(nPF,1); 
    Neurons.Stab(neuron).(Direction) = cell(nPF,1); 
    Neurons.Peak(neuron).(Direction) = cell(nPF,1); 
    Neurons.ActvityInPFvsOutPF(neuron).(Direction) = cell(nPF,1); 
    
    for kPF = 1 : nPF
        summedActivityCorrected = sum(Neurons.PlaceActivity(neuron).(Direction)(cell2mat(Neurons.PassagePFs(neuron).(Direction)(kPF)),:),'omitnan');
        Neurons.Conc(neuron).(Direction)(kPF)= ...
            {sum(summedActivityCorrected(cell2mat(Neurons.PFs(neuron).(Direction)(kPF))))/sum(summedActivityCorrected)};
        Neurons.Stab(neuron).(Direction)(kPF)= ...
            {sum(sum(Neurons.PlaceActivity(neuron).(Direction)(cell2mat(Neurons.PassagePFs(neuron).(Direction)(kPF)),cell2mat(Neurons.PFs(neuron).(Direction)(kPF))),2,'omitnan')>0)/numel(cell2mat(Neurons.PassagePFs(neuron).(Direction)(kPF)))};
        Neurons.Peak(neuron).(Direction)(kPF)= ...
            {max(summedActivityCorrected(:,cell2mat(Neurons.PFs(neuron).(Direction)(kPF))))/numel(cell2mat(Neurons.PassagePFs(neuron).(Direction)(kPF)))};
        if isempty(Neurons.Peak(neuron).(Direction)(kPF)); Neurons.Peak(neuron).(Direction)(kPF)={0}; end

        boundsPF = cell2mat(Neurons.PFs(neuron).(Direction)(kPF));
        activityOutsidePF = (sum(Neurons.SummedPlaceActivity(neuron).(Direction)) - sum(Neurons.SummedPlaceActivity(neuron).(Direction)(boundsPF)))/(size(Neurons.SummedPlaceActivity(neuron).(Direction),2)-size(boundsPF,2));
        Neurons.ActvityInPFvsOutPF(neuron).(Direction)(kPF) = ...
            {mean(Neurons.SummedPlaceActivity(neuron).(Direction)(boundsPF)) / activityOutsidePF};
    end    
end
