function [Neurons]=GetPFs(Experiment,Neurons,Params,Direction)

Occupation=Experiment.(Direction).TimeSpent(:,Experiment.(Direction).MvtZone)>0;
for neuron = 1:Experiment.nNeurons
    %thrPF=Params.critPCs.thrPFs*mean(Neurons.SummedPlaceActivity(neuron).(Direction));
    %bounds=diff(Neurons.SummedPlaceActivity(neuron).(Direction)>thrPF);
    filteredPF = movmean(Neurons.SummedPlaceActivity(neuron).(Direction), Params.smoothActivitySpace);
    thrPF=Params.critPCs.thrPFs*mean(filteredPF);
    bounds=diff(filteredPF>thrPF);
    PFslims=find(bounds);
    if bounds(find(bounds,1,'first'))==-1; PFslims=[1 PFslims]; end
    if bounds(find(bounds,1,'last'))==1; PFslims = [PFslims sum(Experiment.(Direction).MvtZone)]; end
    PFslims=reshape(PFslims,2,[]);

   % Max=find(Neurons.SummedPlaceActivity(neuron).(Direction)==max(Neurons.SummedPlaceActivity(neuron).(Direction)),1,'last');
    Max=find(filteredPF==max(filteredPF),1,'last');
    BinMax=find(PFslims(1,:)<=Max,1,'last');
    Neurons.PFs(neuron).(Direction)=PFslims(1,BinMax):PFslims(2,BinMax);
    
    Neurons.PassagePFs(neuron).(Direction)=find(sum(Occupation(:,Neurons.PFs(neuron).(Direction)),2)>numel(Neurons.PFs(neuron).(Direction))*Params.fractPassagePF);
end


