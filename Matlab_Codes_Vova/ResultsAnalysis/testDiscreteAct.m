sampleAct=[(0:0.01:(Experiment.SizeData-1)*0.01)' results.C_raw'];                 % Data to FMA Sample
FiltActTheta = FilterLFP(sampleAct,'passband',12.5*Params.Theta.boundsTheta); 

%%
nPCDir1=sum([Neurons.isPC.Dir1]);
cellDir1ID=find([Neurons.isPC.Dir1]);
DiscreteActDir1=zeros(Experiment.nNeurons,Experiment.SizeData);

for neuron=cellDir1ID
    Act=results.C_raw(neuron,:);
    thetaAct=FiltActTheta(:,neuron+1);
    thetaRises=thetaAct; thetaRises(~Neurons.Rises(neuron).Matrix)=0;
    Minima = islocalmin(thetaAct,'MinSeparation',5); Minima(~Neurons.Rises(neuron).Matrix)=0;
    DiscreteActDir1(neuron,:)=Minima;
    % MinimaPlusStarts=Minima; MinimaPlusStarts(Neurons.Rises(neuron).Starts)=1;
    % DiscreteAct(neuron,:)=MinimaPlusStarts;
end

DiscreteActDir1=DiscreteActDir1([Neurons.isPC.Dir1],:);
DiscreteActDir1(:,~Experiment.Dir1.Segments)=0;

AllMeanPADir1 = cell2mat({Neurons.SummedPlaceActivity([Neurons.isPC.Dir1]).Dir1}');
[~,MaxDir1]  = max(AllMeanPADir1,[],2); [~,sortedMaxDir1]=sort(MaxDir1);
AllMeanPADir1NormedSorted=AllMeanPADir1(sortedMaxDir1,:)./max(AllMeanPADir1(sortedMaxDir1,:),[],2);
cellDir1ID=cellDir1ID(:,sortedMaxDir1);

DiscreteActDir1Sorted=DiscreteActDir1(sortedMaxDir1,:);
imagesc(DiscreteActDir1Sorted)
imagesc(AllMeanPADir1NormedSorted)

%%








% plot(Act)
% hold on
% plot(thetaRises*50)
% plot(find(MinimaPlusStarts),thetaRises(MinimaPlusStarts)*50,'b*')
% plot(find(Minima),thetaRises(Minima)*50,'r*')
