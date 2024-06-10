% Phase reconstruction
%%
chanLFP=4;
LFP=loadChanLFP(chanLFP,Experiment.path);                       % Load LFP
FiltLFP = FilterLFP(LFP,'passband','theta');
[ThetaPhase,~,uw] = Phase(FiltLFP);
uwCaThetaPhase=uw(TTLstartsLFP,:);
%%
Positions=Experiment.positionXcmSmooth;
Positions=Positions./(2*max(Positions));
Positions=[(1:numel(Positions))'/100 Positions]; 
Positions(logical(Experiment.Dir2.Segments),2)=1-Positions(logical(Experiment.Dir2.Segments),2);

uwCaThetaPhase(:,1)=Positions(:,1);
Positions(~logical(Experiment.Dir2.Segments+Experiment.Dir1.Segments),:)=[];
uwCaThetaPhase(~logical(Experiment.Dir2.Segments+Experiment.Dir1.Segments),:)=[];

% Theta SPikes
thr=0; 

sampleAct=[(0:0.01:(Experiment.SizeData-1)*0.01)' results.C_raw'];                 % Data to FMA Sample
FiltActTheta = FilterLFP(sampleAct,'passband',12.5*Params.Theta.boundsTheta);
DiscreteAct=zeros(Experiment.nNeurons,Experiment.SizeData);
disp('Neuron 001')

for neuron=1:Experiment.nNeurons
    fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',neuron)));
    Minima = islocalmin(FiltActTheta(:,neuron+1),'MinSeparation',5); Minima(~Neurons.Rises(neuron).Matrix)=0;
    DiscreteAct(neuron,:)=Minima;
%     if thr
%         [wt,f]=cwt(results.C_raw(neuron,:),100); wt=abs(wt);
%         thetaBand=find(and(f>4,f<10)); normBand=find(and(f>10,f<20));
%         rapp=sum(wt(thetaBand,:))./sum(wt(normBand,:));
%         threshold=prctile(rapp,90);
%         DiscreteAct(neuron,:)=DiscreteAct(neuron,:).*(rapp>threshold);
%     end
end

DiscreteAct(:,~logical(Experiment.Dir2.Segments+Experiment.Dir1.Segments))=0;

[CellIDs,Timestamps]=ind2sub(size(DiscreteAct),find(DiscreteAct));       % Convert to Time/Cell format
SampleDiscreteAct = [Timestamps/100,CellIDs];
%%
nBins=200;

[~,lambda,Px] = ReconstructPosition(Positions,SampleDiscreteAct,uwCaThetaPhase,'nBins',nBins,'mode','train');
[stats,~,~] = ReconstructPosition(Positions,SampleDiscreteAct,uwCaThetaPhase,'nBins',nBins,'mode','test','lambda',lambda,'Px',Px);

toDel=find(ismember(stats.estimations',(1/nBins)*ones(1,nBins) ,'rows'));

truePos=stats.positions(:,2);truePos(toDel)=nan;
[~,reconstructPos]=max(stats.estimations);reconstructPos(toDel)=nan;

error=abs(truePos-reconstructPos')/nBins;
med=median(error,'omitnan');
mea=mean(error,'omitnan');
%%
figure
histogram(error*100,100)
title('Errors distribution ; median in red ; mean in yellow')
xlabel('Error (cm)')
yl=ylim;
line(100*[med med],[yl(1) yl(2)],'color','r','Linewidth',2)
line(100*[mea mea],[yl(1) yl(2)],'color','y','Linewidth',2)
%%
stillness=logical(ismember(stats.estimations',1/nBins*ones(1,nBins),'rows'));
[~,estims]=max(stats.estimations,[],1);
estims(stillness)=nan;
pos=stats.positions(:,2);pos(stillness)=nan;
figure
%     subplot(2,1,1);
        imagesc(stats.estimations)
        hold on
        plot(pos,'r','Linewidth',2)
%         scatter(1:numel(stats.positions(:,2)),estims,'y','filled')
        plot(estims,'y','Linewidth',2)
        title(['Reconstruction : real position in red, inferred position in yellow' newline 'median error : ' num2str(round(med*100,1)) 'cm ; mean error : ' num2str(round(mea*100,1)) 'cm'])
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InstantThetaFreq=diff(smooth(uw(:,2),100))./(diff(uw(:,1))*2*pi);
CaInstantThetaFreq=InstantThetaFreq(TTLstartsLFP,:);

CaThetaPhase=ThetaPhase(TTLstartsLFP,:);
CaFiltLFP=FiltLFP(TTLstartsLFP,:);

[ActPhase,ampThetaNeur,uwActPhase] = Phase(FiltActTheta);
InstantActFreq=zeros(Experiment.SizeData-1,Experiment.nNeurons);
for n =1:Experiment.nNeurons
    InstantActFreq(:,n)=diff(smooth(uwActPhase(:,n+1),100))./(diff(uwActPhase(:,1))*2*pi);
end

%%
phaseRange=[];
amps=[];
meaPhaseDiff=[]; 
candidates=[];
for neuron=1:Experiment.nNeurons
    thrAmp=prctile(ampThetaNeur(:,neuron+1),75);
    for evt = 1:numel(Neurons.Rises(neuron).Starts)
        if sum(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))~=0
            meaDiff=mean(CaInstantThetaFreq(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))-mean(InstantActFreq(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));
            meaPhaseDiff=[meaPhaseDiff meaDiff];
            amp=median(ampThetaNeur(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),2));
            amps=[amps amp];
            nbSpikes=sum(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));
            if and(amp>thrAmp,meaDiff<0) & nbSpikes>2
                candidates=[candidates [neuron;evt]];
            end
        end
    end
end
figure
histogram(meaPhaseDiff)
figure
histogram(amps)
figure
scatter(amps,meaPhaseDiff)
xlabel('Amplitude of Event intrinsic theta oscillation')
ylabel('Frequency Difference between LFP theta and Event Theta (Hz)')
title('Relation between Theta Frequency Difference (LFP - intrinsic) and Intrinsic theta Amplitude')
%%
todel=[];
for cand=1:numel(candidates(1,:))
    neuron=candidates(1,cand);
    evt=candidates(2,cand);
    NS=sum(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));
    precess=0;
    s=find(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));
    spikesPhase=CaThetaPhase(Neurons.Rises(neuron).Starts(evt)+s-1,2);
    for shift=0:NS
        if issorted(circshift(spikesPhase,shift),'descend')
            precess=1;
        end
    end
    if ~precess || Neurons.Rises(neuron).EvtsCat(evt)==1
       todel=[todel cand]; 
    end
end
candidates(:,todel)=[];
%%
cand=10;
neuron=candidates(1,cand);
evt=candidates(2,cand);
m=max(CaFiltLFP(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),2));
s=find(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));
figure
    subplot(1,2,1)
        plot(CaInstantThetaFreq(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))
            hold on
        plot(10*CaFiltLFP(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),2)./m)
        plot(DiscreteAct(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))
        plot(InstantActFreq(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),neuron+1))
        scatter(s,10*CaFiltLFP((Neurons.Rises(neuron).Starts(evt)+s-1),2)./m,'filled')
        scatter(s,CaThetaPhase((Neurons.Rises(neuron).Starts(evt)+s-1),2),'filled')
    %       plot(FiltActTheta(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),neuron+1))
    %       plot(ActPhase(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),neuron+1))
    %       plot(uwActPhase(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),neuron+1))
          plot(CaThetaPhase(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),2))
    subplot(1,2,2)
        plot(results.C_raw(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))
        hold on
        scatter(s,results.C_raw(neuron,(Neurons.Rises(neuron).Starts(evt)+s-1)),'filled')
        plot(FiltActTheta(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt),neuron+1))