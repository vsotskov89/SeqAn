

%%    Positions   :    Sample of positions in [0..1]
Positions=Experiment.positionXcmSmooth;
Positions=Positions./(2*max(Positions));
Positions=[(1:numel(Positions))'/100 Positions]; 
Positions(logical(Experiment.Dir2.Segments),2)=1-Positions(logical(Experiment.Dir2.Segments),2);
Positions(~logical(Experiment.Dir2.Segments+Experiment.Dir1.Segments),:)=[];
%%    Spikes      :    Sample (t,CellID)
% method='Continuous';
method='thetaSpikes';
% method='Starts';
% method='Uniform';
switch method
    case 'thetaSpikes'
        sampleAct=[(0:0.01:(Experiment.SizeData-1)*0.01)' results.C_raw'];                 % Data to FMA Sample
        FiltActTheta = FilterLFP(sampleAct,'passband',12.5*Params.Theta.boundsTheta);
        DiscreteAct=zeros(Experiment.nNeurons,Experiment.SizeData);
        disp('Neuron 001')
        for neuron=1:Experiment.nNeurons
            fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',neuron)));
            Minima = islocalmin(FiltActTheta(:,neuron+1),'MinSeparation',5); Minima(~Neurons.Rises(neuron).Matrix)=0;
            DiscreteAct(neuron,:)=Minima;
        end
    case 'Continuous'
        DiscreteAct=cell2mat({Neurons.Rises.Matrix}')>0;
    case 'Starts'
        DiscreteAct=zeros(Experiment.nNeurons,Experiment.SizeData);
        disp('Neuron 001')
        for neuron=1:Experiment.nNeurons
            fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',neuron)));
            starts = Neurons.Rises(neuron).Starts;
            DiscreteAct(neuron,starts)=1;
        end
    case 'Uniform'
        freq=8; sampling=100; spacing = sampling/freq;
        DiscreteAct=zeros(Experiment.nNeurons,Experiment.SizeData);
        disp('Neuron 001')
        for neuron=1:Experiment.nNeurons
            fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',neuron)));
            for evt =1:numel(Neurons.Rises(neuron).Starts)
                spikes = round(Neurons.Rises(neuron).Starts(evt):spacing:Neurons.Rises(neuron).Ends(evt));
                DiscreteAct(neuron,spikes)=1;
            end
        end        
end
DiscreteAct(:,~logical(Experiment.Dir2.Segments+Experiment.Dir1.Segments))=0;
[CellIDs,Timestamps]=ind2sub(size(DiscreteAct),find(DiscreteAct));       % Convert to Time/Cell format
SampleDiscreteAct = [Timestamps/100,CellIDs];    

%%
nBins=200;%[100 200 300];%[100 125 150 175 200 250 275 300]; %% 300-350 good for 656 but you must not go over 200for 884 ?
Windows=0.025;%0.125;%[0.025 0.05 0.1 0.15];%[0.1 0.125 0.15 0.175 0.2 0.25 0.5]; % 0.1-0.125 med ; (0.15 or) 0.25 mean ; 0.15 for 884
% --> Optimum 350/0.125 for 656 ; 200/0.15 for 884

medians=zeros(numel(nBins),numel(Windows));
means=zeros(numel(nBins),numel(Windows));

for b = 1:numel(nBins)
    for w = 1:numel(Windows)
        disp([nBins(b) Windows(w)])
        [~,lambda,Px] = ReconstructPosition(Positions,SampleDiscreteAct,'window',Windows(w),'nBins',nBins(b),'mode','train');
        [stats,~,~] = ReconstructPosition(Positions,SampleDiscreteAct,'window',Windows(w),'nBins',nBins(b),'mode','test','lambda',lambda,'Px',Px);

        toDel=find(ismember(stats.estimations',(1/nBins(b))*ones(1,nBins(b)) ,'rows'));

        truePos=stats.positions(:,2);truePos(toDel)=nan;
        [~,reconstructPos]=max(stats.estimations);reconstructPos(toDel)=nan;

        error=abs(truePos-reconstructPos')/nBins(b);
        medians(b,w)=median(error,'omitnan');
        means(b,w)=mean(error,'omitnan');
    end
end
% figure
% imagesc(stats.estimations)
% hold on
% plot(truePos,'r','Linewidth',2)
% plot(reconstructPos,'y','Linewidth',2)


%%
stillness=logical(ismember(stats.estimations',1/nBins*ones(1,nBins),'rows'));
[~,estims]=max(stats.estimations,[],1);
estims(stillness)=nan;
pos=stats.positions(:,2);pos(stillness)=nan;
%%
figure
%     subplot(2,1,1);
        imagesc(stats.estimations)
        hold on
        plot(pos,'r','Linewidth',2)
%         scatter(1:numel(stats.positions(:,2)),estims,'y','filled')
        plot(estims,'y','Linewidth',2)
        title(['Reconstruction : real position in red, inferred position in yellow' newline 'median error : ' num2str(round(medians*100,1)) 'cm ; mean error : ' num2str(round(means*100,1)) 'cm'])
        
%     subplot(2,1,2);
%         imagesc(stats.errors)
%         hold on
%         plot(nBins/2*ones(1,numel(stats.positions(:,2))),'r','Linewidth',2)

%%
figure
histogram(error*100,100)
title('Errors distribution ; median in red ; mean in yellow')
xlabel('Error (cm)')
yl=ylim;
line(100*[median(error,'omitnan') median(error,'omitnan')],[yl(1) yl(2)],'color','r','Linewidth',2)
line(100*[mean(error,'omitnan') mean(error,'omitnan')],[yl(1) yl(2)],'color','y','Linewidth',2)

%%
figure
imagesc(stats.average)
title('Average error over a theta cycle 655 : 23/02')
xlabel('Theta Phase')
xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
xticklabels({'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'})
yticks(0:20:200)
yticklabels({'-100',  '-80',  '-60',  '-40',  '-20',    '0',   '20',   '40',   '60',   '80',  '100'})
ylabel('Average error (cm)')
