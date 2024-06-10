
%%
Dir1=[Experiment.Dir1.Starts Experiment.Dir1.Ends];
Dir2=[Experiment.Dir2.Starts Experiment.Dir2.Ends];

spikes=[];
for n=1:Experiment.nNeurons
   spikes=[spikes ;Neurons.Rises(n).Starts' n*ones(numel(Neurons.Rises(n).Starts),1)];
end

[spikesMvt,~,~] = InIntervals(spikes(:,1),sort([Dir1 ; Dir2])); %on elimine les spikes pendant les pauses
spikes=spikes(spikesMvt,:);

pos=[(1:Experiment.SizeData)' 0.5*fillmissing((Experiment.positionXcmSmooth-min(Experiment.positionXcmSmooth))/max(Experiment.positionXcmSmooth),'previous')]; % on doit normaliser la position entre 0 et 1 (de 0 � 0.5 on a du min au max de dir1 et de 0.5 � 1 le min et le max de Dir2)
[Mvtd2,~,~] = InIntervals(pos(:,1),Dir2);
pos(Mvtd2,2)=1-pos(Mvtd2,2);

[Mvt,~,~] = InIntervals(pos(:,1),sort([Dir1 ; Dir2]));
pos=pos(Mvt,:);
positions=pos;

%%
spikes(:,1)=spikes(:,1)/100;% pour convertir les num�ros de frames en s, approximatif
positions(:,1)=positions(:,1)/100; 
%%
nBinsA=200;%[25 50 100 200]; %nb bins spatiaux
wndwA=0.02;%[0.01 0.02 0.03 0.05 0.1 0.2 0.5];  % fenetre de temps pour la reconstruction 
medMat=zeros(numel(nBinsA),numel(wndwA));
meanMat=zeros(numel(nBinsA),numel(wndwA));
for nb=1:numel(nBinsA)
    nBins=nBinsA(nb);
    for w =1:numel(wndwA)
        wndw=wndwA(w);
        [stats,lambda,Px] = ReconstructPosition(positions,spikes,'mode','train','window',wndw,'nBins',[nBins 1]); % train ; lambda : champs d'activite de tous les neurones, px : r�partition des position (histogramme des positions)
        [~,estim]=max(stats.estimations);
        toKeep=sum(~ismember(stats.estimations,ones(nBins,1)/nBins))>0;
        estim=estim(toKeep)';
        pos=stats.positions(toKeep,2);
        error = min(abs(estim-pos),nBins-abs(estim-pos));
        errorCm=error*200/nBins;
        medMat(nb,w)=mean(errorCm);
        meanMat(nb,w)=median(errorCm);
    end
end
%%
figure
    subplot(1,2,1)
        imagesc(meanMat)
        title('Mean error(cm)'),xlabel('Window Size (s)'),ylabel('Spatial Bin Size (cm)'); xticks(1:numel(wndwA)); yticks(1:numel(nBinsA));
        xticklabels(cellfun(@num2str,num2cell(wndwA),'UniformOutput',false)); yticklabels(cellfun(@num2str,num2cell(200*ones(1,numel(nBinsA))./nBinsA),'UniformOutput',false))
    subplot(1,2,2)
        imagesc(medMat)
        title('Median error(cm)'),xlabel('Window Size (s)'),ylabel('Spatial Bin Size (cm)'); xticks(1:numel(wndwA)); yticks(1:numel(nBinsA));
        xticklabels(cellfun(@num2str,num2cell(wndwA),'UniformOutput',false)); yticklabels(cellfun(@num2str,num2cell(200*ones(1,numel(nBinsA))./nBinsA),'UniformOutput',false))
    
%%
[Mx,MxIdx]=max(squeeze(lambda));
PFs=squeeze(lambda);
[~,sMxIdx]=sort(MxIdx);
figure
imagesc((PFs(:,sMxIdx)./Mx(sMxIdx))')
%%
figure
imagesc(stats.estimations)
hold on
plot(stats.positions(:,2),'r','LineWidth',2)
[~,estim]=max(stats.estimations);
scatter(1:numel(stats.positions(:,2)),estim,'y','filled')

%%
figure
E=stats.estimations;
Emax=max(E);
imagesc(E./Emax)









