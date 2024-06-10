% MainPath='/mnt/cortex-data-312/tmpTheo/Data/Souris154147'; ePhysChan=8; TTLChan=18;
MainPath='/mnt/cortex-data-312/tmpTheo/Data/Souris154148'; ePhysChan=8; TTLChan=18;nbExp=[1,1,1];
% MainPath='/mnt/cortex-data-312/tmpTheo/Data/Souris133656'; ePhysChan=8; TTLChan=19;

%% Load Data

expPaths=GetExpsPaths('NbExp',nbExp,'Path',MainPath);

LoadExtractedData= questdlg('Do you want to load extracted calcium data ?','Extracted data Loading','Yes','No','Yes');
switch LoadExtractedData
    case 'Yes'
        [file,path] =uigetfile([MainPath '/*.mat'],'Select Calcium results file to load');
        ExtractedResults=load(fullfile(path,file));ExtractedResults=ExtractedResults.ExtractedResults;
        nameExtract=strsplit(file,'_');nameExtract=nameExtract{2}(1:end-4);
    case 'No'
        [ExtractedResults,~,nameExtract]=ExtractExpData(expPaths,MainPath);
end

[file,path] = uigetfile([MainPath '/*.dat'],'Select ePhys data to load'); ePhysFile = fullfile(path,file);
ePhysData=loadChanLFP(ePhysChan,[ePhysFile(1:end-3) 'lfp'],1);
TTLData=loadChanDat(TTLChan,ePhysFile); TTLData=single(TTLData);

ePhysDataFilt=FilterLFP(ePhysData,'passband',[80 250]);
[spectrogram,t,f] = MTSpectrogram(ePhysData,'range',[0 500],'window',0.05);
bands = SpectrogramBands(spectrogram,f);
ripplesPower=[t bands.ripples];
%% Extract TTL Timings

TTLstartsTimes=GetTTLtimes(TTLData,ExtractedResults);

% figure
%     yyaxis left
%         plot(TTLstartsTimes,ExtractedResults.C_raw(1,:))
%         hold on
%         plot(TTLstartsTimes,ExtractedResults.C_raw(2,:))
%     yyaxis right
%         plot(ePhysData(:,1),ePhysData(:,2))
%% Get Starts

Experiment.SizeData = size(ExtractedResults.C, 2);                      % Size of Fluo Data (nb of frames)
Experiment.nNeurons = size(ExtractedResults.C, 1);   
Params.Method.ActivityType = 'Spikes';
% Params.Method.SmallEvents = 'RisesDerivative'; 
Params.Method.SmallEvents ='RisesSE';

Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),'Ends',cell(Experiment.nNeurons,1),...
        'Duration',cell(Experiment.nNeurons,1)));
    
Neurons=GetRises(ExtractedResults,Params,Experiment,Neurons);
EvtMat=vertcat(Neurons.Rises.Matrix);EvtMat=full(EvtMat)>0;
EvtStarts=[zeros(Experiment.nNeurons,1) diff(EvtMat,[],2)]>0;
%% Chance level of synchrony
% 100 perms : 3std = 5.5 ; 95th prctile = 7 ; 99th prctile = 8
% 1000 perms : 3std = 5.5 ; 95th prctile = 7 ; 99th prctile = 8
nNeurons = size(ExtractedResults.C, 1);   
window=5;

nbPerm=100;
AllNbStarts=zeros(nbPerm,size(EvtStarts,2)-window);
for p=1:nbPerm
    shifts=randi(size(EvtStarts,2),1,size(EvtStarts,1)) ;
    permEvt=cell2mat(arrayfun(@(x) circshift(EvtStarts(x,:),[1 shifts(x)]),(1:numel(shifts))','un',0));
    for t=1:size(permEvt,2)-window
        actNeurons=sum(permEvt(1:nNeurons,t:t+window-1),2)>0;
        AllNbStarts(p,t)=sum(actNeurons);
    end
end

nbStartsThr=prctile(AllNbStarts(:),99);
%% Get SCEs
% zThreshold=2;
% nbStarts=movsum(sum(EvtStarts),window); 
NeuronIDs=cell(1,size(EvtStarts,2)-window);
nbStarts=zeros(1,size(EvtStarts,2)-window);
for t=1:size(EvtStarts,2)-window
    actNeurons=sum(EvtStarts(1:nNeurons,t:t+window-1),2)>0;
    nbStarts(t)=sum(actNeurons);
    NeuronIDs{t}=find(actNeurons);
end

% nbStartsZscored=(nbStarts-mean(nbStarts))/std(nbStarts);
% putSCE=[nbStartsZscored(1)>zThreshold diff(nbStartsZscored>zThreshold)>0];
% putSCE=and(nbStartsZscored>2,islocalmax(nbStartsZscored)); --> good idea but problem when multi peaks
% putSCE=[nbStarts(1)>=nbStartsThr diff(nbStarts>=nbStartsThr)>0];
putSCE=and(nbStarts>=nbStartsThr,islocalmax(nbStarts,'MinSeparation',10));
putSCEtimes=TTLstartsTimes(logical(putSCE));
putSCEidx=find(putSCE);

for sce = 1:numel(putSCEidx)
    putSCEidx(sce)=putSCEidx(sce)+find(sum(EvtStarts(1:nNeurons,putSCEidx(sce):putSCEidx(sce)+window-1)),1)-1;
end
%% PETD of ripple-band power
[RipPowsync,indices] = Sync(ripplesPower,putSCEtimes,'durations',[-1 1]);
[d,t,p] = SyncHist(RipPowsync,indices,'mode','dist','durations',[-1 1]);%,'bins',[0 500 100]);
    figure;
        PlotColorMap(d,1); clim([0 0.01])
        xticks(1:10:100);xticklabels(round(t(1:10:100),1)); xlabel('Time centered on SCE'); 
        yticks(1:10:100);yticklabels(round(p(1:10:100))); ylabel('Power in the Ripples band')
        title('Peri-SCEs time distribution of Ripples band power')

%% Find Ripples
% noiseChan=1;
% noisyData=loadChanLFP(noiseChan,[ePhysFile(1:end-3) 'lfp'],1);
% noisyDataFilt=FilterLFP(noisyData,'passband',[80 250]);
% [ripples,sd,bad] = FindRipples(ePhysDataFilt,'noise',noisyDataFilt);
[ripples,sd,bad] = FindRipples(ePhysDataFilt);
%% CCG of SCEs and SWRs
RipAndSCEs=[[ripples(:,1) ones(size(ripples,1),1)];[putSCEtimes 2*ones(size(putSCEtimes,1),1)]];
[ccg,t] = CCG(RipAndSCEs(:,1),RipAndSCEs(:,2),'duration',2);
    figure; 
        bar(t,ccg(:,1,2)); %PlotCCG(t,ccg)
        xlabel('Delay (s)'); ylabel('Occurence (#)');title('CCG of Ripples x SCEs') 
 
%% Check SCEs sync with SWRs : 
% PRE48 : base 247 /30.46% ; without redundancy 234/34,77% ; 
% POST48 : base ? ; without redundancy 206 21,64% ; with chance level thr : 258/16,22% ; with local max : 271/19,52%

syncSCEs=zeros(1,numel(putSCEtimes));
 for sce=1:numel(putSCEtimes)
     prevSWR=find(ripples(:,1)<putSCEtimes(sce),1,'last');
     nextSWR=find(ripples(:,1)>putSCEtimes(sce),1,'first');
     if isempty(prevSWR)
         closestSWR=abs(putSCEtimes(sce)-ripples(nextSWR,1));
     elseif isempty(nextSWR)
         closestSWR=abs(putSCEtimes(sce)-ripples(prevSWR,1));
     else
         closestSWR=min(abs(putSCEtimes(sce)-ripples(prevSWR,1)),abs(putSCEtimes(sce)-ripples(nextSWR,1)));
     end
     
     if closestSWR<0.1
         syncSCEs(sce)=1;
     end
 end
 disp(sum(syncSCEs))
 disp(sum(syncSCEs)/numel(syncSCEs))
%% Check SWRs sync with SCEs : 
% PRE48 : base ?? ; without redundancy 214/55,58%
% POST48 : base ? ; without redundancy 189 40,73% ; with chance level thr :% 231/49,78% ; with local max : 267/57,54%
 syncSWRs=zeros(1,size(ripples,1));
 for swr=1:size(ripples,1)
     prevSCE=find(ripples(swr,1)>putSCEtimes(:),1,'last');
     nextSCE=find(ripples(swr,1)<putSCEtimes(:),1,'first');
     if isempty(prevSCE)
         closestSCE=abs(putSCEtimes(nextSCE)-ripples(swr,1));
     elseif isempty(nextSCE)
         closestSCE=abs(putSCEtimes(prevSCE)-ripples(swr,1));
     else
         closestSCE=min(abs(putSCEtimes(nextSCE)-ripples(swr,1)),abs(putSCEtimes(prevSCE)-ripples(swr,1)));
     end
     
     if closestSCE<0.1
         syncSWRs(swr)=1;
     end
 end
 disp(sum(syncSWRs))
 disp(sum(syncSWRs)/numel(syncSWRs))
%% Plot SWR+SCE example 
w=20;
wndw=10;
for sce = 8:9

sceIdx=bigSCEidx(sce);
start=TTLstartsTimes(sceIdx-w);
ending=TTLstartsTimes(sceIdx+wndw-1+w);
figure
yyaxis right
[n,t]=find(EvtStarts(:,sceIdx:sceIdx+wndw-1));
scatter(TTLstartsTimes(sceIdx+t-1),n)
yyaxis left
[status,~,~] = InIntervals(ePhysData(:,1),[start ending]);
plot(ePhysData(status,1),ePhysData(status,2))
ylim([-10000 10000])
end
%% Histogram of neuronal occurences in SCEs
sceNeuronIDs=NeuronIDs(putSCEidx);
allSCEneurons=cell2mat(sceNeuronIDs');
nIDsBins=0.5:1:size(EvtMat,1)+0.5;
syncSCEneurons=cell2mat(sceNeuronIDs(logical(syncSCEs))');
figure
histogram(allSCEneurons,nIDsBins,'Normalization','probability')
hold on
histogram(syncSCEneurons,nIDsBins,'Normalization','probability')

figure
scatter(1:size(EvtMat,1),histcounts(syncSCEneurons,nIDsBins)./histcounts(allSCEneurons,nIDsBins))
lsline

coAct=zeros(size(EvtMat,1),size(EvtMat,1));
for sce = 1:numel(putSCEtimes)
    neurons=sceNeuronIDs{sce};
    for n=1:numel(neurons)
        coAct(neurons(n),neurons(n+1:numel(neurons)))=coAct(neurons(n),neurons(n+1:numel(neurons)))+1;
    end
end

%% Histograms of Peak Power and neuron involment in sync SWRs and SCEs
figure
histogram(ripples(logical(syncSWRs),4),'Normalization','probability')
hold on
histogram(ripples(~logical(syncSWRs),4),'Normalization','probability')
xlabel('Peak Power of the SWR')
ylabel('Probability')
legend({'SWRs sync. with SCEs','SWRs not sync. with SCEs'})

figure
histogram(nbStarts(putSCEidx(logical(syncSCEs))),'Normalization','probability')
hold on
histogram(nbStarts(putSCEidx(~logical(syncSCEs))),'Normalization','probability')
xlabel('Nb of neurons involved in the SCE')
ylabel('Probability')
legend({'SCEs sync. with SWRs','SCEs not sync. with SWRS'})
%% Plot nb starts SCEs
sce=15;
win=25;

figure
imagesc(EvtStarts(:,putSCEidx(sce)-win:putSCEidx(sce)+win))
set(gca,'YDir','normal') 
hold on
plot(nbStarts(putSCEidx(sce)-win:putSCEidx(sce)+win),'r','Linewidth',3)
line([win win],[0 200],'Color','y')

%%
minNbNeuronsSCE=20;
minCommNeurons=8;

bigSCEs=cellfun(@numel,sceNeuronIDs)>=minNbNeuronsSCE;
bigSCEidx=putSCEidx(bigSCEs);
bigSCEnIDs=sceNeuronIDs(bigSCEs);
bigSCEsync=syncSCEs(bigSCEs);
% [Ncomm,rank1,rank2]=intersect(sce1(:,1), sce2(:,1),'stable');

MatCorCoef=zeros(numel(bigSCEnIDs));
MatCorP=ones(numel(bigSCEnIDs));
MatCorN=cell(numel(bigSCEnIDs));
for sce1=1:numel(bigSCEnIDs)-1
    for sce2=sce1+1:numel(bigSCEnIDs)
        [MatCorN{sce1,sce2},rank1,rank2]=intersect(bigSCEnIDs{sce1}, bigSCEnIDs{sce2},'stable');
        if numel(rank1)>=minCommNeurons
            [MatCorCoef(sce1,sce2),MatCorP(sce1,sce2)]=corr(rank1,rank2,'Type','Spearman');
        end
    end
end
MatCorCoef(MatCorP>0.05)=0;
imagesc(MatCorCoef)
%%
[sces1,sces2]=find(MatCorCoef);
MatCorNtest=MatCorN;
MatCorNtest(MatCorP>0.05)=[];

%%
test=MatCorCoef+triu(MatCorCoef)';
test(sum(test,2)==0,:)=[];
test(:,sum(test,1)==0)=[];
imagesc(test)
r = symrcm(test);
imagesc(test(r,r))

%%
test=bigSCEsync(sces1)+bigSCEsync(sces2);
nPerm=1000;
testP=zeros(nPerm,3);
for p =1:nPerm
    tmp1=randi(numel(bigSCEsync),numel(sces1),1);
    tmp2=randi(numel(bigSCEsync),numel(sces1),1);
    tmp=bigSCEsync(tmp1)+bigSCEsync(tmp2);
    testP(p,:)=[sum(tmp==0) sum(tmp==1) sum(tmp==2)];
end

figure    
    bar([0 1 2],[prctile(testP(:,1),95) prctile(testP(:,2),95) prctile(testP(:,3),95)])
    hold on
    histogram(test)
    bar([0 1 2],[prctile(testP(:,1),5) prctile(testP(:,2),5) prctile(testP(:,3),5)]) 

%%
load('Px_148.mat')
load('Lambda_148.mat') 
spikes=[];
for n=1:Experiment.nNeurons
   spikes=[spikes ;Neurons.Rises(n).Starts' n*ones(numel(Neurons.Rises(n).Starts),1)];
end
spikes(:,1)=TTLstartsTimes(spikes(:,1));
nBins=200;
wndw=0.02;
[stats,~,~] = ReconstructPosition([],spikes,'mode','test','window',wndw,'nBins',[nBins 1],'lambda',lambda,'Px',Px);
%%

figure
E=stats.estimations;
Emax=max(E);
E=E./Emax;
imagesc(E)

%%
test=putSCEtimes(logical(bigSCEsync));
test2=putSCEidx(logical(bigSCEsync));

sce=7;
w=10;
start=find(stats.windows(:,1)<test(sce),1,'last');
figure
imagesc(E(:,start-w:start+w))

figure
imagesc(EvtStarts(:,test2(sce)-2*w:test2(sce)+2*w))
set(gca,'YDir','normal') 
hold on
plot(nbStarts(test2(sce)-2*w:test2(sce)+2*w),'r','Linewidth',3)
line(2*[w w],[0 200],'Color','y')

%% Get centers of mass
roifn = single(full(ExtractedResults.A));
pixh = size(ExtractedResults.Cn, 1);
pixw = size(ExtractedResults.Cn, 2);
roifn = reshape(roifn, pixh, pixw,size(roifn,2));
roifnMask = roifn./max(max(roifn));
displayImg = ExtractedResults.PNR;
roifnMask_bin = roifnMask > 0.5;
CoM=zeros(nNeurons,2);
for n=1:nNeurons
    [x,y]=find(squeeze(roifnMask_bin(:,:,n))>0);
    CoM(n,:)=[mean(x) mean(y)];
end
%%
CoMsce=zeros(numel(sceNeuronIDs),2);
dist=zeros(numel(sceNeuronIDs),1);
for sce = 1:numel(sceNeuronIDs)
    CoMsce(sce,:)=[mean(CoM(sceNeuronIDs{sce},1)) mean(CoM(sceNeuronIDs{sce},2))];
    dist(sce)=median(sqrt((CoM(sceNeuronIDs{sce},1)-CoMsce(sce,1)).^2+(CoM(sceNeuronIDs{sce},2)-CoMsce(sce,2)).^2));
end

%%
nShuf=100;
distShuf=zeros(numel(sceNeuronIDs),nShuf);
for shuf=1:nShuf
    for sce = 1:numel(sceNeuronIDs)
        neuronsShuf=randperm(nNeurons,numel(sceNeuronIDs{sce}));
        CoMsce(sce,:)=[mean(CoM(neuronsShuf,1)) mean(CoM(neuronsShuf,2))];
        distShuf(sce,shuf)=median(sqrt((CoM(neuronsShuf,1)-CoMsce(sce,1)).^2+(CoM(neuronsShuf,2)-CoMsce(sce,2)).^2));
    end
end














