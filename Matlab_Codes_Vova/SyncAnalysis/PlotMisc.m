function PlotMisc(nbStarts,sceNeuronIDs,putSCEidx,syncSCEs,syncSWRs,ripples,nNeurons,tCCG,RipSceCCG,nbSyncSCEs,nbPermSyncSCEs,thrSCEs,nbSyncSWRs,nbPermSyncSWRs,thrSWRs,permCCG)
%% Plot CCG
figure; 
    bar(tCCG,RipSceCCG(:,1,2)); %PlotCCG(t,ccg)
    xlabel('Delay (s)'); ylabel('Occurence (#)');title('CCG of Ripples x SCEs') 
    hold on
    plot(tCCG,prctile(permCCG,95),'color','r','LineWidth',2)
%% Histogram of neuronal occurences in SCEs

% allSCEneurons=cell2mat(sceNeuronIDs');
% nIDsBins=0.5:nNeurons+0.5;
% syncSCEneurons=cell2mat(sceNeuronIDs(logical(syncSCEs))');
% figure
% histogram(allSCEneurons,nIDsBins,'Normalization','probability')
% hold on
% histogram(syncSCEneurons,nIDsBins,'Normalization','probability')

%% Scatte of probability of involvment in synchronous vs async. SCEs
% figure
% scatter(1:nNeurons,histcounts(syncSCEneurons,nIDsBins)./histcounts(allSCEneurons,nIDsBins))
% lsline

%% Histograms of Peak Power and neuron involment in sync SWRs and SCEs
figure
histogram(ripples(logical(syncSWRs),4),'Normalization','probability','BinWidth',1)
hold on
histogram(ripples(~logical(syncSWRs),4),'Normalization','probability','BinWidth',1)
xlabel('Peak Power of the SWR')
ylabel('Probability')
legend({'SWRs sync. with SCEs','SWRs not sync. with SCEs'})

figure
histogram(nbStarts(putSCEidx(logical(syncSCEs))),'Normalization','probability','BinWidth',1)
hold on
histogram(nbStarts(putSCEidx(~logical(syncSCEs))),'Normalization','probability','BinWidth',1)
xlabel('Nb of neurons involved in the SCE')
ylabel('Probability')
legend({'SCEs sync. with SWRs','SCEs not sync. with SWRS'})

%% Histogram of statistical signifiance of the SWRs/SCEs synchrony

figure
    subplot(1,2,1)
        histogram(nbPermSyncSCEs, 30, 'Normalization','probability')
        h = gca;
        ylim(h.YLim); ylim('manual')
        line([thrSCEs thrSCEs],h.YLim ,'color','k','LineStyle','--','LineWidth',1.5)
        if nbSyncSCEs>thrSCEs; line([nbSyncSCEs nbSyncSCEs],h.YLim,'color','g','LineWidth',2.5); else; line([nbSyncSCEs nbSyncSCEs],h.YLim,'color','r','LineWidth',2.5); end
        xlabel('Nb of SCEs sync. on SWRs (#)'); ylabel('Nb of occurences (#)'); title('Signifiance of SCEs synchrony with SWRs')
        legend({'Data from random circ. shifts','99th percentile of random dist.','Observed nb of synch. SCEs'})
        set(gca,'FontSize',12)   

    subplot(1,2,2)
        histogram(nbPermSyncSWRs,30, 'Normalization','probability')
        h = gca;
        ylim(h.YLim); ylim('manual')
        line([thrSWRs thrSWRs],h.YLim,'color','k','LineStyle','--','LineWidth',1.5)
        if nbSyncSWRs>thrSWRs; line([nbSyncSWRs nbSyncSWRs],h.YLim,'color','g','LineWidth',2.5); else; line([nbSyncSWRs nbSyncSWRs],h.YLim,'color','r','LineWidth',2.5); end
        xlabel('Nb of SWRs sync. on SCEs (#)'); ylabel('Nb of occurences (#)'); title('Signifiance of SWRs synchrony with SCEs')
        legend({'Data from random circ. shifts','99th percentile of random dist.','Observed nb of synch. SWRs'})
        set(gca,'FontSize',12)   
        
        
        
        
        
        
        
        
        
        