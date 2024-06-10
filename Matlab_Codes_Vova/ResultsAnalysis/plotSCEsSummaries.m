function plotSCEsSummaries(SCEs,Neurons,results,Experiment,Params,FiltAct,Direction)

pixh = size(results.Cn, 1); pixw = size(results.Cn, 2);
roifn = single(full(results.A)); roifn = reshape(roifn, pixh, pixw,size(roifn,2));
roifnMask = roifn./max(max(roifn)); roifnMask_bin = roifnMask > 0.5;
displayImg = results.PNR;

lengthTrack=numel(Experiment.(Direction).MvtZone);

toPlot1=Experiment.positionXcmSmooth; toPlot1(~Experiment.Dir1.Segments)=nan;
toPlot2=Experiment.positionXcmSmooth; toPlot2(~Experiment.Dir2.Segments)=nan;


for sce = 1%:SCEs.(Direction).NbSCE
fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',sce)));    
figure('Visible','on')
    nbPCs=numel(SCEs.(Direction).PCsRanks{sce});
    subplot(3,3,1)
        seqSCE=SCEs.(Direction).StillStarts(SCEs.(Direction).PCsRanks{sce},SCEs.(Direction).Starts(sce):SCEs.(Direction).Ends(sce));
        imagesc(seqSCE)
            yticks(1:nbPCs)
            yticklabels(num2cell(SCEs.(Direction).CellsRanks{sce}))
            xticks(1:2:size(seqSCE,2))
            xticklabels(cellfun(@(c)num2str(10*(str2double(c)-1)),xticklabels,'UniformOutput',0))
            ylabel('Neurons (#)')
            xlabel('Time since event start (ms)')
            title('Activation order')
     subplot(3,3,2)
        pfsSCE=cell2mat({Neurons.NormedSummedPlaceActivity(SCEs.(Direction).CellsRanks{sce}).Dir1}');
        imagesc(pfsSCE);hold on;
        scatter(Experiment.positionXcmSmooth(SCEs.(Direction).ActTime{sce})-find(Experiment.(Direction).MvtZone,1,'first')+1,1:nbPCs,200,'r','+','MarkerFaceColor','r')
            ylim([0.5 nbPCs+0.5])
            xlim([0 lengthTrack]-find(Experiment.(Direction).MvtZone,1,'first')+1)
            yticks(1:nbPCs)
            yticklabels(num2cell(SCEs.(Direction).CellsRanks{sce}))
            xticks((-find(Experiment.(Direction).MvtZone,1,'first')+1):25:lengthTrack)
            xticklabels(cellfun(@(c)num2str(str2double(c)+find(Experiment.(Direction).MvtZone,1,'first')-1),xticklabels,'UniformOutput',0))
            ylabel('Neurons (#)')
            xlabel('Position (cm)')
            title(['Corresponding PFs' newline '(red + : mouse position' newline 'at activation time)'])
    subplot(3,3,3)
        pfsAlign=(seqSCE'*pfsSCE)';
        imagesc(pfsAlign./max(pfsAlign))
            xticks(1:2:size(seqSCE,2))
            xticklabels(cellfun(@(c)num2str(10*(str2double(c)-1)),xticklabels,'UniformOutput',0))
            ylim([0 lengthTrack]-find(Experiment.(Direction).MvtZone,1,'first')+1)
            yticks((-find(Experiment.(Direction).MvtZone,1,'first')+1):25:lengthTrack)
            yticklabels(cellfun(@(c)num2str(str2double(c)+find(Experiment.(Direction).MvtZone,1,'first')-1),yticklabels,'UniformOutput',0))
            ylabel('Neurons (#)')
            xlabel('Time since event start (frame)')
            ylabel('Position (cm)')
            title("PF Alignement")
    subplot(3,3,4)
        imagesc(displayImg); hold on;
            for cell = 1 : nbPCs
                contour(roifnMask_bin(:,:,SCEs.(Direction).CellsRanks{sce}(cell)),[1,1],'LineColor','g', 'linewidth', 1);
                [x,y]=find(roifnMask_bin(:,:,SCEs.(Direction).CellsRanks{sce}(cell)));
                coor=median([x,y]);
                text(coor(2),coor(1),num2str(SCEs.(Direction).CellsRanks{sce}(cell)),'Color','y','FontSize',12,'FontWeight','bold');
            end
            xticks([])
            yticks([])
    subplot(3,3,5:6)  
        p=plot(FiltAct(SCEs.(Direction).CellsRanks{sce},SCEs.(Direction).Starts(sce)-19:SCEs.(Direction).Starts(sce)+40)');
            xticks([0 20 40 60])
            xticklabels([-200 0 200 400])
            xlabel('Time relative to SCE onset (ms)')
            ylabel('Fluorescence')
            title("Fluo traces ; Dashed red line : SCE's onset)")
            line([20 20],[min(ylim),max(ylim)],'LineWidth',1,'Color','r','LineStyle',':')
            legend(p,string(SCEs.(Direction).CellsRanks{sce}),'location','eastoutside')
    subplot(3,3,7)  % Correlations
        imagesc(Neurons.CorrMatrix(SCEs.(Direction).CellsRanks{sce},SCEs.(Direction).CellsRanks{sce}))
%             set(gca,'xaxisLocation','top')
%             x = repmat(1:numel(nIdx),numel(nIdx),1); y = x';
%             t = num2cell(round(correlCraw(nIdx,nIdx),2)); t = cellfun(@num2str, t, 'UniformOutput', false); 
% %             text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
            xticks(1:nbPCs); xticklabels(SCEs.(Direction).CellsRanks{sce}) ; xlabel('Neurons #') 
            yticks(1:nbPCs); yticklabels(SCEs.(Direction).CellsRanks{sce}) ; ylabel('Neurons #')
            clim([0.05 0.35])
        subplot(3,3,8:9) % Position of the mouse during event 
            wndw=1000;
            plot(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1,Experiment.positionXcmSmooth(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1))
            hold on
            plot(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1,toPlot1(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1))
            plot(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1,toPlot2(SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1):SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1))
            scatter(SCEs.(Direction).Starts(sce),Experiment.positionXcmSmooth(SCEs.(Direction).Starts(sce)),200,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor','k')
                xlim([SCEs.(Direction).Starts(sce)-min(wndw,SCEs.(Direction).Starts(sce)-1) SCEs.(Direction).Starts(sce)+min(wndw,numel(Experiment.positionXcmSmooth)-SCEs.(Direction).Starts(sce)-1)-1])
                ylim([-10 110])
                legend({'Position','Dir.1','Dir.2','Event'},'location','eastoutside')
                ylabel('Position (cm)')
                xlabel('Time (frame)')
                title('Event timing relative to position')
    sgtitle(['PCs with PFs in ' Direction ' , Event ' num2str(sce)])
    if Params.SCEs.PrintPDFs
        print(gcf,'-dpdf','-bestfit', [Params.PDFsFolder 'tmp\' num2str(sprintf('%03d',neuron))]);
        close
    end
end

if Params.SCEs.PrintPDFs
    names=dir([Params.PDFsFolder 'tmp/*.pdf']); names =fullfile({names(:).folder}, {names(:).name});                           
    append_pdfs([Params.PDFsFolder 'AllPDFs'], names{:})
    delete *.pdf
    movefile([Params.PDFsFolder 'AllPDFs'],[Params.PDFsFolder Experiment.file '_AllSCEs_Win' ...
        num2str(Params.SCEs.Window) 'Thr' num2str(Params.SCEs.Thr) 'maxInt' num2str(Params.SCEs.maxInterval) '.pdf'])  
end
    