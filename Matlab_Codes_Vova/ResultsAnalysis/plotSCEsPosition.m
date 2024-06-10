function plotSCEsPosition(SCEs,Experiment,Params)

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

    toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
    toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
    
    plot(Experiment.positionXcmSmooth); plot(toPlot1); plot(toPlot2)

    scatter(SCEs.Dir1.Starts,Experiment.positionXcmSmooth(SCEs.Dir1.Starts),100,'p','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','k')%[0.9290 0.6940 0.1250])
    scatter(SCEs.Dir2.Starts,Experiment.positionXcmSmooth(SCEs.Dir2.Starts),100,'p','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor','k')%[0.9290 0.6940 0.1250])
    title('Position and Ripples Events Candidates'); xlabel('Time (#Frame)'); ylabel('Position (cm')
    
if Params.SCEs.PrintPDFs
    print(gcf,'-dpdf','-bestfit', [Params.PDFsFolder 'tmp\000']);                                   % Print pdf
    close
end