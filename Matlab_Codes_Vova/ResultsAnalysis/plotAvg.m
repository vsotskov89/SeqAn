function plotAvg(Experiment,Neurons,Params,results)

SPAdir1=cell2mat({Neurons.SummedPlaceActivity(cell2mat({Neurons.isPC2.Dir1})).Dir1}');
for neuron = 1:size(SPAdir1, 1)
    SPAdir1(neuron, :) = movmean(SPAdir1(neuron, :), Params.smoothActivitySpace); %pas tres elegant de faire une boucle, a changer un jour...
end
[Maxes1,Maxes1Idx]  = max(SPAdir1,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); SPAdir1=SPAdir1./Maxes1;
SPAdir2=cell2mat({Neurons.SummedPlaceActivity(cell2mat({Neurons.isPC2.Dir2})).Dir2}');
for neuron = 1:size(SPAdir2, 1)
    SPAdir2(neuron, :) = movmean(SPAdir2(neuron, :), Params.smoothActivitySpace); %pas tres elegant de faire une boucle, a changer un jour...
end
[Maxes2,Maxes2Idx]  = max(SPAdir2,[],2); [~,sMAxes2Idx]=sort(Maxes2Idx); SPAdir2=SPAdir2./Maxes2;

figure
    s1=subplot(1,4,1);
        imagesc(SPAdir1(sMAxes1Idx,:))
        title(s1, ['Place Cells Dir1 : ' num2str(sum(cell2mat({Neurons.isPC2.Dir1})))])
        ylabel('Neurons #'); xlabel('Position (cm)')
    s2=subplot(1,4,2);
        imagesc(SPAdir2(sMAxes2Idx,:))
        title(s2, ['Place Cells Dir2 : ' num2str(sum(cell2mat({Neurons.isPC2.Dir2})))])
        ylabel('Neurons #'); xlabel('Position (cm)')
    s3=subplot(1,4,3:4); hold on
        pixh = size(results.Cn, 1); pixw = size(results.Cn, 2);
        roifn = single(full(results.A)); roifn = reshape(roifn, pixh, pixw,size(roifn,2));
        roifnMask = roifn./max(max(roifn)); roifnMask_bin = roifnMask > 0.5;
        imagesc(results.PNR)
        for neuron=1:Experiment.nNeurons
            if or(Neurons.isPC2(neuron).Dir1,Neurons.isPC2(neuron).Dir2)
                contour(roifnMask_bin(:,:,neuron),[1,1],'LineColor','g', 'linewidth', 1);
            else                    
                contour(roifnMask_bin(:,:,neuron),[1,1],'LineColor','r', 'linewidth', 1);
            end
        end
        xticks([]); yticks([]); title(s3,'Position of ROIs in the FOV (PCs in green)')
    sgtitle(['PCs repartition in the track' newline 'Number of neurons : ' num2str(Experiment.nNeurons) newline ...
        'Number of PCs : ' num2str(sum(cell2mat({Neurons.isPC2.Dir1}))+sum(cell2mat({Neurons.isPC2.Dir2}))-sum(and(cell2mat({Neurons.isPC2.Dir1}),cell2mat({Neurons.isPC2.Dir2})))) ...
        ' (' num2str(sum(and(cell2mat({Neurons.isPC2.Dir1}),cell2mat({Neurons.isPC2.Dir2})))) ' bidir.)'])

if Params.PrintPDFs
    print(gcf,'-dpdf','-bestfit', [Params.PDFsFolder 'tmp/000']);                                   
    close
    names=dir([Params.PDFsFolder 'tmp/*.pdf']); names =fullfile({names(:).folder}, {names(:).name});                           
    append_pdfs([Params.PDFsFolder 'AllPDFs'], names{:})
    delete([Params.PDFsFolder 'tmp/*.pdf'])
   % movefile([Params.PDFsFolder 'AllPDFs'],[Params.PDFsFolder Experiment.file '_AllCells_' Params.Method.SmallEvents '.pdf'])
    movefile([Params.PDFsFolder 'AllPDFs'],[Params.PDFsFolder '_AllCells_' Params.Method.SmallEvents '.pdf'])   
end
end
