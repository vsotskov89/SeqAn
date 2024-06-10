
roifn = single(full(results.A));
pixh = size(results.Cn, 1);
pixw = size(results.Cn, 2);
roifn = reshape(roifn, pixh, pixw,size(roifn,2));
roifnMask = roifn./max(max(roifn));
displayImg = results.PNR;
roifnMask_bin = roifnMask > 0.5;

figure;
subplot(1,2,1)
    imagesc(displayImg)
    xticks([]);xticklabels([])
    yticks([]);yticklabels([])
    colormap gray
    title('Direction 1')
    hold on
    
    N=100;
    col=jet(N);
    for neuron=find(cell2mat({Neurons.isPC.Dir1}))
        contour(roifnMask_bin(:,:,neuron),[1,1],'LineColor',col(round(mean(Neurons.PFs(neuron).Dir1)/numel(Neurons.SummedPlaceActivity(1).Dir1)*N),:), 'linewidth', 2);
    end

subplot(1,2,2)
    imagesc(displayImg)
    xticks([]);xticklabels([])
    yticks([]);yticklabels([])
    colormap gray
    title('Direction 2')
    hold on
    
    for neuron=find(cell2mat({Neurons.isPC.Dir2}))
        contour(roifnMask_bin(:,:,neuron),[1,1],'LineColor',col(round(mean(Neurons.PFs(neuron).Dir2)/numel(Neurons.SummedPlaceActivity(1).Dir2)*N),:), 'linewidth', 2);
    end
    
     sgtitle([strrep(Experiment.file,'_','-') newline 'All PCs ROIs : PF position coded by color' newline 'Blue : 0cm --> Red : 100cm'])