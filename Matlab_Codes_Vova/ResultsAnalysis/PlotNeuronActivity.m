figure('Position', [200 120 900 600],'Visible','on');
        subplot(2,10,2:5);
            I1=imagesc(Neurons.PlaceActivity(neuron).Dir1); hold on;
                set(I1,'AlphaData',~isnan(Neurons.PlaceActivity(neuron).Dir1))
                set(gca,'YDir','normal')
                clim([0 max(prctile(Neurons.PlaceActivity(neuron).Dir1(:),99),0.8)])
            P1=plot(Neurons.SummedPlaceActivity(neuron).Dir1,'LineWidth',2,'Color','r');
            yticklabels({})
            xlabel('Position (cm)')
            if Neurons.isPC(neuron).Dir1
                line([Neurons.PFs(neuron).Dir1(1) Neurons.PFs(neuron).Dir1(1)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                line([Neurons.PFs(neuron).Dir1(end) Neurons.PFs(neuron).Dir1(end)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                title(['{\color{green} PC} for Direction 1' newline ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir1)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir1,3)) '%' newline ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir1,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir1,3)) '%'])
                subplot(2,10,1)
                    I2=imagesc(sum(Neurons.PlaceActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0);
                        set(I2,'AlphaData',~isnan(sum(Neurons.PlaceActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0))
                    xticks({})
                    xlabel('Stability')
                    ylabel('move # (Direction 1)')
                    colormap(subplot(2,10,1),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
                    set(gca,'YDir','normal')
            elseif ~isempty(Neurons.PFs(neuron).Dir1)
                line([Neurons.PFs(neuron).Dir1(1) Neurons.PFs(neuron).Dir1(1)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                line([Neurons.PFs(neuron).Dir1(end) Neurons.PFs(neuron).Dir1(end)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                title(['{\color{red} Not a PC} for Direction 1' newline ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir1)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir1,3)) '%' newline ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir1,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir1,3)) '%'])
            end
            
        subplot(2,10,7:10);
            I3=imagesc(Neurons.PlaceActivity(neuron).Dir2); hold on;
                set(I3,'AlphaData',~isnan(Neurons.PlaceActivity(neuron).Dir2))
                set(gca,'YDir','normal')
                clim([0 max(prctile(Neurons.PlaceActivity(neuron).Dir2(:),99),0.8)])
            P2=plot(Neurons.SummedPlaceActivity(neuron).Dir2,'LineWidth',2,'Color','r');
            yticklabels({})
            xlabel('Position (cm)')
            if Neurons.isPC(neuron).Dir2
                line([Neurons.PFs(neuron).Dir2(1) Neurons.PFs(neuron).Dir2(1)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                line([Neurons.PFs(neuron).Dir2(end) Neurons.PFs(neuron).Dir2(end)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                title(['{\color{green} PC} for Direction 2' newline ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir2)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir2,3)) '%' newline ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir2,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir2,3)) '%'])
                subplot(2,10,6)
                    I4=imagesc(sum(Neurons.PlaceActivity(neuron).Dir2(:,Neurons.PFs(neuron).Dir2),2,'omitnan')>0);
                        set(I4,'AlphaData',~isnan(sum(Neurons.PlaceActivity(neuron).Dir2(:,Neurons.PFs(neuron).Dir2),2,'omitnan')>0))
                    xticks({})
                    xlabel('Stability')
                    ylabel('move # (Direction 2)')
                    colormap(subplot(2,10,6),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
                    set(gca,'YDir','normal')
            elseif ~isempty(Neurons.PFs(neuron).Dir2)
                line([Neurons.PFs(neuron).Dir2(1) Neurons.PFs(neuron).Dir2(1)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                line([Neurons.PFs(neuron).Dir2(end) Neurons.PFs(neuron).Dir2(end)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                title(['{\color{red} Not a PC} for Direction 2' newline ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir2)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir2,3)) '%' newline ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir2,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir2,3)) '%'])
            end
            
        subplot(2,10,11:17); hold on;
            toPlot1=results.C_raw(neuron,:); toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
            toPlot2=results.C_raw(neuron,:); toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
            plot(results.C_raw(neuron,:)); plot(toPlot1); plot(toPlot2)
            xlabel('Time (s)')
            title(['Activity Raw (Blue)' newline 'Movement periods in Direction 1/2 in Red/Yellow'])
                
        subplot(2,10,18:20)
            roifn = single(full(results.A));
            pixh = size(results.Cn, 1);
            pixw = size(results.Cn, 2);
            roifn = reshape(roifn, pixh, pixw,size(roifn,2));
            roifnMask = roifn./max(max(roifn));
            displayImg = results.PNR;
            roifnMask_bin = roifnMask > 0.5;
            imagesc(displayImg)
            hold on
            contour(roifnMask_bin(:,:,neuron),[1,1],'LineColor','r', 'linewidth', 1);
            xticks([])
            yticks([])
            title('Position of the ROI in the FOV')
       
    sgtitle(['Neuron #' num2str(neuron)]) 