% on plotte les summaries normaux mais on veut aussi la position de la
% souris et sa vitesse; 

for neuron = 1 : 121
 figure('Position', [100 100 1800 1200],'Visible','on'); % 'off' [200 120 900 600]
    
    normFact = max(max(Neurons.SummedPlaceActivity(neuron).Dir1),max(Neurons.SummedPlaceActivity(neuron).Dir2))/Experiment.Dir1.NbPassage*1.05;
    
        subplot(3,10,2:5);
            I1=imagesc(Neurons.PlaceActivity(neuron).Dir1); hold on;
                set(I1,'AlphaData',~isnan(Neurons.PlaceActivity(neuron).Dir1))
                set(gca,'YDir','normal')
                set(gca, 'color', 'k');
                clim([0 max(prctile(Neurons.PlaceActivity(neuron).Dir1(:),99),0.8)])
            %P1=plot(5*Neurons.SummedPlaceActivity(neuron).Dir1,'LineWidth',2,'Color','r');
            %hold on;  
            plot(Neurons.SummedPlaceActivity(neuron).Dir1/normFact,'LineWidth',2,'Color','m');
            %plot(25*movmean(Neurons.SummedPlaceActivity(neuron).Dir1,Params.smoothActivitySpace)./max(movmean(Neurons.SummedPlaceActivity(neuron).Dir1,Params.smoothActivitySpace)),'LineWidth',2,'Color','m'); 
            yticklabels({})
            xlabel('Position (cm)')

            nPF =  size(Neurons.PFs(neuron).Dir1,1);
            reliability = zeros(nPF, Experiment.Dir1.NbPassage);
            if sum(cell2mat(Neurons.isPC(neuron).Dir1)) > 0
                text1 = '{\color{green} PC} for Direction 1'; 
            else
                text1 = '{\color{red} Not a PC} for Direction 1';
            end
            title([  text1  newline  'Repetability : ' mat2str(100*round(cell2mat(Neurons.Stab(neuron).Dir1), 3)') ' %'])
     
            for kPF = 1 : nPF
                if cell2mat(Neurons.isPC(neuron).Dir1(kPF)) > 0
                    lineColor = 'g';
                else
                    lineColor = 'r';
                end
            xBounds = cell2mat(Neurons.PFs(neuron).Dir1(kPF));
            line([ xBounds(1) xBounds(1)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color',lineColor,'LineStyle',':')
            line([ xBounds(end) xBounds(end)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color',lineColor,'LineStyle',':')
            reliability(kPF,:)= sum(Neurons.PlaceActivity(neuron).Dir1(:,cell2mat(Neurons.PFs(neuron).Dir1(kPF))),2,'omitnan')>0;
            end
            
       subplot(3,10,1)
            I2=imagesc(reliability');
            %  set(I2,'AlphaData',~isnan(sum(Neurons.PlaceActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0))
            xticks({})
            xlabel('Reliability')
            ylabel('move # (Direction 1)')
            colormap(subplot(3,10,1),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
            set(gca,'YDir','normal')
           
        
       subplot(3,10,7:10);
            I3=imagesc(Neurons.PlaceActivity(neuron).Dir2); hold on;
                set(I3,'AlphaData',~isnan(Neurons.PlaceActivity(neuron).Dir2))
                set(gca,'YDir','normal')
                set(gca, 'color', 'k');
                clim([0 max(prctile(Neurons.PlaceActivity(neuron).Dir2(:),99),0.8)])
            %P2=plot(5*Neurons.SummedPlaceActivity(neuron).Dir2,'LineWidth',2,'Color','r');
            %hold on; 
            plot(Neurons.SummedPlaceActivity(neuron).Dir2/normFact,'LineWidth',2,'Color','m');
            %plot(25*movmean(Neurons.SummedPlaceActivity(neuron).Dir2,Params.smoothActivitySpace)./max(movmean(Neurons.SummedPlaceActivity(neuron).Dir2,Params.smoothActivitySpace)),'LineWidth',2,'Color','m'); 
            yticklabels({})
            xlabel('Position (cm)')

            nPF =  size(Neurons.PFs(neuron).Dir2,1);
            reliability = zeros(nPF, Experiment.Dir2.NbPassage);
            if sum(cell2mat(Neurons.isPC(neuron).Dir2)) > 0
                text1 = '{\color{green} PC} for Direction 2'; 
            else
                text1 = '{\color{red} Not a PC} for Direction 2';
            end
            title([  text1  newline  'Repetability : ' mat2str(100*round(cell2mat(Neurons.Stab(neuron).Dir2), 3)') ' %'])
     
            for kPF = 1 : nPF
                if cell2mat(Neurons.isPC(neuron).Dir2(kPF)) > 0
                    lineColor = 'g';
                else
                    lineColor = 'r';
                end
            xBounds = cell2mat(Neurons.PFs(neuron).Dir2(kPF));
            line([ xBounds(1) xBounds(1)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color',lineColor,'LineStyle',':')
            line([ xBounds(end) xBounds(end)], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color',lineColor,'LineStyle',':')
            reliability(kPF,:)= sum(Neurons.PlaceActivity(neuron).Dir2(:,cell2mat(Neurons.PFs(neuron).Dir2(kPF))),2,'omitnan')>0;
            end
            
       subplot(3,10,6)
            I4=imagesc(reliability');
            %  set(I2,'AlphaData',~isnan(sum(Neurons.PlaceActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0))
            xticks({})
            xlabel('Reliability')
            ylabel('move # (Direction 2)')
            colormap(subplot(3,10,6),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
            set(gca,'YDir','normal')       
            
        ax1 = subplot(3,10, 11:17); hold on;
            toPlot1=results.C_raw(neuron,:); toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
            toPlot2=results.C_raw(neuron,:); toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
            plot(results.C_raw(neuron,:)); plot(toPlot1); plot(toPlot2)
            xlabel('Time (s)')
            title('Activity during rest (blue) and movement in direction 1/2 (red/yellow)')
                
        ax2 = subplot(3,10, 21:27); hold on;
            toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
            toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
            plot(Experiment.positionXcmSmooth); hold on; plot(toPlot1); plot(toPlot2)  
            
        linkaxes([ax1,ax2],'x');
             
        subplot(3,10,18:20)
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
    set(gcf, 'color', 'w');
    set(gcf, 'InvertHardCopy','off')
    pause; 
end