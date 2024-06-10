warning('off', 'all'); 

for neuron = 1 : Experiment.nNeurons
    for kDirection = 1 : 2
        
        Direction = sprintf('Dir%d',kDirection);

%         if kDirection == 1 
%             Direction = 'Dir1'; % not elegant...
%         else
%             Direction = 'Dir2';
%         end
        
        MvtZone = sum(Experiment.(Direction).TimeSpent>0,'omitnan')>0.4*Experiment.(Direction).NbPassage; %MvtZone represente les portions du track que la souris a visité pendant au moins 40% des passages

                % first copy what we had in the function GetActivity

        PlaceActivity(neuron).(Direction) = zeros(Experiment.(Direction).NbPassage,ceil(max(Experiment.positionXcmSmooth)));          % Matrix of Place-related Activities ; dims : # neuron, # mouvement, Activity

            for nSeg = 1:Experiment.(Direction).NbPassage                                                                                      % Loop on mouvement periods
                for time = Experiment.(Direction).Starts(nSeg):Experiment.(Direction).Ends(nSeg)                    % Loop on times within mouvement periods
                    PlaceActivity(neuron).(Direction)(nSeg,ceil(Experiment.positionXcmSmooth(time))) = ...
                        PlaceActivity(neuron).(Direction)(nSeg,ceil(Experiment.positionXcmSmooth(time)))+Neurons.Rises(neuron).Matrix(time);       % Add activity to the corresponding neuron, mouvement period, spatial bin
                end
            end    
        PlaceActivity(neuron).(Direction)=PlaceActivity(neuron).(Direction)(:,MvtZone)./Experiment.(Direction).TimeSpent(:,MvtZone);    % normalize activity by time spend in each in each time bin

        if kDirection == 2
            PlaceActivity(neuron).(Direction) = fliplr(PlaceActivity(neuron).(Direction));
        end
        %Then we want to compute something similar with the 'Still' period that follows this mouvement. 

        % first find the time stamps of the beginning and end of the still periods.
        StillStarts = Experiment.(Direction).Ends;
        StillEnds = zeros(size(StillStarts));

        for nSeg = 1 : numel(StillStarts)
             ii = find(Experiment.('Still').Starts == StillStarts(nSeg));
             if numel(ii) >0
                 StillEnds(nSeg) = Experiment.('Still').Ends(ii);
             else
                  StillEnds(nSeg) = StillStarts(nSeg);
             end
        end
        %[StillStarts, StillEnds]
        StillDuration = StillEnds-StillStarts+1;
        tavg = 20; 
        N = ceil(max(StillDuration)/ tavg);

        StillActivity(neuron).(Direction) = nan(Experiment.(Direction).NbPassage,N);

        for nSeg = 1 : numel(StillStarts)
                A = Neurons.Rises(neuron).Matrix(StillStarts(nSeg):StillEnds(nSeg)); %// data
                B=reshape([A(:); nan(mod(-numel(A),tavg),1)],tavg,[]);
                out = nanmean(B);
                StillActivity(neuron).(Direction)(nSeg, 1: ceil((StillDuration(nSeg))/ tavg)) = out;
        end

        % now we want to concatenate both matrices and normalize by the max of each line. 

        PlaceAndStillActivity(neuron).(Direction) = [PlaceActivity(neuron).(Direction), StillActivity(neuron).(Direction)];
     %   maxes=max(PlaceAndStillActivity(neuron).(Direction),[],2); maxes(maxes==0)=-1;                                  % Take maxes of each passage
     %   PlaceAndStillActivity(neuron).(Direction)=PlaceAndStillActivity(neuron).(Direction)./repmat(maxes,1,size(PlaceAndStillActivity(neuron).(Direction),2));   % Normalize activity by the maxes


    end
    
end

% then plot 

fs = 9;

for neuron =1:Experiment.nNeurons
    fprintf([repmat('\b',1,4) '%s\n'],num2str(sprintf('%03d',neuron)));
    
    figure('Position', [200 120 900 600],'Visible','on');
        subplot(3,10,2:10);
            I1=imagesc(PlaceAndStillActivity(neuron).Dir1); hold on;
                set(I1,'AlphaData',~isnan(PlaceAndStillActivity(neuron).Dir1))
                set(gca,'YDir','normal')
                imax = max(prctile(PlaceAndStillActivity(neuron).Dir1(:),99), prctile(PlaceAndStillActivity(neuron).Dir2(:),99)); 
                imax = max(imax, 1);
                clim([0 imax])
                %colorbar;
            P1=plot(Neurons.SummedPlaceActivity(neuron).Dir1,'LineWidth',2,'Color','r');
            yticklabels({})
           % xlabel('Position (cm)','FontSize', fs)
            if Neurons.isPC(neuron).Dir1
                line([Neurons.PFs(neuron).Dir1(1) Neurons.PFs(neuron).Dir1(1)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                line([Neurons.PFs(neuron).Dir1(end) Neurons.PFs(neuron).Dir1(end)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                title(['{\color{green} PC} for Direction 1'  ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir1)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir1,3)) '%'  ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir1,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir1,3)) '%'],'FontSize', fs)
                subplot(3,10,1)
                    I2=imagesc(sum(PlaceAndStillActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0);
                        set(I2,'AlphaData',~isnan(sum(PlaceAndStillActivity(neuron).Dir1(:,Neurons.PFs(neuron).Dir1),2,'omitnan')>0))
                    xticks({})
                    xlabel('Stability','FontSize', fs)
                    ylabel('move # (Direction 1)','FontSize', fs)
                    colormap(subplot(3,10,1),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
                    set(gca,'YDir','normal')
            elseif ~isempty(Neurons.PFs(neuron).Dir1)
                line([Neurons.PFs(neuron).Dir1(1) Neurons.PFs(neuron).Dir1(1)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                line([Neurons.PFs(neuron).Dir1(end) Neurons.PFs(neuron).Dir1(end)], [0 Experiment.Dir1.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                title(['{\color{red} Not a PC} for Direction 1'  ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir1)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir1,3)) '%'  ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir1,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir1,3)) '%'],'FontSize', fs)
            end
            
        subplot(3,10,12:20);
            I3=imagesc(PlaceAndStillActivity(neuron).Dir2); hold on;
                set(I3,'AlphaData',~isnan(PlaceAndStillActivity(neuron).Dir2))
                set(gca,'YDir','normal')
                clim([0 imax])
                %colorbar;
            P2=plot(fliplr(Neurons.SummedPlaceActivity(neuron).Dir2),'LineWidth',2,'Color','r');
            yticklabels({})
           % xlabel('Position (cm)','FontSize', fs)
           Noc = size(PlaceActivity(neuron).(Direction),2 ); %size of occupied zone (used to flip data)
            if Neurons.isPC(neuron).Dir2
                line([Noc-Neurons.PFs(neuron).Dir2(1)+1 Noc-Neurons.PFs(neuron).Dir2(1)+1], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                line([Noc-Neurons.PFs(neuron).Dir2(end)+1 Noc-Neurons.PFs(neuron).Dir2(end)+1], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','g','LineStyle',':')
                title(['{\color{green} PC} for Direction 2'  ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir2)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir2,3)) '%'  ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir2,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir2,3)) '%'],'FontSize', fs)
                subplot(3,10,11)
                    I4=imagesc(sum(PlaceAndStillActivity(neuron).Dir2(:,Neurons.PFs(neuron).Dir2),2,'omitnan')>0);
                        set(I4,'AlphaData',~isnan(sum(PlaceAndStillActivity(neuron).Dir2(:,Neurons.PFs(neuron).Dir2),2,'omitnan')>0))
                    xticks({})
                    xlabel('Stability','FontSize', fs)
                    ylabel('move # (Direction 2)','FontSize', fs)
                    colormap(subplot(3,10,11),[[1 0.3 0.3];[1 0.3 0.3];[0.3 1 0.3]])
                    set(gca,'YDir','normal')
            elseif ~isempty(Neurons.PFs(neuron).Dir2)
                line([Noc-Neurons.PFs(neuron).Dir2(1)+1 Noc-Neurons.PFs(neuron).Dir2(1)+1], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                line([Noc-Neurons.PFs(neuron).Dir2(end)+1 Noc-Neurons.PFs(neuron).Dir2(end)+1], [0 Experiment.Dir2.NbPassage+1],'LineWidth',2,'Color','r','LineStyle',':')
                title(['{\color{red} Not a PC} for Direction 2'  ...
                'PF width : ' num2str(numel(Neurons.PFs(neuron).Dir2)) 'cm ; Activity in PF : ' num2str(100*round(Neurons.Conc(neuron).Dir2,3)) '%'  ...
                'Peak at ' num2str(100*round(Neurons.Peak(neuron).Dir2,3)) '% of theor. max ; Stability : ' num2str(100*round(Neurons.Stab(neuron).Dir2,3)) '%'],'FontSize', fs)
            end
            
        subplot(3,10,21:27); hold on;
            data1 = results.C_raw(neuron,:);
            %data1 = Neurons.Rises(neuron).Matrix;
            toPlot1=data1; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
            toPlot2=data1; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
            plot(data1); plot(toPlot1); plot(toPlot2); %plot(results.AvgRoiNoBaseline(neuron,:)/20+100)
            xlabel('Time (s)','FontSize', fs)
            title(['Activity Raw (Blue)'  'Movement periods in Direction 1/2 in Red/Yellow'],'FontSize', fs)
                
        subplot(3,10,28:30)
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
            title('Position of the ROI in the FOV','FontSize', fs)
       
    sgtitle(['Neuron #' num2str(neuron)])    
    print(gcf,'-dpdf','-bestfit', [Params.PDFsFolder 'tmp\' num2str(sprintf('%03d',neuron))]);         % Print pdf
    close
                 
end


names=dir([Params.PDFsFolder 'tmp/*.pdf']); names =fullfile({names(:).folder}, {names(:).name});                           
append_pdfs([Params.PDFsFolder 'AllPDFs'], names{:})
delete([Params.PDFsFolder 'tmp/*.pdf'])
movefile([Params.PDFsFolder 'AllPDFs'],[Params.PDFsFolder '_MoveAndStill_' Params.Method.SmallEvents '.pdf'])   

