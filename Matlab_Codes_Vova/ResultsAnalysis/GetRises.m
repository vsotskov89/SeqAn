function Neurons=GetRises(results,Params,Experiment,Neurons)
%% Define Activity

% I did a lot of experiment with the smoothing so feel free to change the method and the window because
% they are just the latest values I tested and not calibrated parameters
switch Params.Method.ActivityType
    case 'Calcium'
        Act=smoothdata(results.C_raw(:,1:Experiment.SizeData),2,'sgolay',50);                                       % Activity is the smoothed raw trace
    case 'Derivative'
        Act=smoothdata(diff(smoothdata(results.C_raw(:,1:Experiment.SizeData),2,'sgolay',50),1,2),2,'sgolay',50);   % Activity is the smoothed derivative of smoothed calcium traces        
        Act(Act<0)=0; Act = [zeros(Experiment.nNeurons,1) Act];
    case 'Spikes'
        Act=smoothdata(movsum(results.S,10,2),2);                                                                   % Activity is the smoothed moving sum of BSD spikes over 100ms (10 --> [10 0] if "true start" is wanted)
    case 'Binary'
        Act=ones(Experiment.nNeurons,Experiment.SizeData);                                                          % Activity is 1
end

%% Find Rises and restrict activity to rises
switch Params.Method.SmallEvents

   
    %_________ Only Rises, With small evts _________ 
    case 'RisesSE'        
        for neuron = 1:Experiment.nNeurons
            Rises=zeros(1,Experiment.SizeData);                                                             % Initialize matrix
            % run BSD again, find pPhys.threshold and recompute results.C using only what is above this threshold.     
%             [Ninf,~,Palg,Pphys,Oalg] = BSD(results.C_raw(neuron, :), Params.Obsd, Params.Pbsd); % Performs deconvolution.
%             delta = exp(-Oalg.dt/Palg.tauRise - Oalg.dt/Palg.tauDecay);
%             gamma = exp(-Oalg.dt/Palg.tauRise)+exp(-Oalg.dt/Palg.tauDecay)-exp(-Oalg.dt/Palg.tauRise - Oalg.dt/Palg.tauDecay) ;
%             eta = (Palg.tauRise/Palg.tauDecay)^(Palg.tauDecay/(Palg.tauDecay - Palg.tauRise))*( Palg.tauDecay/Palg.tauRise - 1)./(exp(-Oalg.dt/Palg.tauDecay)-exp(-Oalg.dt/Palg.tauRise));
%             UnderTresh = Ninf < Palg.threshold ;
%             s = Ninf; 
%             s(UnderTresh) = 0; 
%             c=filter(1,[eta, -eta * (gamma+delta)  ,  eta * delta  ],s);
                     
            noise = std(results.C_raw(neuron, :)-results.C(neuron, :));
            evts=results.C(neuron,1:Experiment.SizeData)>noise;   %0.1                                             % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
            RProps=regionprops(evts','area','PixelList');                                                   % Get continuous regions of events (i.e : separate the events)
            [~,maxIDs]=arrayfun(@(S)max(results.C(neuron,S.PixelList(:,2))),RProps);                        % Find maximum within each events
            onlyRises=cell2mat(arrayfun(@(S,I)S.PixelList(1:I,2),RProps,maxIDs,'UniformOutput',0));         % Rises are defined as the activity from event's start to maximum : discard the decay
            Rises(onlyRises)=Act(neuron,onlyRises);                                                        % Activity is restricted to this rises only
            dRises=diff(Rises>0); starts=find(dRises==1)+1; ends=find(dRises==-1);                          % Find starts/ends
            if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
            if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end    
            Neurons.Rises(neuron).Matrix=Rises;                                                             % Save values in structure
            Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;                           % Save values in structure
       
        
%             figure; 
%                 time = Experiment.timeFluo-Experiment.timeFluo(1);
%                 ax1 = subplot(3,1,1);
%                 plot(time, results.C_raw(neuron, :)); hold on;   plot(time,results.C(neuron, :)); plot(time, fResNorm(:,2)); %plot(time, c); % plot(sAct(neuron, :)./max(sAct(neuron, :)).*max(results.C(neuron,:))); plot(Neurons.Rises(neuron).Matrix./max(sAct(neuron, :)).*max(results.C(neuron,:)));
%                 ax2 = subplot(3,1,2);
%                plot(time, evts);ylim([-0.2 1.2])
%                 %plot(time, Neurons.Rises(neuron).Matrix); ylim([-0.2 1.2]) % hold  on;  plot(time, results.S(neuron,:)/10); % plot(time, evts); %plot(time, s); % plot(evts);  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));
%                 ax3 = subplot(3,1,3);
%                % plot(dFiltAct);  hold  on;   plot(dC); plot(sdFiltAct); plot(RisesIdx);% plot(results.S(neuron,:));
%                 toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
%                 toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
%                 plot(time, Experiment.positionXcmSmooth); hold on; plot(time, toPlot1); plot(time,toPlot2)    
%                 linkaxes([ax1,ax2, ax3],'x');
  
                
        end
        
     %_________ Rises Using Derivative Calculation
     
    case 'RisesDerivative'
        for neuron = 1:Experiment.nNeurons
            
         %   Rises=zeros(1,Experiment.SizeData);                                     % Initialize matrix
          %  noise = std(results.C_raw(neuron, :)-results.C(neuron, :));
           % evts=results.C(neuron,1:Experiment.SizeData)>noise;   %0.1              % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
 
            dataC = results.C_raw(neuron,:);
           % [dataC, ~] = remove_baseline_sliding_window(dataC, 10, 1000);    % we try a shorter sliding window : does not change much        
            
            % Test of a different method using filtered version of C_raw            
%             tdataC=[time dataC']; % Data to FMA samples
%             ftdataC=FilterLFP(tdataC,'passband',[0 40]); % this is not
%             what theo told me to do but I don't understand the
%             normalisation stuff

%             % calcul du bruit sur fResNorm(:,2);
%             [b, noise] = estimate_baseline_noise2(ftdataC(:,2));                % find baseline and noise
%             ftdataC(:,2) = ftdataC(:,2) - b;                                   % subtract residual baseline
%             evts1=ftdataC(:,2)>2*noise;                                         % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
%             gdtdataC = gradient(ftdataC(:,2));
%             [~, noiseD] = estimate_baseline_noise2(gdtdataC);                % find baseline and noise
%             Rises1 = Act(neuron,:);
%             Rises1(gdtdataC<noiseD)=0; 
%             Rises1(not(evts1)) = 0;
             
             
            % third method : median filter
            fdataC2 =  medfilt1(dataC,10);
            [~, noise2] = estimate_baseline_noise2(fdataC2);  
            evts2=fdataC2>2*noise2;                                            % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
            gfdataC2 = gradient(fdataC2);  %gfdataC3 = [0 diff(fdataC2)];
            noiseD2 = std(gfdataC2);
            Rises2 = Act(neuron,:);
            Rises2(gfdataC2<noiseD2)=0; 
            Rises2(not(evts2)) = 0;
            Neurons.Rises(neuron).Matrix=Rises2;                                                             % Save values in structure

            dRises2=diff(Rises2>0); starts=find(dRises2==1)+1; ends=find(dRises2==-1); 
            if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
            if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end 
            Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;

            
%             RProps=regionprops(evts','area','PixelList');                                                   % Get continuous regions of events (i.e : separate the events)
%             [~,maxIDs]=arrayfun(@(S)max(results.C(neuron,S.PixelList(:,2))),RProps);                        % Find maximum within each events
%             onlyRises=cell2mat(arrayfun(@(S,I)S.PixelList(1:I,2),RProps,maxIDs,'UniformOutput',0));         % Rises are defined as the activity from event's start to maximum : discard the decay
%             Rises(onlyRises)=Act(neuron,onlyRises);                                                        % Activity is restricted to this rises only
%             dRises=diff(Rises>0); starts=find(dRises==1)+1; ends=find(dRises==-1);                          % Find starts/ends
%             if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
%             if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end    
%             Neurons.Rises(neuron).Matrix=Rises2;                                                             % Save values in structure
%             Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;                           % Save values in structure


%             figure; 
%                 %time = Experiment.timeFluo-Experiment.timeFluo(1);
%                 time = 0: 0.01 : (size(results.C, 2)-1)*0.01;
%                 ax1 = subplot(4,1,1);
%                 plot(time, results.C_raw(neuron, :)); hold on;  plot(time,  fdataC2); % plot(time, ftdataC(:,2)); 
%                 x = [time(1) time(end)]; y = [0 0]; line(x,y,'Color','black','LineStyle','--')
%                 ax2 = subplot(4,1,2);
%                 plot(time, evts2);ylim([-0.2 1.2]); %hold on;  plot(time,1+ evts1); plot(time, 2+evts2);
%                 %plot(time, Neurons.Rises(neuron).Matrix);    hold  on;  plot(time, results.S(neuron,:)); %plot(time, SfmaFilt(:,2)); plot(time, Sthresh) % ylim([-0.2 1.2]) % plot(time, evts); %plot(time, s); % plot(evts);  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));
%                 ax3 = subplot(4,1,3);
%                 plot(time, gfdataC2); % hold on; plot(time, gfdataC3); 
%                 ax4 = subplot(4,1,4);
%                 plot(time, Rises2); %hold on; plot(time, Rises2);
%                 
% %                                 toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
% %                 toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
% %                 plot(time, Experiment.positionXcmSmooth); hold on; plot(time, toPlot1); plot(time,toPlot2)    
%                linkaxes([ax1,ax2, ax3, ax4],'x');
% %            
        end
        
        
     %_________ Rises and Decays, With small evts _________
     case 'RisesAndDecaysSE'
      
         for neuron = 1:Experiment.nNeurons
            Rises=zeros(1,Experiment.SizeData);                                                             % Initialize matrix
            noise = std(results.C_raw(neuron, :)-results.C(neuron, :));
            evts=results.C(neuron,1:Experiment.SizeData)>noise/3;   %0.1                                            % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
            %RProps=regionprops(evts','area','PixelList');                                                   % Get continuous regions of events (i.e : separate the events)
            %[~,maxIDs]=arrayfun(@(S)max(results.C(neuron,S.PixelList(:,2))),RProps);                        % Find maximum within each events
            %onlyRises=cell2mat(arrayfun(@(S,I)S.PixelList(1:I,2),RProps,maxIDs,'UniformOutput',0));         % Rises are defined as the activity from event's start to maximum : discard the decay
            Rises(evts)=Act(neuron,evts);                                                        % Activity is restricted to this rises only
            dRises=diff(Rises>0); starts=find(dRises==1)+1; ends=find(dRises==-1);                          % Find starts/ends
            if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
            if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end    
            Neurons.Rises(neuron).Matrix=Rises;                                                             % Save values in structure
            Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;                           % Save values in structure
     
            
%           figure; 
%                 time = Experiment.timeFluo-Experiment.timeFluo(1);
%                 ax1 = subplot(3,1,1);
%                 plot(time, results.C_raw(neuron, :)); hold on;   plot(time,results.C(neuron, :)); % plot(sAct(neuron, :)./max(sAct(neuron, :)).*max(results.C(neuron,:))); plot(Neurons.Rises(neuron).Matrix./max(sAct(neuron, :)).*max(results.C(neuron,:)));
%                 ax2 = subplot(3,1,2);
%                 plot(time, Neurons.Rises(neuron).Matrix);  hold  on;  plot(time, results.S(neuron,:)/10); plot(time, evts); % plot(evts);  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));
%                 ax3 = subplot(3,1,3);
%                % plot(dFiltAct);  hold  on;   plot(dC); plot(sdFiltAct); plot(RisesIdx);% plot(results.S(neuron,:));
%                 toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
%                 toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
%                 plot(time, Experiment.positionXcmSmooth); hold on; plot(time, toPlot1); plot(time,toPlot2)    
%                 linkaxes([ax1,ax2, ax3],'x');
        
        
        end 
        
        
%_________Only rises, With an Improved method _________ (work in progress)
    case 'RisesImprovedSE'
        
        sampleAct=[(1:Experiment.SizeData)' results.C_raw'];                                                        % Data to FMA Sample
        FiltAct = FilterLFP(sampleAct,'passband',12.5*Params.Method.filtRange); FiltAct=FiltAct(:,2:end)';          % Filter then Sample to Data
        
        for neuron=1:Experiment.nNeurons                                                                            % Loop on neurons
            dFiltAct=[0 diff(FiltAct(neuron,:))];                                                                   % Takes the derivative of the Activation
%             sdFiltAct=smooth(dFiltAct,2,'sgolay',10); 
            sdFiltAct=dFiltAct;
            sdposFiltAct=sdFiltAct; sdposFiltAct(sdposFiltAct<0)=0;                                                 % Smoothes the derivative and rectifies it
            dC=[0 diff(results.C(neuron,:))]; dposC=dC>0;                                                           % Creates a rectified C matrix
            
            RisesIdx=zeros(1,Experiment.SizeData);                                                      % Init Rises indexes matrix
            RisesIdx(results.C(neuron,:)>0.5 & sdFiltAct>0 & dposC>0)=1;                                % Events = Fluo is rising simultaneously with a CNMFE-detected events 
            init=find([0 diff(RisesIdx)==1]);                                                           % Find the starts of these events
            
            for evt = 2:numel(init)-1                                                                   % Loop on events (exclude first and last : NEED TO ADDRESS THESE CASES)
                StartEvt=find(sdposFiltAct(init(evt-1):init(evt))==0,1,'last')+init(evt-1)-1;           % True start of the event is defined as when the fluo traces starts rising
                if isempty(StartEvt); StartEvt=init(evt-1); end                                         % if it's rising since last event : ~merge the events

                EndEvt=find(sdposFiltAct(init(evt):init(evt+1))==0,1,'first')+init(evt)-1;              % True end is defined as when the fluo traces stops rising       
                if isempty(EndEvt); EndEvt=init(evt+1); end                                             % if it does not stop rising until last event : ~merge the events (probably redundant)
                
                RisesIdx(StartEvt:EndEvt)=1;                                                            % Correct the indexes matrix with new starts/ends
            end
            
            dRisesIdxN=diff(RisesIdx);                                                                  % Take the derivative of the indexes matrix
            starts=find(dRisesIdxN==1)+1; ends=find(dRisesIdxN==-1);                                    % Find starts/ends and durations
            if ~isempty(starts)
                if starts(1)>ends(1); starts=[1 starts]; end                                                % Boundaries checks
                if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end
                durations=ends-starts; 
                shortEvts=find(durations<Params.Method.minDur);                                             % Find short events
                for evt =1:numel(shortEvts)                                                                 % Loop on the short events
                    RisesIdx(starts(shortEvts(evt)):ends(shortEvts(evt)))=0;                                % Delete them
                end
                starts(shortEvts)=[]; ends(shortEvts)=[];
                
                Neurons.Rises(neuron).Matrix=Act(neuron,:);Neurons.Rises(neuron).Matrix(~RisesIdx)=0;                    % TO BE REMOVED
%                 Neurons.Rises(neuron).Matrix=sdposFiltAct; Neurons.Rises(neuron).Matrix(~RisesIdx)=0;                     % Take only the activity corresponding to the rises indexes
                
                Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;
            else
                Neurons.Rises(neuron).Matrix=zeros(1,size(sdposFiltAct,2));
            end
        end
        
%_________ Without small evts _________ NOT WORKING because of modifications made in new programs
%     case 'noSE'
%         Rises=zeros(Experiment.nNeurons,Experiment.SizeData);                                          % Initialize
%         for neuron = 1:Experiment.nNeurons                                                             % Loop on neurons
%             act=results.C_raw(neuron,1:Experiment.SizeData);
%             sAct = smooth(results.C(neuron,1:Experiment.SizeData),20);                                 % Smooth
%             DsAct= diff(sAct);                                                              % Take deriv
%             Peak=DsAct>0;                                                                   % Find periods of rises
% 
%             negPartIx=act<0;                                                                % Take the negative part of the trace
%             thresholdAct=abs(prctile(act(negPartIx),5));                                    % Threshold of activity is defined as the 95th percentile of the noise (negative part is noise only)
% 
%             RProps=regionprops(Peak,'area','PixelList');                                    % Get contiguous region of rising fluorescence
%             smallRegions=arrayfun(@(S)S.Area <30,RProps);                                   % Delete short periods (<30 frames)
%             RProps(smallRegions)=[];
% 
%             HighPx=find(act>thresholdAct);                                                  % Find regions that are above activity threshold
%             lowRegions=arrayfun(@(S)numel(intersect(S.PixelList(:,2),HighPx))<15,RProps);   % Find rise periods which do not cross this threshold (for >15 frames) and delete them
%             RProps(lowRegions)=[];
% 
%             RProps=struct2cell(RProps); RProps=cell2mat(RProps(2,:)');                      % Format Data
%             if ~isempty(RProps)
%                 RProps=RProps(:,2)';                        
%                 Rises(neuron,RProps)=act(RProps);                                           % Get a trace with only Rise events 
%             end
%         end
            
end