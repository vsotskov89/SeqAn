function Neurons=GetRises(results,Params,Experiment,Neurons)
% Find events in each neuron's activity

switch Params.Method.SmallEvents
%_________ Without small evts _________ NOT WORKING because of modifications made in new programs
    case 'noSE'
        Rises=zeros(Experiment.nNeurons,Experiment.SizeData);                                          % Initialize
        for neuron = 1:Experiment.nNeurons                                                             % Loop on neurons
            act=results.C_raw(neuron,1:Experiment.SizeData);
            sAct = smooth(results.C(neuron,1:Experiment.SizeData),20);                                 % Smooth
            DsAct= diff(sAct);                                                              % Take deriv
            Peak=DsAct>0;                                                                   % Find periods of rises

            negPartIx=act<0;                                                                % Take the negative part of the trace
            thresholdAct=abs(prctile(act(negPartIx),5));                                    % Threshold of activity is defined as the 95th percentile of the noise (negative part is noise only)

            RProps=regionprops(Peak,'area','PixelList');                                    % Get contiguous region of rising fluorescence
            smallRegions=arrayfun(@(S)S.Area <30,RProps);                                   % Delete short periods (<30 frames)
            RProps(smallRegions)=[];

            HighPx=find(act>thresholdAct);                                                  % Find regions that are above activity threshold
            lowRegions=arrayfun(@(S)numel(intersect(S.PixelList(:,2),HighPx))<15,RProps);   % Find rise periods which do not cross this threshold (for >15 frames) and delete them
            RProps(lowRegions)=[];

            RProps=struct2cell(RProps); RProps=cell2mat(RProps(2,:)');                      % Format Data
            if ~isempty(RProps)
                RProps=RProps(:,2)';                        
                Rises(neuron,RProps)=act(RProps);                                           % Get a trace with only Rise events 
            end
        end
    
%_________ With small evts _________ 
    case 'SE'        
         filtSize = 30;
    %     sAct=smoothdata(results.C_raw(:,1:SizeData),2,'sgolay',50);                                                   % Activity is smoothed calcium traces
    %    sAct=(diff(smoothdata(results.C_raw(:,1:Experiment.SizeData),2,'sgolay',50),1,2));     
         sAct=smoothdata(diff(smoothdata(results.C_raw(:,1:Experiment.SizeData),2,'lowess',filtSize),1,2),2,'lowess',filtSize);      % Activity is the smoothed derivative of smoothed calcium traces        
        sAct(sAct<0)=0; sAct = [zeros(Experiment.nNeurons,1) sAct];                                                     % negative to 0 because we do not want "negative activity"
        
        for neuron = 1:Experiment.nNeurons
            Rises=zeros(1,Experiment.SizeData);                                                             % Initialize matrix
            evts=results.C(neuron,1:Experiment.SizeData)>0.1;                                               % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
            RProps=regionprops(evts','area','PixelList');                                                   % Get continuous regions of events (i.e : separate the events)
            [~,maxIDs]=arrayfun(@(S)max(results.C(neuron,S.PixelList(:,2))),RProps);                        % Find maximum within each events
            onlyRises=cell2mat(arrayfun(@(S,I)S.PixelList(1:I,2),RProps,maxIDs,'UniformOutput',0));         % Rises are defined as the activity from event's start to maximum : discard the decay
            Rises(onlyRises)=sAct(neuron,onlyRises);                                                        % Activity is restricted to this rises only
            dRises=diff(Rises>0); starts=find(dRises==1)+1; ends=find(dRises==-1);                          % Find starts/ends
            if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
            if starts(end)>ends(end); ends=[ends Experiment.SizeData]; end    
            Neurons.Rises(neuron).Matrix=Rises;                                                             % Save values in structure
            Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;   % Save values in structure
          
            
            noise = std(results.C_raw(neuron, :)-results.C(neuron, :));
    
              figure; 
                ax1 = subplot(2,1,1);
                plot(results.C_raw(neuron, :)); hold on;   plot(results.C(neuron, :));  plot(sAct(neuron, :)./max(sAct(neuron, :)).*max(results.C(neuron,:))); plot(Neurons.Rises(neuron).Matrix./max(sAct(neuron, :)).*max(results.C(neuron,:)));
                ax2 = subplot(2,1,2);
                plot(Neurons.Rises(neuron).Matrix);  hold  on;  plot(results.S(neuron,:)/10); plot(evts);%  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));
               % ax3 = subplot(3,1,3);
               % plot(dFiltAct);  hold  on;   plot(dC); plot(sdFiltAct); plot(RisesIdx);% plot(results.S(neuron,:));
                linkaxes([ax1,ax2],'x');
            
                figure; 
                plot(results.C_raw(neuron, :)/noise); hold on; plot(Neurons.Rises(neuron).Matrix);
        end
        
%_________ With an Improved method _________ (work in progress)
    case 'improvedSE'
        filtSize = 50;
        sAct=smoothdata(diff(smoothdata(results.C_raw(:,1:Experiment.SizeData),2,'sgolay',filtSize),1,2),2,'sgolay',filtSize);      % Activity is the smoothed derivative of smoothed calcium traces        
        sAct(sAct<0)=0; sAct = [zeros(Experiment.nNeurons,1) sAct];
        
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
                
                Neurons.Rises(neuron).Matrix=sAct(neuron,:);Neurons.Rises(neuron).Matrix(~RisesIdx)=0;                    % TO BE REMOVED
%                 Neurons.Rises(neuron).Matrix=sdposFiltAct; Neurons.Rises(neuron).Matrix(~RisesIdx)=0;                     % Take only the activity corresponding to the rises indexes
                
                Neurons.Rises(neuron).Starts=starts; Neurons.Rises(neuron).Ends=ends; Neurons.Rises(neuron).Duration=ends-starts;
            else
                Neurons.Rises(neuron).Matrix=zeros(1,size(sdposFiltAct,2));
            end
            
%                        % Check data
%                 figure; 
%                 ax1 = subplot(3,1,1);
%                 plot(results.C_raw(neuron, :)); hold on;   plot(results.C(neuron, :)); plot(FiltAct(neuron,:)); plot(sAct(neuron, :)./max(sAct(neuron, :)).*max(FiltAct(neuron,:)));
%                 ax2 = subplot(3,1,2);
%                 plot(Neurons.Rises(neuron).Matrix);  hold  on;  plot(results.S(neuron,:)/10); %  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));
%                 ax3 = subplot(3,1,3);
%                 plot(dFiltAct);  hold  on;   plot(dC); plot(sdFiltAct); plot(RisesIdx);% plot(results.S(neuron,:));
%                 linkaxes([ax1,ax2,ax3],'x');

        end
        
end