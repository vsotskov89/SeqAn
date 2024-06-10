function Neurons=FindTheta(Neurons,Experiment,results,Params)

Neurons.Theta = struct('Matrix',repmat({zeros(1,Experiment.SizeData)},Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),'Ends',cell(Experiment.nNeurons,1));

sampleAct=[(0:0.01:(Experiment.SizeData-1)*0.01)' results.C_raw'];                                                      % Data to FMA Sample
FiltActTheta = FilterLFP(sampleAct,'passband',12.5*Params.Theta.boundsTheta);                                           % Filter (Theta Band) then Sample to Data
if Params.Theta.Normalization;  FiltActFreqNorm = FilterLFP(sampleAct,'passband',12.5*Params.Theta.boundsFreqNorm);     % Filter (Normalization band) then Sample to Data
else;                           FiltActFreqNorm=zeros(Experiment.SizeData,Experiment.nNeurons+1); end

for neuron=1:Experiment.nNeurons                                                                                                        % Loop on neurons
    ThetaEvts = GetThetaEvts(FiltActTheta(:,[1 neuron+1]),FiltActFreqNorm(:,[1 neuron+1]),Params.Theta.Normalization,'durations',Params.Theta.Durations,'thresholds',Params.Theta.Thresholds);          % Use FMA FindRipples to find theta periods
    Neurons.Theta(neuron).Starts=zeros(1,size(ThetaEvts,1)); Neurons.Theta(neuron).Starts=zeros(1,size(ThetaEvts,1));                        % init starts/ends matrices
    for period =1:size(ThetaEvts,1)                                                                                                          % Loop on theta periods
        Neurons.Theta(neuron).Matrix(floor(ThetaEvts(period,1)*100+1):floor(ThetaEvts(period,3)*100))=1;                                                        % Keep periods timings
        Neurons.Theta(neuron).Starts(period)=ThetaEvts(period,1)*100+1; Neurons.Theta(neuron).Ends(period)=ThetaEvts(period,3)*100;               % Keep starts and ends
    end
    Neurons.Rises(neuron).HasTheta=zeros(1,numel(Neurons.Rises(neuron).Starts));                                                            % init "Theta In Event" matrix
    for evt=1:numel(Neurons.Rises(neuron).Starts)                                                                                           % Loop on events
        Overlap=sum(Neurons.Theta(neuron).Matrix(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))/Neurons.Rises(neuron).Duration(evt); % Calculate fraction on event with theta
        if Overlap > Params.Theta.minOverlap                                                                                                    % if overlap over a threshold
            Neurons.Rises(neuron).HasTheta(evt)=1;                                                                                                  % event is defined as having theta
        end
    end
end

