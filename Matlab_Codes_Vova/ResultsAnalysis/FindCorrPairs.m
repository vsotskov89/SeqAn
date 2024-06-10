function Neurons=FindCorrPairs(Neurons,Experiment,results)

pixh = size(results.Cn, 1); pixw = size(results.Cn, 2);
roifn = single(full(results.A)); roifn = reshape(roifn, pixh, pixw,size(roifn,2));
roifnMask = roifn./max(max(roifn)); roifnMask_bin = roifnMask > 0.5;

% Calculating Neurons Coordinates
for neuron=1:Experiment.nNeurons
    [x,y]=find(roifnMask_bin(:,:,neuron));
    Neurons.Coordinates(neuron).x=mean(x);
    Neurons.Coordinates(neuron).y=mean(y);
end

% Calculating Distances between Neurons
Neurons.Distances=zeros(Experiment.nNeurons);
for neuron1=1:size(results.C_raw,1)
    for neuron2=1:size(results.C_raw,1)
        Neurons.Distances(neuron1,neuron2)= ...
            sqrt((Neurons.Coordinates(neuron1).x-Neurons.Coordinates(neuron2).x)^2 + (Neurons.Coordinates(neuron1).y-Neurons.Coordinates(neuron2).y)^2);
    end
end

% Compute a Correlation Matrix
AllStillRises=cell2mat({Neurons.Rises.Matrix}'); AllStillRises=AllStillRises(:,logical(Experiment.Still.Segments));
Neurons.CorrMatrix=corrcoef(AllStillRises');
% CorrMatrix=triu(CorrMatrix,1);


%%%%% SUPPRESSION NOT IMPLEMENTED