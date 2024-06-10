function [Neurons,Experiment]=GetPCstatus(Experiment,Neurons,Params,Direction)

% old criteria
% for neuron=1:Experiment.nNeurons                                                                    % loop on neurons
%     Neurons.isPC(neuron).(Direction)= ...
%         Neurons.Peak(neuron).(Direction) > ...
%             Params.critPCs.thrPeak - Params.critPCs.StabDiscFactor*Neurons.Stab(neuron).(Direction);        % Cell is PC if it pass this criterion
% end

for neuron=1:Experiment.nNeurons                                                                    % loop on neurons
    Neurons.isPC(neuron).(Direction)= ...
        (Neurons.Stab(neuron).(Direction) > Params.critPCs.thrPeak) ...
    * ( Neurons.ActvityInPFvsOutPF(neuron).(Direction) > Params.critPCs.thrAct ) ; 
end

Experiment.PCs.(Direction) = find(cell2mat({Neurons.isPC.(Direction)}));
