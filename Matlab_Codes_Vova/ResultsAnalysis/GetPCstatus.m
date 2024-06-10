function [Neurons,Experiment]=GetPCstatus(Experiment,Neurons,Params,Direction)

% old criteria
% for neuron=1:Experiment.nNeurons                                                                    % loop on neurons
%     Neurons.isPC(neuron).(Direction)= ...
%         Neurons.Peak(neuron).(Direction) > ...
%             Params.critPCs.thrPeak - Params.critPCs.StabDiscFactor*Neurons.Stab(neuron).(Direction);        % Cell is PC if it pass this criterion
% end

% for neuron=1:Experiment.nNeurons                                                                    % loop on neurons
%     Neurons.isPC(neuron).(Direction)= ...
%         (Neurons.Stab(neuron).(Direction) > Params.critPCs.thrStab) ...
%     & ( Neurons.ActvityInPFvsOutPF(neuron).(Direction) > Params.critPCs.thrAct ) ; 
% end
% 
% Experiment.PCs.(Direction) = find(cell2mat({Neurons.isPC.(Direction)}));


Experiment.PCs.(Direction) = zeros(Experiment.nNeurons, 1); 

for neuron=1:Experiment.nNeurons                                                                    % loop on neurons
    nPF = size(Neurons.PFs(neuron).(Direction),1);
    Neurons.isPC(neuron).(Direction) = cell(nPF,1); 
    for kPF = 1 : nPF 
        isPC = cell2mat(Neurons.Stab(neuron).(Direction)(kPF))> Params.critPCs.thrStab;     % stability criteria
        Experiment.PCs.(Direction) = Experiment.PCs.(Direction) + isPC;
        Neurons.isPC(neuron).(Direction)(kPF)= { isPC } ;
    end
    Neurons.isPC2(neuron).(Direction) = (sum(cell2mat(Neurons.isPC(neuron).(Direction))) > 0);
end

Experiment.PCs.(Direction) = find(Experiment.PCs.(Direction));
