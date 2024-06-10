function [nbStarts,NeuronIDs]=GetNbStarts(EvtStarts,window)

NeuronIDs=cell(1,size(EvtStarts,2)-window);
nbStarts=zeros(1,size(EvtStarts,2)-window);
%nbStarts2=zeros(1,size(EvtStarts,2)-window);


for t=1:size(EvtStarts,2)-window
    actNeurons=sum(EvtStarts(:,t:t+window-1),2)>0; %neurons that have an event starting during the considered window
    nbStarts(t)=sum(actNeurons); %number of neurons that have an event starting during the considered window
     NeuronIDs{t}=find(actNeurons); % with this expression, when a neuron  has 2 events starting during the considered window it is counted once
%    [NeuronIDs{t},~]=find(EvtStarts(:,t:t+window-1)); % with this expression, when a neuron  has 2 events starting during the considered window it is counted twice
%    nbStarts2(t) = numel(NeuronIDs{t});
end

%figure; plot(nbStarts(1,:)); hold on; plot(nbStarts2(1,:)); 


