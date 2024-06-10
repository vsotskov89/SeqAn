function coAct=GetCoactivationMat(putSCEtimes,sceNeuronIDs,nNeurons)

coAct=zeros(nNeurons);
for sce = 1:numel(putSCEtimes)
    neurons=sceNeuronIDs{sce};
    for n=1:numel(neurons)
        coAct(neurons(n),neurons(n+1:numel(neurons)))=coAct(neurons(n),neurons(n+1:numel(neurons)))+1;
    end
end

figure
    imagesc(coAct)
