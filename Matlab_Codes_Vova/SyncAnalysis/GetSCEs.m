function [putSCE,putSCEidx,putSCEtimes]=GetSCEs(nbStarts,nbStartsThr,EvtStarts,wndw,TTLstartsTimes)

putSCE=and(nbStarts>=nbStartsThr,islocalmax(nbStarts,'MinSeparation',10));
putSCEidx=find(putSCE);

for sce = 1:numel(putSCEidx)
    putSCEidx(sce)=putSCEidx(sce)+find(sum(EvtStarts(:,putSCEidx(sce):putSCEidx(sce)+wndw-1)),1)-1;
end
putSCE=Unfind(putSCEidx,numel(putSCE));
putSCEtimes=TTLstartsTimes(putSCEidx);


