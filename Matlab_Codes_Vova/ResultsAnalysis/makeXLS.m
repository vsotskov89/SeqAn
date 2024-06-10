function Params=makeXLS(Neurons,Experiment,Params,Line)

if Line==1
    Params.xlsName=[Params.XLSsFolder datestr(now,'ddmm_HHMM_') 'SummaryTransients.xlsx'];         % Name
    xlsHead={'Manip','Tot. Cell', 'Tot. PCs','PCs Dir. 1','PCs Dir. 2',...                  % Header
        'Totat Nb Transients','PCs Nb Transients','Non-PCs Nb Transients',...
        'Total Transients cat. 1','PCs Transients cat. 1','Non-PCs Transients cat. 1',...
        'Total Transients cat. 2','PCs Transients cat. 2','Non-PCs Transients cat. 2',...
        'Total Transients cat. 3','PCs Transients cat. 3','Non-PCs Transients cat. 3'};
    xlswrite(Params.xlsName,xlsHead,1,'A1'); % Initialize XLS
end

PCs=cell2mat({Neurons.isPC.Dir2})|cell2mat({Neurons.isPC.Dir1});

starts={Neurons.Rises.Starts};
nbStarts=sum(cellfun(@(C)numel(C),starts)); nbStartsPCs=sum(cellfun(@(C)numel(C),starts(PCs))); nbStartsNonPCs=sum(cellfun(@(C)numel(C),starts(~PCs)));

Categories={Neurons.Rises.EvtsCat};
Cat1=sum(cellfun(@(C)sum(C==1),Categories)); Cat1PCs=sum(cellfun(@(C)sum(C==1),Categories(PCs))); Cat1NonPCs=sum(cellfun(@(C)sum(C==1),Categories(~PCs)));
Cat2=sum(cellfun(@(C)sum(C==2),Categories)); Cat2PCs=sum(cellfun(@(C)sum(C==2),Categories(PCs))); Cat2NonPCs=sum(cellfun(@(C)sum(C==2),Categories(~PCs)));
Cat3=sum(cellfun(@(C)sum(C==3),Categories)); Cat3PCs=sum(cellfun(@(C)sum(C==3),Categories(PCs))); Cat3NonPCs=sum(cellfun(@(C)sum(C==3),Categories(~PCs)));

xlsSumm ={Experiment.file, Experiment.nNeurons, sum(PCs), numel(Experiment.PCs.Dir1) , numel(Experiment.PCs.Dir2),...
    nbStarts,nbStartsPCs,nbStartsNonPCs,...
    Cat1,Cat1PCs,Cat1NonPCs,...
    Cat2,Cat2PCs,Cat2NonPCs,...
    Cat3,Cat3PCs,Cat3NonPCs};
xlswrite(Params.xlsName,xlsSumm,1,['A' num2str(Line+1)]);