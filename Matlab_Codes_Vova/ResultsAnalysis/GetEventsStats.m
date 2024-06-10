function [Neurons,Stats]=GetEventsStats(Neurons,Experiment,results,Params)

shift1=find(Experiment.Dir1.MvtZone,1,'first'); % Shift because Position are analyzed only in mvt zone
shift2=find(Experiment.Dir2.MvtZone,1,'first'); % so : 0-1cm is not the first bin of activity matrices

for neuron=1:Experiment.nNeurons 
    Neurons.Rises(neuron).EvtsCat=zeros(1,numel(Neurons.Rises(neuron).Starts));     % Init variables
    Neurons.Rises(neuron).inPF=zeros(1,numel(Neurons.Rises(neuron).Starts)); 
    noiseLvl= abs(prctile(results.C_raw(neuron,results.C_raw(neuron,:)<0),1));      % noise level is defined as the 1st percentile of the negative part (no event) of the fluo
    
    for evt=1:numel(Neurons.Rises(neuron).Starts)                                                                       % Loop on events
        sizeEvt=max(results.C_raw(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)))...             % Size of event is min - max (not really accurate for now)
            -min(results.C_raw(neuron,Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));                  %%%%%% FAIRE CA AVEC LE SIGNAL FILTR2 (plus pr�cise
        if  sizeEvt < Params.EvtsCategories(1)*noiseLvl                                                                     % Event categories depends of the ratio size/noiseLvl
            Neurons.Rises(neuron).EvtsCat(evt) = 1;
        elseif  sizeEvt > Params.EvtsCategories(2)*noiseLvl
            Neurons.Rises(neuron).EvtsCat(evt) = 3; 
        else
            Neurons.Rises(neuron).EvtsCat(evt) = 2;
        end
        [~,Peak] = max(Neurons.Rises(neuron).Matrix(Neurons.Rises(neuron).Starts(evt):Neurons.Rises(neuron).Ends(evt)));    % Find index of the peak 
        PosPeak=Experiment.positionXcmSmooth(Neurons.Rises(neuron).Starts(evt)+Peak-1);                                     % Find animal position at peak time
        if Neurons.isPC(neuron).Dir1 && Experiment.Dir1.Segments(Neurons.Rises(neuron).Starts(evt)+Peak-1) &&...              % Check if Peak happends in the PF
            and(PosPeak>Neurons.PFs(neuron).Dir1(1)+shift1,PosPeak<Neurons.PFs(neuron).Dir1(end)+shift1)
             Neurons.Rises(neuron).inPF(evt)=1;
        elseif Neurons.isPC(neuron).Dir2 && Experiment.Dir2.Segments(Neurons.Rises(neuron).Starts(evt)+Peak) &&...
            and(PosPeak>Neurons.PFs(neuron).Dir2(1)+shift2,PosPeak<Neurons.PFs(neuron).Dir2(end)+shift2)
             Neurons.Rises(neuron).inPF(evt)=1;
        end
    end
end

nbCats=3;

Stats.All.NbEvt=arrayfun(@(S)numel(S.EvtsCat),Neurons.Rises);
Stats.All.NbHasTheta=arrayfun(@(S)sum(S.HasTheta),Neurons.Rises);
Stats.All.FractHasTheta=Stats.All.NbHasTheta./Stats.All.NbEvt;

Stats.AllPF.NbEvt=arrayfun(@(S)numel(S.inPF),Neurons.Rises);
Stats.AllPF.FractEvt=Stats.AllPF.NbEvt./Stats.All.NbEvt;
% Stats.AllPF.NbHasTheta=arrayfun(@(S)sum(S.HasTheta & S.inPF),Neurons.Rises);
% Stats.AllPF.FractHasTheta=Stats.AllPF.NbHasTheta./Stats.AllPF.NbEvt;

for cat=1:nbCats
    Stats.(['Cat' num2str(cat)]).NbEvt=arrayfun(@(S)sum(S.EvtsCat==cat),Neurons.Rises);
    Stats.(['Cat' num2str(cat)]).FractEvt=Stats.(['Cat' num2str(cat)]).NbEvt./Stats.All.NbEvt;
%     Stats.(['Cat' num2str(cat)]).NbHasTheta=arrayfun(@(S)sum(S.HasTheta(S.EvtsCat==cat)),Neurons.Rises);
%     Stats.(['Cat' num2str(cat)]).FractHasTheta=Stats.(['Cat' num2str(cat)]).NbHasTheta./Stats.(['Cat' num2str(cat)]).NbEvt;
end

% figure
% subplot(2,2,1);hold on
%         histogram(Stats.Cat1.NbEvt,0:10:500,'FaceAlpha',0.4)
%         histogram(Stats.Cat2.NbEvt,0:10:500,'FaceAlpha',0.4)
%         histogram(Stats.Cat3.NbEvt,0:10:500,'FaceAlpha',0.4)
%         xlabel('Nb of evts in cat.'); ylabel('Nb of neurons') ; title('Nb of Evts within each categories')
%     subplot(2,2,2);hold on
%         histogram(Stats.Cat1.FractEvt*100,0:5:100,'FaceAlpha',0.4)
%         histogram(Stats.Cat2.FractEvt*100,0:5:100,'FaceAlpha',0.4)
%         histogram(Stats.Cat3.FractEvt*100,0:5:100,'FaceAlpha',0.4)
%         xlabel('% of evts in cat.'); ylabel('Nb of neurons') ; title('% of Evts within each categories')
%     subplot(2,2,3);hold on
%         histogram(Stats.Cat1.NbHasTheta,0:1:50,'FaceAlpha',0.4)
%         histogram(Stats.Cat2.NbHasTheta,0:1:50,'FaceAlpha',0.4)
%         histogram(Stats.Cat3.NbHasTheta,0:1:50,'FaceAlpha',0.4)
%         xlabel('Nb of evts with theta'); ylabel('Nb of neurons') ; title('Nb of Evts with theta for each cat. of evts')
%     subplot(2,2,4);hold on
%         histogram(Stats.Cat1.FractHasTheta*100,0:5:100,'FaceAlpha',0.4)
%         histogram(Stats.Cat2.FractHasTheta*100,0:5:100,'FaceAlpha',0.4)
%         histogram(Stats.Cat3.FractHasTheta*100,0:5:100,'FaceAlpha',0.4)
%         xlabel('% of evts with theta'); ylabel('Nb of neurons') ; title('% of Evts with theta for each cat. of evts')
% sgtitle(['Recap ' Experiment.file ' : Event Categories and Theta' newline 'Parameters : Theta freq. : ' num2str(Params.Theta.boundsTheta(1)) '-' num2str(Params.Theta.boundsTheta(2))...
%     'Hz ; Thresholds : Peak >' num2str(Params.Theta.Thresholds(2)) 'std, Lims >' num2str(Params.Theta.Thresholds(1)) 'std'])


