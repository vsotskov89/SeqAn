function AnalyseSequences(MainPath,expPaths,DeepCAD)

% PCs have been analysed elsewhere. Here just load and use the corresponding matrices 

NeuronsBehavior = load([ cell2mat(expPaths.Behavior.Position) '\Neurons.mat' ]) ;

SPAdir1=cell2mat({NeuronsBehavior.Neurons.SummedPlaceActivity(cell2mat({NeuronsBehavior.Neurons.isPC.Dir1})).Dir1}');
for neuron = 1:size(SPAdir1, 1)
    SPAdir1(neuron, :) = movmean(SPAdir1(neuron, :), 9);
end
[Maxes1,Maxes1Idx]  = max(SPAdir1,[],2); [~,sMAxes1Idx]=sort(Maxes1Idx); SPAdir1=SPAdir1./Maxes1;
figure
imagesc(SPAdir1(sMAxes1Idx,:));

SPAdir2=cell2mat({NeuronsBehavior.Neurons.SummedPlaceActivity(cell2mat({NeuronsBehavior.Neurons.isPC.Dir2})).Dir2}');
for neuron = 1:size(SPAdir2, 1)
    SPAdir2(neuron, :) = movmean(SPAdir2(neuron, :), 9);
end
[Maxes2,Maxes2Idx]  = max(SPAdir2,[],2); [~,sMAxes2Idx]=sort(Maxes2Idx); SPAdir2=SPAdir2./Maxes2;
figure
imagesc(SPAdir2(sMAxes2Idx,:));


% Now we create the matrix to plot. 
% We choose Dir2 : more promising. 

% syncSCEs : vecteur indiquant si chaque SCE est synchrone ou pas. 

nSyncSCEs = sum(syncSCEs);
SyncSCEsIdx = find(syncSCEs);
% plottons deja la premiere

nwinplot = 20;

EvtOrdered = true(size(SPAdir2, 1), nSyncSCEs*2*nwinplot);
ActMat = vertcat(Neurons.Rises.Matrix);
ActOrdered = max(max(ActMat))*ones(size(SPAdir2, 1), nSyncSCEs*2*nwinplot);

for kSCE = 1 : nSyncSCEs

    timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
    [~,iSCE] = min(abs(TTLstartsTimes-timeSCE));
    EvtOrdered(:,(2*kSCE-1)*nwinplot - nwinplot+1 : (2*kSCE-1)*nwinplot + nwinplot-1) = EvtStarts(sMAxes2Idx, iSCE - nwinplot+1  : iSCE + nwinplot-1); 
    ActOrdered(:,(2*kSCE-1)*nwinplot - nwinplot+1 : (2*kSCE-1)*nwinplot + nwinplot-1) = ActMat(sMAxes2Idx, iSCE - nwinplot+1  : iSCE + nwinplot-1); 
end

figure; imagesc(EvtOrdered); 
figure; imagesc(ActOrdered);


% moyenne centrée sur les ripples/
nwinplot2 = 30;
SumActOrdered =  zeros(size(SPAdir2, 1), 2*nwinplot2);
notPC = and (not(cell2mat({NeuronsBehavior.Neurons.isPC.Dir1})), not(cell2mat({NeuronsBehavior.Neurons.isPC.Dir2})));
figure; 
for kSCE = 1 : nSyncSCEs
    kSCE
    timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
    [~,iSCE] = min(abs(TTLstartsTimes-timeSCE));
    SumActOrdered = SumActOrdered + ActMat(sMAxes2Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2);
    subplot(2,4,1);imagesc(EvtStarts(sMAxes2Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,5);imagesc(ActMat(sMAxes2Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,2);imagesc(EvtStarts(sMAxes1Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,6);imagesc(ActMat(sMAxes1Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,3);imagesc(EvtStarts(notPC, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,7);imagesc(ActMat(notPC, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,4);imagesc(EvtStarts(:, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,4,8);imagesc(ActMat(:, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))

    pause; 
end
figure; imagesc(SumActOrdered);


% moyenne centrée sur les ripples/
nwinplot2 = 30;
figure; 
for kSCE = 1 : nSyncSCEs
    kSCE
    timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
    [~,iSCE] = min(abs(TTLstartsTimes-timeSCE));
    subplot(2,1,1);imagesc(EvtStarts(sMAxes1Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,1,2);imagesc(ActMat(sMAxes1Idx, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    pause; 
end


%%
% on veut faire la figure suivante. 
% on considere toutes les cellules. On les ordonnes selon le PF dans l'une
% des 2 directions, et ceux indépendamment de si elles sont PC ou pas. On
% plotte ensuite les SCEs a chaque fois qu'un ripple est synchrone. 
% On attribue une valeur differente a l'evenement selon si la cellules est
% PC seulement pour cette direction, PC pour les 2 directions ou pas PC du
% tout. 
% appelons ces matrices EvtStartDir1 et EvtStartDir2

% j'ai la flemme de faire ca elegamment, je fais une boucle. 
EvtStartsDir1 = zeros(size(EvtStarts));
for neuron =1 : nNeurons
    if NeuronsBehavior.Neurons.isPC(neuron).Dir1 == 1
        if NeuronsBehavior.Neurons.isPC(neuron).Dir2 == 1
            EvtStartsDir1(neuron,:) = 3*double(EvtStarts(neuron,:));
        else
            EvtStartsDir1(neuron,:) = 2*double(EvtStarts(neuron,:));
        end
    else
    EvtStartsDir1(neuron,:) = double(EvtStarts(neuron,:));
    end
end

EvtStartsDir2 = zeros(size(EvtStarts));
for neuron =1 : nNeurons
    if NeuronsBehavior.Neurons.isPC(neuron).Dir2 == 1
        if NeuronsBehavior.Neurons.isPC(neuron).Dir1 == 1
            EvtStartsDir2(neuron,:) = 3*double(EvtStarts(neuron,:));
        else
            EvtStartsDir2(neuron,:) = 2*double(EvtStarts(neuron,:));
        end
    else
    EvtStartsDir2(neuron,:) = double(EvtStarts(neuron,:));
    end
end

% maintenant je calcule aussi l'ordre des neurones comme en haut mais pour
% tous les neurones

SPAdir1all=cell2mat({NeuronsBehavior.Neurons.SummedPlaceActivity.Dir1}');
for neuron = 1:size(SPAdir1all, 1)
    SPAdir1all(neuron, :) = movmean(SPAdir1all(neuron, :), 9);
end
[Maxes1,Maxes1Idxall]  = max(SPAdir1all,[],2); [~,sMAxes1Idxall]=sort(Maxes1Idxall); SPAdir1all=SPAdir1all./Maxes1;
figure
%imagesc(SPAdir1all(sMAxes1Idxall,:));

SPAdir2all=cell2mat({NeuronsBehavior.Neurons.SummedPlaceActivity.Dir2}');
for neuron = 1:size(SPAdir2all, 1)
    SPAdir2all(neuron, :) = movmean(SPAdir2all(neuron, :), 9);
end
[Maxes2,Maxes2Idxall]  = max(SPAdir2all,[],2); [~,sMAxes2Idxall]=sort(Maxes2Idxall); SPAdir2all=SPAdir2all./Maxes2;
figure
%imagesc(SPAdir2all(sMAxes2Idxall,:));

figure;
for kSCE = 1 : nSyncSCEs
    kSCE
    timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
    [~,iSCE] = min(abs(TTLstartsTimes-timeSCE));
    subplot(2,2,1);imagesc(EvtStartsDir1(sMAxes1Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,2,3);imagesc(ActMat(sMAxes1Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,2,2);imagesc(EvtStartsDir2(sMAxes2Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,2,4);imagesc(ActMat(sMAxes2Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    pause; 
end


% now we plot the stability in color.
EvtStartsStabDir1 = zeros(size(EvtStarts));
EvtStartsStabDir2 = zeros(size(EvtStarts));

for neuron =1 : nNeurons
    EvtStartsStabDir1(neuron,:) = NeuronsBehavior.Neurons.Stab(neuron).Dir1 *double(EvtStarts(neuron,:));
    EvtStartsStabDir2(neuron,:) = NeuronsBehavior.Neurons.Stab(neuron).Dir2 *double(EvtStarts(neuron,:));
end

figure;
for kSCE = 1 : nSyncSCEs
    kSCE
    timeSCE = putSCEtimes(SyncSCEsIdx(kSCE)); %timing dans la ref ephy
    [~,iSCE] = min(abs(TTLstartsTimes-timeSCE));
    subplot(2,2,1);imagesc(EvtStartsStabDir1(sMAxes1Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1), [0, 1])
    subplot(2,2,3);imagesc(ActMat(sMAxes1Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    subplot(2,2,2);imagesc(EvtStartsStabDir2(sMAxes2Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1), [0, 1])
    subplot(2,2,4);imagesc(ActMat(sMAxes2Idxall, iSCE - nwinplot2+1  : iSCE + nwinplot2-1))
    pause; 
end



