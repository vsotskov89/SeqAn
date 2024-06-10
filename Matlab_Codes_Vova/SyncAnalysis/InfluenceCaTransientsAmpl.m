function InfluenceCaTransientsAmpl(ExtractedResults,TTLstartsTimes,EvtStarts,ripples, nNeurons )


% First we observe the data. We plot for the first neurons, C_raw,
% EvtStarts et la position des ripples (bandes vertes). 
figure; 
NneuronFig = 60;
for neuron = 1 : NneuronFig
    plot(TTLstartsTimes, ExtractedResults.C_raw(neuron,:)+100*(neuron-1)); 
    hold on; 
    plot(TTLstartsTimes, EvtStarts(neuron,:)*100+100*(neuron-1), 'k');
end
for kripple = 1 : size(ripples, 1)
     xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
     yBox = [0, 100*NneuronFig, 100*NneuronFig, 0, 0];
     patch(xBox, yBox, 'green', 'FaceAlpha', 0.2, 'FaceColor', 'green', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Then we calculte the Ca Evt amplitude for each Evt. 
EvtAmpl = double(EvtStarts);
[raw, col, v] = find(EvtStarts);

dataFilt = zeros(size(ExtractedResults.C_raw));
for k = 1: size(dataFilt,1)
        dataFilt(k,:) = medfilt1(ExtractedResults.C_raw(k,:),10);
end

%first idea, we compare with the amplitude of the local transient. This
%actually does not work, but maybe it is normal because the amplitude of a
%transient is always the same ?
%countZeros = 0;
% for k = 1:numel(raw)
%     EvtAmpl(raw(k),col(k)) = 1/2*sum((dataFilt(raw(k), col(k)-1 : col(k)+2) .* [-1 -1 1 1])) + 1e-8;
%     if EvtAmpl(raw(k),col(k)) == 0
%         countZeros = countZeros +1;
%     end
% end
% second idea. We take the DeltaF at the evt start. 
% for k = 1:numel(raw)
%     EvtAmpl(raw(k),col(k)) = dataFilt(raw(k), col(k)) + 1e-8;
% end
% third idea, we take the maximum of the fltered signal in the following
% 300ms
for k = 1:numel(raw)
    EvtAmpl(raw(k),col(k)) = max(dataFilt(raw(k), col(k): min(col(k)+30,end)));
end


[raw, col, Ampl] = find(EvtAmpl);
%Ampl = EvtAmpl([raw, col]);
figure; histogram(Ampl)

% create a matrix indicating proximity of Evt to Ripple
% EvtRippleProx = 1 if a ripple is close to the event, -1 if not. 

minDist = 0.1;  
EvtRippleProx = -EvtStarts;
for k = 1:numel(raw)
    DistClosestRip = min(abs(ripples(:,2)-TTLstartsTimes(col(k))));
    if DistClosestRip < minDist
        EvtRippleProx(raw(k),col(k)) = 1;
    end
end

% verification graphique que le calcul est OK
figure; [yPoints,xPoints] = find(EvtRippleProx==1); plot(TTLstartsTimes(xPoints),yPoints,'.k');
hold on;
for kripple = 1 : size(ripples, 1)
    %kripple = 1;
    xBox = [ripples(kripple, 1), ripples(kripple, 1), ripples(kripple, 3), ripples(kripple, 3), ripples(kripple, 1)];
    yBox = [0, nNeurons, nNeurons, 0, 0];
    patch(xBox, yBox, 'green', 'FaceAlpha', 0.1, 'FaceColor', 'green', 'FaceAlpha', 0.1);
end

[raw, col, Prox] = find(EvtRippleProx);

%figure; scatter(Prox, Ampl);
AmplClose = Ampl(Prox==1);
AmplFar = Ampl(Prox == -1);
figure; 
h1 = histogram(AmplClose); 
hold on; 
h2 = histogram(AmplFar);
h1.Normalization = 'probability';  
h1.BinWidth = 5;
%h1.FaceColor = 'b';
h2.Normalization = 'probability'; 
h2.BinWidth = 5;
%h2.FaceColor = 'r';
