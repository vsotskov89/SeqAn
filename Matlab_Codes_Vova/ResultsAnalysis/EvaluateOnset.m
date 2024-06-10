
MainPath='G:\Souris154148\22-06-30'; ePhysChan=8; TTLChan=18;nbExp=[1,1,1];

DeepCAD='';                                                                 % If you want to use DeepCAD or not. In the concatenated result folder, the matrices must
nDownsample = 1;                                                               % N_downsample>1 To analyse data at lower resolution. Requires a file results_XXHz.mat computed using CNMFE
% N_downsample = 1 to analyse data at 100Hz
% Find path to files
expPaths=GetExpsPaths('NbExp',nbExp,'Path',MainPath);                                                       % Load the paths to the differents file (or ask the user for them if it wasn't done before)

% Load or extract the calcium results file
results=LoadExtractedCalciumData(MainPath,expPaths,DeepCAD,nDownsample,'');                                       % Load Calcium data (extract it from the concatenated file if it wasn't done before)

SizeData = size(results.C, 2); Experiment.SizeData = SizeData;                                     % Build false Experiment Struct.,
nNeurons = size(results.C, 1); Experiment.nNeurons = nNeurons;
Neurons = struct('Rises',struct('Matrix',cell(Experiment.nNeurons,1),'Starts',cell(Experiment.nNeurons,1),...   % and false Neurons Struct. that are required for GetRises
    'Ends',cell(Experiment.nNeurons,1),'Duration',cell(Experiment.nNeurons,1)));
Act=ones(nNeurons,SizeData);

neuron = 1;

dataC = results.C_raw(neuron,:);
fdataC2 =  dataC; %medfilt1(dataC,7); %GCaMP8 10; GCaMP6 :
[~, noise2] = estimate_baseline_noise2(fdataC2);
evts2=fdataC2>1.5*noise2;   %2*  1.5*                                       % Find region where there are events according to CNMFE (this is the most accurate method to find small events until now)
gfdataC2 = gradient(fdataC2);
gfdataC2(2:end) = gfdataC2(1:end-1);
noiseD2 = std(gfdataC2);
Rises2 = Act(neuron,:);
evts3 = gfdataC2>1.5*noiseD2;
Rises2(gfdataC2<1*noiseD2)=0; %2* ; Rises2(gfdataC2<2*noiseD2)=0;
%evts = evts2;
evts = and(evts2, evts3);%
Rises2(not(evts)) = 0;
dRises2=diff(Rises2>0); starts=find(dRises2==1)+1; ends=find(dRises2==-1);
if starts(1)>ends(1); starts=[1 starts]; end                    % Boundaries checks
if starts(end)>ends(end); ends=[ends SizeData]; end

% calcul de la position de l'evenenemnt comme le max de la derivee dans la
% fenetre de rises;



f=figure;
time = 0: 0.01 : (size(results.C, 2)-1)*0.01; % time = Experiment.timeFluo-Experiment.timeFluo(1);
ax1 = subplot(4,1,1);
plot(time, results.C_raw(neuron, :)); hold on; plot(time, results.S(neuron, :)); plot(time,  fdataC2); % plot(time, ftdataC(:,2));
x = [time(1) time(end)]; y = [0 0]; line(x,y,'Color','black','LineStyle','--')
ax2 = subplot(4,1,2);
plot(time, Rises2); ylim([-0.2 1.2]); %hold on; plot(time, Act(neuron,:));
ax3 = subplot(4,1,3);
plot(time, gfdataC2); %plot(time, gdtdataC); hold on;
ax4 = subplot(4,1,4);
plot(time, evts2); hold on; plot(time, evts3); plot(time, evts); ylim([-0.2 1.2]); %hold on;  plot(time,1+ evts1); plot(time, 2+evts2);
%plot(time, Neurons.Rises(neuron).Matrix);    hold  on;  plot(time, results.S(neuron,:)); %plot(time, SfmaFilt(:,2)); plot(time, Sthresh) % ylim([-0.2 1.2]) % plot(time, evts); %plot(time, s); % plot(evts);  plot(dC); plot(sdFiltAct); plot(RisesIdx); %plot(results.S(neuron,:));

%                 toPlot1=Experiment.positionXcmSmooth; toPlot1(setdiff(1:end,find(Experiment.Dir1.Segments)))=nan;
%                 toPlot2=Experiment.positionXcmSmooth; toPlot2(setdiff(1:end,find(Experiment.Dir2.Segments)))=nan;
%                 plot(time, Experiment.positionXcmSmooth); hold on; plot(time, toPlot1); plot(time,toPlot2)
linkaxes([ax1,ax2, ax3, ax4],'x');


% align all the events for this specific neuron;
figure; 
%hold on; 
for k = 1 : numel(starts)
    k
    ax1 = subplot(2,1,1);
    plot(-200:10:+500, results.C_raw(neuron, starts(k)-20:starts(k)+50)-results.C_raw(neuron, starts(k)));
    ax2 = subplot(2,1,2);
    plot(-200:10:+500, gfdataC2(starts(k)-20:starts(k)+50));
    linkaxes([ax1,ax2],'x');
    pause; 
end

