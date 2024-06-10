function [Experiment,results,Params]=GetPositionSpeed(Params,Experiment,results)
% Get instant position and speed of the animal %

% Load Position Data
try
    load(Experiment.PositionFile);                                                                % Position variable name is mousePos1
    positionX=mousePos1(1,:);                                                            % We are only interested in mouvement in the X direction
catch
    disp('  Unable to recover Position file : Manual selection required')
    tmp=[]; PosFiles=[];
    while ~isequal(tmp,0)
        tmp=uigetfile('*.csv','Select All Position files',Experiment.path,'MultiSelect','on');
        PosFiles=[PosFiles {tmp}];
    end
    mousePos1=[];
    for file=1:numel(PosFiles)-1
        tmp=csvread([Experiment.path PosFiles{file}],3);
        mousePos1=[mousePos1 tmp(:,2:3)'];
    end
    if min(mousePos1(1,:))<=0; mousePos1(1,:)=mousePos1(1,:)+2*abs(floor(min(mousePos1(1,:)))); end
    % sometimes DLC gets the position wrong and attributes a position at
    % the begining of the track. Let's remove those aberrant points
    %figure; plot(mousePos1(1,:)); hold on; plot(gradient(mousePos1(1,:))); 
    mousePos1x = mousePos1(1,:);
    mousePos1x(mousePos1x<20) = NaN; %20
%    mousePos1x(mousePos1x>1235) = NaN; 
    mouseSpeed1x = gradient(mousePos1x);
    artefacts = mouseSpeed1x' - smooth(mouseSpeed1x', 'rlowess');
  
    figure; 
    ax1 = subplot(2,1,1);
    plot(mousePos1x);
    ax2 = subplot(2,1,2);
    plot(gradient(mousePos1x));  hold on; plot(artefacts); 
    linkaxes([ax1,ax2],'x');
    
    % with DLC there point that are mislabeled. We try an easy solution :
    % if the speed is too high, we remove points.  
    speedTh = 40;
    k1  = 1;
    kseg = 7; 
    while k1 < size(mousePos1, 2)-kseg -1
        %k1 
        k2max = 0;
        if abs(artefacts(k1))>speedTh
            s = sign(artefacts(k1));
            mousePos1x(k1+1) = NaN;
            for k2= 1:2*kseg
                if s*artefacts(k1+k2)<-speedTh
                    mousePos1x(k1+2 : k1+k2+1) = NaN(1, k2);
                    k2max = k2;
                    %break;
                end
            end
            if k2max == 0
                    for k2=1:kseg
                        if artefacts(k1+k2)>speedTh
                            mousePos1x(k1+2 : k1+k2+1) = NaN(1, k2);
                        end
                    end
            end
        end
        k1 = k1 + k2max + 1;
    end
    mousePos1x = fillmissing(mousePos1x,'linear','EndValues','nearest');
 
    figure; 
    ax1 = subplot(2,1,1);
    plot(mousePos1x);
    ax2 = subplot(2,1,2);
    plot(gradient(mousePos1x)); % hold on; plot(artefacts); 
    linkaxes([ax1,ax2],'x');

    mousePos1(1,:) = mousePos1x;
    save([Experiment.path 'mousePos1.mat'], 'mousePos1')
    positionX=mousePos1(1,:);
end

positionBaslerPix = positionX;
%positionX = smooth(positionX,7,'sgolay');              %smoothing
%positionX= interp1(Experiment.timeBasler,positionX,Experiment.timeFluo,'spline') ; % Rescale Behavioral Data to Fluo Data (20 points/s --> 100points/s)

positionX= interp1(Experiment.timeBasler,positionX,Experiment.timeFluo) ; % Rescale Behavioral Data to Fluo Data (20 points/s --> 100points/s)

% Remove some frames
Params.nStart = 5;                                                                      % Remove the 4 first frames (that were removed for CNMFE analysis)
Experiment.timeFluo = Experiment.timeFluo(Params.nStart : end);                         % OLD check : timeFluo = timeFluo(nStart :find(timeFluo==max(timeFluo(:))));
positionX = positionX(Params.nStart : end);                       % OLD check : positionX = positionX(nStart :find(timeFluo==max(timeFluo(:))));
% positionX = fillmissing(positionX, 'linear');                     % Fill any missing data
if Params.test10
    Experiment.timeFluo=Experiment.timeFluo(1:10:numel(Experiment.timeFluo));
    positionX = positionX(1:10:numel(positionX));
    c_raw=zeros(size(results.C_raw,1),size(Experiment.timeFluo,1));
    c=zeros(size(results.C,1),size(Experiment.timeFluo,1));
    for r= 1:size(results.C_raw,1)
        c_raw(r,:)=interp1(1:numel(results.C_raw(1,:)),results.C_raw(r,:),1:numel(positionX));
        c(r,:)=interp1(1:numel(results.C(1,:)),results.C(r,:),1:numel(positionX));
    end
    results.C_raw=c_raw;
    results.C=c;
    Experiment.SizeData=size(Experiment.timeFluo, 1); 
end
%Convert position from px to cm
f=figure('Visible','off');
    binWidth = 10; 
    h = histogram(positionX,'BinWidth',binWidth); 
    threshold =floor(((max(positionX)-min(positionX))/2+min(positionX))/binWidth);          % ~Middle of the track
    
    [fit_norm, ~]=fit((1:threshold)',(h.Values(1:threshold))','gauss1');                    % Fit au gaussian on the 2 halves of the track
    posReward1 = fit_norm.b1*binWidth +  h.BinEdges(1);                                     % position of the reward is assumed to be the mean of the gaussian
    [fit_norm, ~]=fit((threshold:size(h.Values, 2))',(h.Values(threshold:end))','gauss1');
    posReward2 = fit_norm.b1*binWidth +  h.BinEdges(1);
close(f)
Experiment.px2cm = Params.rewardDist/(posReward2-posReward1);
positionXcm = positionX * Experiment.px2cm;

Experiment.positionXcmSmooth = abs(smooth(positionXcm,'sgolay'));              %smoothing
%Experiment.positionXcmSmooth = abs(smooth(positionXcm,5, 'lowess'));              %smoothing

% Compute speed
% deltaT = (Experiment.timeFluo(1000)-Experiment.timeFluo(1))/1000;          % dt ; % OLD dt : deltaT = (timeFluo(end)-timeFluo(1))/SizeData;               
% Experiment.mouseSpeed = gradient(Experiment.positionXcmSmooth)/deltaT;     % v=dx/dt --> speed can be negative (dir. 2)
%Experiment.mouseSpeed = smooth(Experiment.mouseSpeed,'sgolay');              %smoothing

% Compute speed from initial trace

positionBaslerCmSmooth = abs(smooth(positionBaslerPix, 5, 'lowess'))* Experiment.px2cm;              %smoothing and conversion to cm
%figure; plot(positionBaslerPix); hold on; plot(positionBaslerPixSmooth);
deltaTBasler = (Experiment.timeBasler(1000)-Experiment.timeBasler(1))/1000;          % dt ; % OLD dt : deltaT = (timeFluo(end)-timeFluo(1))/SizeData;               
mouseSpeedBasler = gradient(positionBaslerCmSmooth)/deltaTBasler;
mouseSpeedFluo= interp1(Experiment.timeBasler,mouseSpeedBasler,Experiment.timeFluo) ; % Rescale Behavioral Data to Fluo Data (20 points/s --> 100points/s)
Experiment.mouseSpeed = mouseSpeedFluo(Params.nStart : end);  
Experiment.positionXcmSmooth = Experiment.positionXcmSmooth(Params.nStart : end);  % added by cathie on 14/12/2022 
 
% figure; 
% ax1 = subplot(2,1,1);
% plot(Experiment.positionXcmSmooth);
% ax2 = subplot(2,1,2);
% plot(Experiment.mouseSpeed); hold on; yline(2); yline(-2); 
% linkaxes([ax1,ax2],'x');



