function [CaLFP,LFP]=AlignLFP(results,Params,Experiment)

%% Loading data and finding timings
CaTTL=loadChanDat(Params.LFP.TTLchan,Experiment.path); CaTTL = single(CaTTL);   % Load the TTL data
dTTL = diff(CaTTL);                                                             % Time derivation of the TTL signal to catch TTL rises

dTTL(dTTL<max(dTTL)/50) = 0;                                                    % Set to 0 everything below a half of the max
TTLstarts = find(dTTL);                                                         % Find the position of every remaining non zero values (strong rises)

distStarts = diff(TTLstarts);                               % find distance between putative rises (should be around 200 because TTL is 100Hz so 1 every 10ms=200pts at 20kHz)
TTLstarts(distStarts < 0.5*median(distStarts)) = [];        % remove points that are too close (<medianDist/2 ?) (Happen when rises is on 3 points --> 2 points kept)

distStarts = diff(TTLstarts);                               % Update distances
TTLBreaks = find(distStarts > 3*median(distStarts));        % Find gap between TTLs, i.e when TTLs are spaced by more than 3 times the median
TTLBreaks=[TTLBreaks(1) TTLBreaks(end)];                    % Keep the first and the last value (just in case there is a gap in the middle, hopefully there is no gap before the start and after the end

TTLstarts((TTLBreaks(2)+1):end) = [];                       % remove TTL after second break (so after the end of the experiment)
TTLstarts(1:TTLBreaks(1)) = [];                             % remove TTL before first break (so before the start of the experiment)

if Params.LFP.CheckPlot                                     % Plot to check if Experiments start/stop are well placed on the TTL recording
    TTLBreaksTimes = [(TTLstarts(TTLBreaks(1))+TTLstarts(TTLBreaks(1)+1))/2 (TTLstarts(TTLBreaks(2))+TTLstarts(TTLBreaks(2)+1))/2]; % Timings of the breaks
    figure; plot(CaTTL); hold on; yl = [min(CaTTL) max(CaTTL)];
        line([TTLBreaksTimes(1) TTLBreaksTimes(1)],yl,'Color','red'); line([TTLBreaksTimes(2) TTLBreaksTimes(2)],yl,'Color','red');
end
clear CaTTL dTTL                                            % Clear heavy variables

%% Upsample Calcium to 1250Hz
disp(['4 firsts TTL were removed ; ' num2str(numel(TTLstarts)-4-numel(results.C_raw(1,:))) ' were removed at the end'])
TTLstarts=TTLstarts(5 : 4+numel(results.C_raw(1,:)));                                   % Remove the 4 firsts TTLs (these images are not kept and keep the timings for calcium data points
TTLstartsLFP=round(TTLstarts/16);                                                       % timings 20kHz --> timings 1250Hz

timingsLFP=TTLstartsLFP(1):TTLstartsLFP(end);                                           % range of timings at 1250Hz    

%% Load LFP and filter and get theta phase

LFP=loadChanLFP(Params.LFP.chanLFP,Experiment.path);                       % Load LFP
if Params.LFP.CheckPlot                                         % Plot to check if how calcium and LFP are aligned
    upCaTrace=interp1(TTLstartsLFP,results.C_raw(1,:),timingsLFP);  % interpolate 1 Calcium traces based on LFP timings of TTL to 1250Hz
    figure   
        subplot(2,1,1)                                               
            plot(LFP(:,2)./std(LFP(:,2))); hold on; plot(timingsLFP,upCaTrace./std(upCaTrace))
        subplot(2,1,2)
            plot(LFPrestrict(:,2)./std(LFPrestrict(:,2))); hold on; plot(upCaTrace./std(upCaTrace))
end 

FiltLFP = FilterLFP(LFP,'passband','theta');
[ThetaPhase,~,~] = Phase(FiltLFP);
CaThetaPhase=ThetaPhase(TTLstartsLFP,:);
CaLFP=LFP(TTLstartsLFP,:);
% clear ThetaPhase LFP FiltLFP timingsLFP                                             % Clear heavy variables
%% DO RECONSTRUCTION HERE

%% 
% CaThetaPhase(:,1)=CaThetaPhase(:,1)-CaThetaPhase(1,1)+0.01;
% CaThetaPhase=Interpolate(CaThetaPhase,stats.windows(:,1));
% 
% %% Fast plot
% err=reconstructPos-stats.positions(:,2)';
% figure
% scatter(CaThetaPhase(:,2),err)
% 

%plot(uw(1:end-1,1),diff(smooth(uw(:,2),100))./(diff(uw(:,1))*2*pi))












