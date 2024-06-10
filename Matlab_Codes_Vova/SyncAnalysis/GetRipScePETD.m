function GetRipScePETD(ripplesPower,putSCEtimes,durations)

[RipPowsync,indices] = Sync(ripplesPower,putSCEtimes,'durations',durations);
[d,t,p] = SyncHist(RipPowsync,indices,'mode','dist','durations',durations,'smooth',1);%,'bins',[0 4000 200]);
    figure;
        PlotColorMap(d,1); clim([0 0.01]); colormap parula ;clim([0.00 0.007]) %clim([0.001 0.0075])
        xticks(1:10:100);xticklabels(round(t(1:10:100),1)); xlabel('Time centered on SCE (s)'); 
        yticks(1:10:100);yticklabels({round(p(1:10:100))}); ylabel('Power in the Ripples band (ua)')
        title('Peri-SCEs time distribution of Ripples band power')