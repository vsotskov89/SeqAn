function [RipSceCCG,tCCG]=GetRipSceCCG(ripples,putSCEtimes,durations,binSize)
StartOrPeak=2; % [start_t peak_t end_t peakNormalizedPower] 1: Start 2: Peak
dur=diff(durations);
RipAndSCEs=[[ripples(:,StartOrPeak) ones(size(ripples,1),1)];[putSCEtimes 2*ones(size(putSCEtimes,1),1)]];
[RipSceCCG,tCCG] = CCG(RipAndSCEs(:,1),RipAndSCEs(:,2),'duration',dur,'binSize',binSize);

        
        
        
        
        
        
        