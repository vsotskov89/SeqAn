function timestamps = Unshift(timestamps, intervals),

% The opposite of the option 'shift' in Restrict

shiftedIntervals(:,2) = CumSum(diff(intervals,[],2));
shiftedIntervals(:,1) = [0; shiftedIntervals(1:end-1,2)];
toAdd = intervals(:,1) - shiftedIntervals(:,1);
[ok,w] = InIntervals(timestamps,shiftedIntervals); 
timestamps(~ok) = nan;
timestamps(ok) = timestamps(ok) + toAdd(w(ok));
