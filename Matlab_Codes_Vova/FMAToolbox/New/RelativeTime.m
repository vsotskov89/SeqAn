function [RT, whichIntervals] = RelativeTime(samples, intervals)

% Returns a vector of the relative time within an interval that a sample is
% in. For example, a sample with a timestamp in the middle of an interval
% would have relative time 0.5. A sample that falls outside of the
% intervals given is NaN.

samples = samples(:,1); %only timestamps are taken


RT = nan(size(samples));

if size(intervals,2)==3,
    whichIntervals = zeros(length(samples), size(intervals,2)-1);
    for i=1:size(intervals,2)-1,
        theseIntervals = intervals(:,i:i+1);
        [inIntervals,whichIntervals(:,i)] = InIntervals(samples, theseIntervals);
        RT(inIntervals) = RelativeTime(samples(inIntervals), theseIntervals)+i-1;
    end
    return
end
[inIntervals whichIntervals] = InIntervals(samples, intervals);
RT(inIntervals,1) = (samples(inIntervals) - intervals(whichIntervals(inIntervals),1))./...
    (intervals(whichIntervals(inIntervals),end) - intervals(whichIntervals(inIntervals),1));