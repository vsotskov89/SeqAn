function logical = IntervalsIntersect(intervals1, intervals2)

% In case you have, say, windows (intervals1), and events with a certain 
% duration, say, ripples (intervals2), and you want to have the identity of
% the windows which contain ripples, you may not want to limit yourself to
% the windows containing the ripple peak.
% This function returns a logical vector the size of intervals1, with ones 
% if any portion of intervals2 is present in the respective interval1.


% % start of intervals1 is included within an interval2. This should be counted.
% l1 = InIntervals(intervals1(:,1), intervals2);
% 
% % end of intervals1 is included within an interval2. This should be counted.
% l2 = InIntervals(intervals1(:,2), intervals2);
% 
% % intervals2 is completely contained within an interval1. This is the last option that should be counted. 
% % Note that InIntervals(intervals2(:,2), intervals1) would not bring more cases.
% l3 = CountInIntervals(intervals2(:,1), intervals1)>0;



% start or end of intervals2 is included within an interval1. This should be counted.
l1 = CountInIntervals(intervals2(:,1), intervals1)>0;
l2 = CountInIntervals(intervals2(:,2), intervals1)>0;
% intervals2 completely incorporates ones of intervals1.
l3 = InIntervals(intervals1(:,2), intervals2);

logical = l1 | l2 | l3;