function r = RelativeTimes(t,intervals,values),

% This function calls RelativeTime, but it takes in multiple-column intervals. 
% Provide values for each column of the interval (if t coincides with a an element
% of, say, the 3rd column, it takes the third value from the list 'values'). If
% values are not provided, the default values are from 0 to 1 for timestamps between
% the first and second columns, between 1 and 2 for timestapts between the second
% and the third columns, etc.
% Example, 
% RelativeTimes(t, [spindles(:,1)-0.5 spindles(:,[1 2 3]) spindles(:,3)+0.5], [-2 -1 0 1 2])
% will return a vector depicting where each timestamp was found in the sequence of
% events: 0.5s before the spindle starts, spindlestart,spindlepeak,spindlestop,0.5s after.
% i.e. a value of 0.25 will reflect the timestamp is just after a spindle peak, 25% of the
% between the spindle peak and the spindle end.

if nargin<3,
    values = 0:size(intervals,2)-1;
end

r = nan(length(t),size(intervals,2)-1);
for i=1:size(intervals,2)-1.
    r(:,i) = RelativeTime(t,intervals(:,[i i+1]));
end

% In case a timestamp is present in more than one event, we will select the
% event that is closest to the timestamp, and ignore the timestamps's position relative to
% other events. (In the above example, if a timestamp is 0.2 seconds after a spindle end, 
% but also 0.1 seconds before a spindle start, the 0.2 value will be ignored).
[~,columnsToKeep] = max(abs(0.5-r),[],2);
indicesToKeep = sub2ind(size(r),(1:length(t))',columnsToKeep);

for i=1:size(intervals,2)-1,
    r(:,i) = r(:,i)*(values(i+1)-values(i)) + values(i);
end

r = r(indicesToKeep);