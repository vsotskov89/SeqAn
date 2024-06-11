function cMap=  ColorMap(axis,varargin)

% This function changes the color map of the given axis to a custom color map,
% a gradient of the desired colors. The desired colors should be 2 at the minimum 
% [start stop], and can be any number above that.
% EXAMPLE:
% ColorMap(gca,[1 1 1],[1 0 0]); % changes the colormap to a gradient from white ([1 1 1]) to red ([1 0 0])
% ColorMap(gca,[0 0 0.2],[0 0.5 0.3],[1 1 0]); % will change the colormap to a map gradiating from dark blue to dark green to yellow

nMilestones = length(varargin);
if nMilestones<2,
    error('Please provide at least a starting and an ending color');
end


for i=1:nMilestones,
    if ~isvector(varargin{i}) || length(varargin{i})~=3 || sum(varargin{i}>1 | varargin{i}<0)>0,
        error('arguments should be in number format, e.g. [1 0 0] for red');
    else
        colors(i,:) = varargin{i};
    end
end

cMap = nan(nMilestones,3);
cMap(1,:) = colors(1,:);
cMap(end,:) = colors(nMilestones,:);

milestones = [1 round(([1:nMilestones-1])*72/(nMilestones-1))];

try
    for i=2:length(milestones),
        cMap(milestones(i),:) = colors(i,:);
        for column=1:3,
            if cMap(milestones(i),column)==cMap(milestones(i-1),column),
                cMap((milestones(i-1):milestones(i)),column) = cMap(milestones(i-1),column);
            else
                cMap((milestones(i-1):milestones(i)),column) = cMap(milestones(i-1),column):(cMap(milestones(i),column)-cMap(milestones(i-1),column))/(milestones(i)-milestones(i-1)):cMap(milestones(i),column);
            end
        end
    end
catch
    keyboard
end

cMap(cMap>1) = 1;cMap(cMap<0) = 0;

colormap(axis,cMap);

