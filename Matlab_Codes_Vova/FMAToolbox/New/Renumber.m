function x = Renumber(x,mode)

%
%   This function takes a vector and renumbers all the entries so there are
%   no gaps between them and they start from one.
%   
%   EXAMPLE: 
%   x = [30 4 4 5 30 6];
%   x = Renumber(x);
%   
%   This returns [1 2 2 3 1 4].
%

%% Checking if numbers are integers as they should: 

if nansum(rem(x,1)) ~= 0,
    error('Input should be an array of integers');
end
if ~exist('mode','var')
    mode = 'stable';
end


%% Going through the list changing the numbers

% [uniquelist index] = unique(x, mode);
% [~, index] = sort(index);
% uniquelist = uniquelist(index); %the unique numbers in the order that they appear in
% uniquelist = uniquelist(~isnan(uniquelist));
% howmany = numel(uniquelist);
% 
% for i=1:howmany,
%     x(x==uniquelist(i)) = -i;
% end
% x = -x;

[~,~,x] = unique(x,mode);
