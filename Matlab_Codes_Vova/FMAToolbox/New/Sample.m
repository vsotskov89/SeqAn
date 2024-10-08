function samples = Sample(histogram, n)

% sample from distributions.
% EXAMPLE
% [h,t] = hist(spikes(:,1),resolution);
% samplespikes = t(Sample(h,1000));
% NOTE: multiple histograms may be provided in different columns of h.
%
% By Raly Todorova (2016), inspired by Dahua Lin's discretesample (2008)

if nargin<2,
    error('Please provide 2 arguments: pdf(s)/histogram(s) and n of desired samples');
end

if isvector(histogram),
    histogram = histogram(:);
end

nDistributions = size(histogram,2);
nBins = size(histogram,1);

histogram = bsxfun(@rdivide,histogram,nansum(histogram));

% The main reasoning of this code was inspired by Dahua Lin's discretesample (2008)
% The reasoning is as follows:
% We can sample from a uniform distribution using rand and then check where it falls relative to some intervals
% We have nBins intervals. If rand(i) falls within the 16th bin, sample(i)=16.
% Now, to adjust the probability of obtaining each possible value, we can adjust the ineterval durations to reflect out distribution:
% The i-th interval duration should be proportional to the pdf (histogram) value at the i-th bin provided.
% This is achieved by creating intervals that start at 0 and build up with the values of the pdf (histogram):
cdf = [zeros(1,nDistributions); cumsum(histogram)];
intervals = [reshape(cdf(1:end-1,:),[],1) reshape(cdf(2:end,:),[],1)];
% the intervals need to be sorted, or else CountInIntervals skips them
[intervals,order] = sortrows(intervals); 
count = CountInIntervals(sort(rand(n,1)),intervals);
% get a translation vector that would reverse the ordering of the intervals
inverse(order,1) = 1:length(order);
count = count(inverse,1);

% The generalisation to multiple distributions is achieved by having them all stacked below one another counts = [countsFor1; countsFor2;...;countsForLast]; 
% Sampling is achieved by repeating the i-th id count(i) times (e.g. rand(i) falls within the 16th bin => count(16) is increased by 1, and 16 is sampled
samples = repelem(rem((1:nBins*nDistributions)-1,nBins)'+1,count);
% Reshape the final result so that the samples stacked below one another are ordered by respective distributions (one per column)
samples = reshape(samples,n,nDistributions);
samples = Shuffle(samples')';

end

%% HELPER FUNCTIONS

function shuffled = Shuffle(matrix,Logical)

% row by row
% if logical is provided, it should indicate which elements are to shuffle (0s will remain as given)

if nargin>1,
    [~,indices] = sort(rand(size(matrix')).*(-Logical'));
else
    [~,indices] = sort(rand(size(matrix')));
end
shuffled = matrix(sub2ind(size(matrix),meshgrid(1:size(matrix,1),1:size(matrix,2))',indices'));
end


function out = repelem(vector,repetitions)

% repeat each element separately [matlab 2015+ should have a function 'repelem' that does this]
% EXAMPLE: repelem([5;2],[3;4]) = [5;5;5;2;2;2;2];
% Only vector input is supported.
% If you need to perform the action with a matrix, use repelem on indices:
% ind = 1:size(matrix,1);
% repeated = matrix(repelem(ind(:),repetitions),:);

% In case a single value was provided as repetitions, this value applies for all the elements of 'vector'
if length(repetitions)==1,
    repetitions = repetitions*ones(size(vector));
end

% This function works with column inputs
repetitions = repetitions(:);
vector = vector(:);

% Prepare the output vector with as many elements as the total of the counts (repetitions) expected
out = zeros(sum(repetitions),1);

% For each element of 'vector', get the first respective index of 'out' that should be assigned this element
% Example: vector = [10 20 30]; repetitions = [1 2 1]; starts = [1 2 4], 
% as 10 should be first seen in the 1st index of 'out', 20 in the 2nd, and 30 in the 4th.
% This is computed the following way: the first index we observe any value is, by default, 1. The following index
% starts(2)  = (1 + the number of repetitions of vector(1)). The following one, 
% starts(3) = (starts(2) + the number of repetitions of vector(2)), and so on: 
starts = cumsum([1; repetitions(1:end-1)]);

% Assign the values of 'vector' to their first respective indices
out(starts) = vector;

reset = false(size(out));
reset(starts) = true;

% The vector 'out' contains zeros where the previous value should be repeated.
out = CumSum(out,reset);

% Reasoning as to why the function works when certain 'repetition' values are set to 0.
% This function is centered around a curious behaivour of Matlab:
% when assigning multiple values of the same index, matlab only assigns the last value given.
% Example: a = zeros(1,2); a([1 1 1 2 2]) = [1 2 3 4 5]; gives a= 3 5. 
% (Matlab overrides a(1)=1 and a(1)=2 by a(1)=3 as it comes last. Same for a(2)=4 overriden by a(2)=5)
% We can use this as an advantage when certain values of 'repetitions' are zero:
% When, say, repetitions(1) = 0, starts = [1;1...], and the first index of out being 
% set to the first value of 'vector' (i.e. out(1) = vector(1)), is overriden by the next assignment
% (out(1) = vector(2)). Only the last unique value of 'starts' remains valid, so all the
% assignements of 'out' by 'vector' indices of 0 repetition will be overriden.
% The only exception to this is the last index of 'vector', i.e. out(n) = vector(end);
% There is no future index to override this, therefore a solution is to check for this problem and 
% remove this aberrant value of 'out':
if repetitions(end) == 0, out(end) = []; end

end