function ppc = PPC(phases)


% Pairwise Phase Consistency as defined by  
% Vinck, Wingerden, Womelsdorf, Fries, Pennartz (2010, NeuroImage)

phases = phases(:);
difference = bsxfun(@minus,phases,phases');
difference(eye(size(difference))==1) = nan; % exclude angle comparisons with itself
% as the diagonal is now composed of nans, there are n*(n-1) non-nan elements in the matrix 'difference' 
% i.e. taking the mean is equilalent to the sum, normalised by n*(n-1) as in Vinck et al (2010)
ppc = nanmean(cos(difference(:))); % as cos(a-b) = cos(a)*cos(b)+sin(a)*sin(b), Vinck et al's definition of PPC
