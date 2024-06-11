function r = RandomIntervals(n, sumDuration)

% generates n random intervals from 0 to 1, covering a total of 'sumDuration'

if ~InIntervals(sumDuration,[0 1]),
    error('Overall duration of the intervals should be between 0 and 1');
end
r = zeros(n,2);
r(:,1) = sortrows(rand(n,1));
r(:,2) = r(:,1)+sumDuration.*([r(2:end,1); 1]-r(:,1)); % on avarage, that's 0.1
r(1) = r(1)-sumDuration*(r(1));