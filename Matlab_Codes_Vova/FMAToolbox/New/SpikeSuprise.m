function [z,p] = SpikeSuprise(count, firingRates, binSize)
%SpikeSuprise Estimates the surprise of observing a number of spikes ('count')
% counted over a wondow of a given duration ('binSize')
% given a neuron's average firing rate ('firingRates')
% assuming Poisson spiking
% The surprise is given in z units where z>1.96
% corresponds to spike counts significantly higher than
% what would be expected from a Poisson neuron, and
% z<-1.96 corresponds to significantly lower spike
% counts than expected.
% 'count' and 'firingRates' can be vectors corresponding to
% multiple neurons

r = firingRates;
k = count;
dt = binSize;

% Code for p(x==k), later replaced by code for p(x>k)
% p = exp(k.*log(r.*dt)-logfactorial(k)-r.*dt); 
% smaller = k<r.*dt; % for the sign of z

lambda = r.*dt;
p = poisscdf(k,lambda);
smaller = p<0.5;

p(~smaller) = poisscdf(k(~smaller),lambda(~smaller),'upper');

z = p2z(p*2); % make it two-tailed. p=0.98 is just as unlikely as p=0.02;
z(abs(z)==Inf) = 40; % matlab rounds p-values smaller than z=40 as 0, i.e. z=Inf
z(smaller) = -z(smaller);
