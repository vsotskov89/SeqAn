function [b, sn] = estimate_baseline_noise2(y)

% estimate the baseline and noise level for single calcium trace 
% use only left part of histogram
% inputs: 
%   y:  T*1 vector, calcium trace 
 
% outputs: 
%   b: scalar, baseline 
%   sn: scalar, sigma of the noise 

temp = quantile(y, 0:0.1:1); 
dbin = (max(temp)-min(temp))/1000; 
f=figure('Visible','off');
h = histogram(y,'BinWidth',dbin);
threshold = find(h.Values==max(h.Values))+3 ; % on prend le max de l'histogram auquel on rajoute 3 bin pour avoir bien le pic entier; 3 est un peu arbitraire, on affine plus bas. 
[fit_norm, ~]=fit((1:threshold)',(h.Values(1:threshold))','gauss1');
% General model Gauss1:
% fit_norm(x) =  a1*exp(-((x-b1)/c1)^2) 
% sn = fit_norm.c1 *dbin /sqrt(2);
% b = fit_norm.b1 * dbin + h.BinLimits(1);

threshold = find(h.Values==max(h.Values))+ floor(fit_norm.c1/2) ; 
[fit_norm, ~]=fit((1:threshold)',(h.Values(1:threshold))','gauss1');
% General model Gauss1:
% fit_norm(x) =  a1*exp(-((x-b1)/c1)^2) 
sn = fit_norm.c1 *dbin /sqrt(2);
b = fit_norm.b1 * dbin + h.BinLimits(1);
close(f)



