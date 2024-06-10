function normalizedSquaredSignal=testGetNSS(signal)

frequency = 1250;

windowLength = round(frequency/1250*11);
window = ones(windowLength,1)/windowLength;
keep = [];
sd = [];

squaredSignal = signal.^2;

[normalizedSquaredSignal,~] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);



function y = Filter0(b,x)

if size(x,1) == 1,
	x = x(:);
end

if mod(length(b),2) ~= 1,
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;
