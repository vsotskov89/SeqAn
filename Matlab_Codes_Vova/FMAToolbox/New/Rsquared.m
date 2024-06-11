function [rsq, SSres] = Rsquared(model, data);

% Example rsq = Rsquared(polyfit(data(:,1), data(:,2),1), data(:,[1 2]);

yfit = polyval(model, data(:,1));
yresid = data(:,2) - yfit;
SSres = sum(yresid.^2);
SStotal = (length(data(:,2))-1)*var(data(:,2));
rsq = 1-SSres/SStotal;