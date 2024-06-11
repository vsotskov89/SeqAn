function smoothed = Shrink(matrix,rowBinSize,columnBinSize)

% This is a very fast version of Smooth. Instead of applying a kernel, this
% function would return a matrix 'smoothed' which is produced by taking the
% mean of [columnBinSize rowBinSize]-sized blocks within the matrix.
% From this follows that columnBinSize needs to be a factor of size(matrix,2),
% and rowBinSize needs to be a factor of size(matrix,1). If not, nans will be
% be appended (as equally spaced as possible).
%%
n0 = size(matrix);
n = ceil(n0./[rowBinSize columnBinSize]).*[rowBinSize columnBinSize];

% This is a heuristic to make the code run faster: 
if columnBinSize>rowBinSize, 
    smoothed = Shrink(matrix',columnBinSize,rowBinSize)';
    return;
end

if n0(1)~=n(1),
    if n0(1)==1, matrix = [matrix; nan(n(1)-1,n0(2))]; else
    nanRows = round(linspace(1,n(1),2+n(1)-n0(1))); nanRows([1 end]) = []; % Equally spaced rows of nans
    matrix0 = matrix;
    matrix = nan(n(1),n0(2));
    matrix(~ismember(1:n(1),nanRows),:) = matrix0;
    end
end

if n0(2)~=n(2),
    if n0(2)==1, matrix = [matrix nan(n(1),n(2)-1)]; else
	nanColumns = round(linspace(1,n(2),2+n(2)-n0(2))); nanColumns([1 end]) = []; % Equally spaced columns of nans
    matrix0 = matrix;
    matrix = nan(n);
    matrix(:,~ismember(1:n(2),nanColumns),:) = matrix0; % fill non-nan values
   end
end


if columnBinSize == 1 && rowBinSize == 1,
    smoothed = matrix;
    return
end

order = [];
for i=1:columnBinSize,
    if i<columnBinSize,
        order = [order find(rem(1:size(matrix,2),columnBinSize)==i)];
    else
        order = [order find(rem(1:size(matrix,2),columnBinSize)==0)];
    end
end

% This code was produced by trial and error. Too complicated to understand, let alone annotate.
smoothed = matrix(:,order);
smoothed = reshape(smoothed,rowBinSize,[]);
smoothed = reshape(smoothed,[],columnBinSize)';
smoothed = reshape(smoothed,1,[])';
smoothed = reshape(smoothed,columnBinSize*rowBinSize,[]);
smoothed = nanmean(smoothed);
smoothed = smoothed(:);
smoothed = reshape(smoothed,size(matrix)./[rowBinSize columnBinSize]);

