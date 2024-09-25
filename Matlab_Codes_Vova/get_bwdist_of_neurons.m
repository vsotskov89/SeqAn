function bw_dists = get_bwdist_of_neurons(A, dim1, dim2)
%Get the distribution of half-"thickness" of neurons which is the maximum
%of the euclidean distance transform
As = reshape(full(A'), size(A,2), dim1, dim2);
bw_dists = zeros(size(As,1),1);
for i = 1:size(As,1)
    bw_dists(i) = max(max(bwdist(~squeeze(As(i,:,:)) > 0)));
end

