function filter_neurons_by_bwdist(neuron, min_bwdist)
%filter neurons by minimal value of their half-"thickness" which is the maximum
%of the euclidean distance transform
    tags_ = zeros(size(neuron.A,2), 1, 'like', uint16(0));
    bw_dists = get_bwdist_of_neurons(neuron.A, neuron.options.d1, neuron.options.d2);
    tags_ = tags_ + uint16(bw_dists < min_bwdist);
    ids = find(tags_);
    if ~isempty(ids)
        neuron.delete(ids);
    end
end