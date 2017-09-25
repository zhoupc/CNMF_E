function A_ = post_process_spatial(obj, A_new)
% postprocess spatial components of all neurons

if ~exist('A_new', 'var')
    A_new = obj.reshape(obj.A, 2);
end

%             A_new = threshold_components(A_new, obj.options);
spatial_constraints = obj.options.spatial_constraints;
circular_shape = spatial_constraints.circular;
connected_shape = spatial_constraints.connected;

[d1, d2, K] = size(A_new);
d = d1*d2;
A_new = reshape(A_new, d, K);

A_new = mat2cell(A_new, d, ones(1,K));

parfor m=1:K
    ai = reshape(full(A_new{m}), d1, d2);
    
    % remove isolated pixels
    if connected_shape
        ai = connectivity_constraint(ai);
    end
    
    if circular_shape
        ai = circular_constraints(ai); % assume neuron shapes are spatially convex
    end
    
    A_new{m} = ai(:);
end
A_ = sparse(cell2mat(A_new));
end