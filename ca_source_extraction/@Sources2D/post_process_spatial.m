function A_ = post_process_spatial(obj, A_new)
% postprocess spatial components of all neurons 

if ~exist('A_new', 'var')
    A_new = obj.A;
end

%             A_new = threshold_components(A_new, obj.options);
spatial_constraints = obj.options.spatial_constraints;
circular_shape = spatial_constraints.circular;
connected_shape = spatial_constraints.connected;
if ismatrix(A_new)
    d1 = obj.options.d1;
    d2 = obj.options.d2;
    d = d1*d2;
    K = size(A_new,2);
else
    [d1, d2, K] = size(A_new);
    d = d1*d2;
    A_new = reshape(A_new, d, K);
end

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