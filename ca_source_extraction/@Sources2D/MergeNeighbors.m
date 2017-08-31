function  [merged_ROIs, newIDs] = MergeNeighbors(obj, dmin, method)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   dmin: scalar or 1*2 vector. If it's a scalar, it's the minimum distance
%   between two pixels; elsewise, it means [minimum distance, maximum
%   difference between two time constances of decaying]

% output:
%   merged_ROIs: cell arrarys, each element contains indices of merged
%   components
%   newIDs: vector, each element is the new index of the merged neurons

%% Author: Pengcheng Zhou, Carnegie Mellon University.
%  The basic idea is proposed by Eftychios A. Pnevmatikakis: high temporal
%  correlation + spatial overlap
%  reference: Pnevmatikakis et.al.(2016). Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data. Neuron

%% variables & parameters
A = obj.A;          % spatial components
if isempty(obj.C_raw)
    obj.C_raw = obj.C;
end
C_raw = obj.C_raw;
C = obj.C;
deconv_options_0 = obj.options.deconv_options;

%% find neuron pairs to merge
if ~exist('method', 'var') || isempty(method)
    method = 'center'; 
end
if strcmp(method, 'center')
    ctr = obj.estCenter();
    yy = ctr(:,1);
    xx = ctr(:,2);
else
    [~,temp] = max(obj.A, [], 1);
    [yy, xx] = ind2sub([obj.options.d1, obj.options.d2], temp);
end
dist_v = sqrt(bsxfun(@minus, xx, xx').^2 + bsxfun(@minus, yy, yy').^2); 

%% using merging criterion to detect paired neurons
flag_merge = (dist_v<=dmin(1));

% neurons should have similar decaying time constant if we want to merge
% them. 
if length(dmin)>1
    max_decay_diff = dmin(2);
else
    max_decay_diff = 2;
end
K = size(obj.C,1); 
taud = zeros(K, 1);
for m=1:K
    temp = ar2exp(obj.P.kernel_pars(m));
    taud(m) = temp(1);
end
decay_diff = abs(bsxfun(@minus, taud, taud'));
flag_merge = flag_merge & (decay_diff<max_decay_diff);

[l,c] = graph_connected_comp(sparse(flag_merge));     % extract connected components
MC = bsxfun(@eq, reshape(l, [],1), 1:c);
MC(:, sum(MC,1)==1) = [];
if isempty(MC)
    fprintf('All pairs of neurons are below the merging criterion!\n\n');
    merged_ROIs = [];
    newIDs = [];
    return;
else
    fprintf('%d neurons will be merged into %d new neurons\n\n', sum(MC(:)), size(MC,2));
end

% %% start merging
[nr, n2merge] = size(MC);
ind_del = false(nr, 1 );    % indicator of deleting corresponding neurons
merged_ROIs = cell(n2merge,1);
newIDs = zeros(nr, 1);
for m=1:n2merge
    IDs = find(MC(:, m));   % IDs of neurons within this cluster
    merged_ROIs{m} = IDs;
    
    % determine searching area
    active_pixel = (sum(A(:,IDs), 2)>0);
    
    % update spatial/temporal components of the merged neuron
    data = A(active_pixel, IDs)*C_raw(IDs, :);
    ci = C_raw(IDs(1), :);
    for miter=1:10
        ai = data*ci'/(ci*ci');
        ci = ai'*data/(ai'*ai);
    end
    % normalize ai to make its maximum to be 1
    sn = get_noise_fft(ci);
    A(active_pixel, IDs(1)) = ai*sn;
    obj.C_raw(IDs(1), :) = ci/sn;
    %     [obj.C(IDs(1), :), obj.S(IDs(1), :), tmp_kernel] = deconvCa(ci, obj.kernel, 3, true, false);
    try
        [obj.C(IDs(1), :), obj.S(IDs(1),:), deconv_options] = deconvolveCa(ci, deconv_options_0);
        obj.P.kernel_pars(IDs(1), :) = deconv_options.pars;
        newIDs(IDs(1)) = IDs(1);
        % remove merged elements
        ind_del(IDs(2:end)) = true;
    catch
        ind_del(IDs) = true;
    end
end
newIDs(ind_del) = [];
newIDs = find(newIDs);

% remove merged neurons and update obj
obj.delete(ind_del);
% obj.A(:, ind_del) = [];
% obj.C_raw(ind_del, :) = [];
% obj.C(ind_del, :) = [];
% obj.S(ind_del, :) = [];
% obj.P.kernel_pars(ind_del, :) = [];
