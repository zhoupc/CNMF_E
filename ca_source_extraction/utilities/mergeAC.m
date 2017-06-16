function  [merged_ROIs, newIDs,A,C_raw] = mergeAC(A,C_raw,merge_thr)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%   on correlations of spatial shapes ('A'),  calcium traces ('C').

% output:
%   merged_ROIs: cell arrarys, each element contains indices of merged
%   components
%   newIDs: vector, each element is the new index of the merged neurons

%% Author: Shijie Gu, ShanghaiTech University.
%  The basic idea is proposed by Eftychios A. Pnevmatikakis: high temporal
%  correlation + spatial overlap
%  reference: Pnevmatikakis et.al.(2016). Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data. Neuron

%% variables & parameters

if ~exist('merge_thr', 'var') || isempty(merge_thr) || numel(merge_thr)~=3
    merge_thr = [0.5, 0.7, 0];
end

A_thr = merge_thr(1);
C_thr = merge_thr(2);
K = size(C_raw,1);   % number of neurons


%% find neuron pairs to merge
% compute spatial correlation
% temp = bsxfun(@times, A, 1./sum(A.^2,1));
temp = bsxfun(@times, A>0, 1./sqrt(sum(A>0)));
A_overlap = temp'*temp;

% compute temporal correlation
C_corr = corr(C_raw')-eye(K);


%% using merging criterion to detect paired neurons
flag_merge = (A_overlap>A_thr)&(C_corr>C_thr);

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
    C_raw(IDs(1), :) = ci/sn;
    
    ind_del(IDs(2:end)) = true;

end

newIDs(ind_del) = [];
newIDs = find(newIDs);

% remove merged neurons and update

A(:, ind_del) = [];
C_raw(ind_del, :) = [];

end
