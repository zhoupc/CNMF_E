function  [merged_ROIs,Ait,Cit,ind_del] = mergeAC(Amask,File,merge_thr)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   Amask: concatnated Amask from neurons from all files.
%   Concatnated A and C from neurons from all files.
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%   on correlations of spatial shapes ('A'),  calcium traces ('C').

% output:
%   merged_ROIs: cell arrarys, each element contains indices of merged components.
%   newIDs: vector, each element is the new index of the merged neurons
%   merged A and C

%% Author: Shijie Gu, techel@live.cn, modified from quickMerge() by Pengcheng Zhou
%  The basic idea is proposed by Eftychios A. Pnevmatikakis: high temporal
%  correlation + spatial overlap
%  reference: Pnevmatikakis et.al.(2016). Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data. Neuron

%% variables & parameters
if ~exist('merge_thr', 'var') || isempty(merge_thr) || numel(merge_thr)~=2
    merge_thr = [0.7, 0.7];
end

A_thr = merge_thr(1);
C_thr = merge_thr(2);
C=cat(2,File.Cin);
A=cat(2,File.Ain);
K = size(C,1);   % number of neurons
%% find neuron pairs to merge
% compute spatial correlation
temp = bsxfun(@times, Amask, 1./sqrt(sum(Amask)));
A_overlap = temp'*temp;

% compute temporal correlation
C_corr = corr(C')-eye(K);

% using merging criterion to detect paired neurons
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

[nr, n2merge] = size(MC);
ind_del = false(nr, 1);    % indicator of deleting corresponding neurons
merged_ROIs = cell(n2merge,1);
newIDs = zeros(nr, 1);
Ait=zeros(size(Amask,1),n2merge);
Cit=zeros(n2merge,size(C,2));

% start merging
for m=1:n2merge
    IDs = find(MC(:, m));   % IDs of neurons within this cluster
    merged_ROIs{m} = IDs;
    
    % determine searching area
    active_pixel = (sum(A(:,IDs), 2)>0);
    
    % update spatial/temporal components of the merged neuron
    
    %data = A(active_pixel, IDs)*C(IDs, :);
    data=[];
    for i=1:numel(File)
        FileA=File(i).Ain;
        FileC=File(i).Cin;
        data = [data FileA(active_pixel, IDs)*FileC(IDs, :)];
    end
    
    ci = max(C(IDs, :),[],1);
    for miter=1:10
        ai = data*ci'/(ci*ci');
        ci = ai'*data/(ai'*ai);
    end
    
    Ait(active_pixel,m) = ai;
    Cit(m, :) = ci;    
    ind_del(IDs) = true;
end

% % remove merged neurons and update
% 
% A(:, ind_del) = [];
% C(ind_del, :) = [];

end
