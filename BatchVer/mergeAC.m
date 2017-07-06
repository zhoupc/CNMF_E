function  [ACS,MC,newIDs,merged_ROIs] = mergeAC(Amask,ACS,merge_thr)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   Amask: concatnated Amask from neurons from one or many files.
%   ACS:  structure, having fields of Ain, Cin (concatnated A and C from
%         neurons from all files), and std of Cin(STD).
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%         on correlations of spatial shapes ('A'),  calcium traces ('C').
% output:
%   Merged-component reduced ACS. For those neurons that are merged, they all have the same A as well
%       as STD in each file's ACS. For example, if neuron 1,3,5 are merged,
%       neuron 1 in each file's ACS's A and STD will be the same while neuron 3,5 are
%       deleted in all files's ACS's A and STD. Since C is not used in
%       later steps, C in ACS is not updated.
%   Other outputs are the same as the original quickMerge().

% Author: Shijie Gu, techel@live.cn, modified from quickMerge() by Pengcheng Zhou
%  The basic idea is proposed by Eftychios A. Pnevmatikakis: high temporal
%  correlation + spatial overlap
%  reference: Pnevmatikakis et.al.(2016). Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data. Neuron
%% variables & parameters
if ~exist('merge_thr', 'var') || isempty(merge_thr) || numel(merge_thr)~=2
    merge_thr = [0.7, 0.7];
end

A_thr = merge_thr(1);
C_thr = merge_thr(2);

A=cat(2,ACS.Ain);
C=cat(2,ACS.Cin);

K = size(C,1);   % number of neurons
%% find neuron pairs to merge
% compute spatial correlation
temp = bsxfun(@times, Amask, 1./sqrt(sum(Amask)));
A_overlap = temp'*temp;

% compute temporal correlation
C_corr = corr(C')-eye(K);

% using merging criterion to detect paired neurons
flag_merge = (A_overlap>A_thr)&(C_corr>C_thr);
mergegroups={};
for i=1:size(flag_merge,1)
    ind_temp=find(flag_merge(i,:));
    mergegroups_intersect = cellfun(@(x) intersect(x,ind_temp),mergegroups,'UniformOutput', false);
    mergegroups_idx = find(~cellfun('isempty',mergegroups_intersect));
    if mergegroups_idx>0;
        mergegroups{mergegroups_idx}=union(ind_temp,mergegroups{mergegroups_idx});
    else        
        mergegroups{end+1}=ind_temp;
    end
end

allneurons=1:size(flag_merge,1);
MC=cellfun(@(x) ismember(allneurons,x),mergegroups,'UniformOutput',false);
MC=cat(1,MC{:});
MC=MC';

%%%%
% [l,c] = graph_connected_comp(sparse(flag_merge));     % extract connected components
% MC = bsxfun(@eq, reshape(l, [],1), 1:c);
% MC(:, sum(MC,1)==1) = [];
%%%%
if isempty(MC)
    fprintf('All pairs of neurons are below the merging criterion!\n\n');
    newIDs=[]; merged_ROIs=[];
    return;
else
    fprintf('%d neurons will be merged into %d new neurons\n\n', sum(MC(:)), size(MC,2));
end

[nr, n2merge] = size(MC);
merged_ROIs = cell(n2merge,1);
newIDs = zeros(nr, 1);
ind_del=false(nr,1);

% start merging
for m=1:n2merge
    %oldIDs=IDs;
    IDs = find(MC(:, m));  % IDs of neurons within this cluster
%    IDs=setdiff(IDs,oldIDs);
    merged_ROIs{m} = IDs;
    
    % determine searching area
    active_pixel = sum(A(:,IDs), 2)>0;

    
    % update spatial/temporal components of the merged neuron   
    % data = A(active_pixel, IDs)*C(IDs, :);
    data=[];
    for i=1:numel(ACS)
        FileA=File(i).Ain;
        FileC=File(i).Cin;
        data = [data FileA(active_pixel, IDs)*FileC(IDs, :)];
    end
    
    data=data./length(IDs);
    [~,I] = max(std(C(IDs, :))); % choose the most confident(with biggest std) ci.
    ci=C(IDs(I),:);
    for miter=1:10
        ai = data*ci'/(ci*ci');
        ci = ai'*data/(ai'*ai);
    end
    ind_del(IDs(2:end))=true;
    newIDs(IDs(1)) = IDs(1);
    % making ai nicer.
%     temp = ai>quantile(ai, 0.3, 1);
%     ai(~temp(:)) = 0;
   
    %Amask(:,IDs(1)) = ai>0;
    %C(IDs(1), :) = ci;
    for i=1:numel(ACS)
        ACS(i).Ain(:,IDs(1))=ai;
        ACS(i).STD(IDs(1))=std(ci);
        %FileSTD=ACS(i).STD; FileSTD(IDs(1))=std(ci);   ACS(i).STD=FileSTD;
    end    
end
% 
% Amask(:,ind_del)=[];
% C(ind_del, :) = [];
for i=1:numel(ACS)
    ACS(i).Ain(:,ind_del)=[];
    ACS(i).STD(ind_del)=[];
%     FileA=ACS(i).Ain;   FileA(:,ind_del)=[]; ACS(i).Ain=FileA;
%     FileSTD=ACS(i).STD; FileSTD(ind_del)=[]; ACS(i).STD=FileSTD;
end   
newIDs(ind_del) = [];
newIDs = find(newIDs);
end
