function  [Amask,C,ACS] = mergeAC(Amask,C,ACS,merge_thr)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   Amask: concatnated Amask from neurons from one or many files.
%   C: concatnated temporal trace from neurons from one or many files.
%   ACS:  structure, having fields of Ain, Cin (concatnated A and C from
%         neurons from all files), and STD of Cin.
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%         on correlations of spatial shapes ('A'),  calcium traces ('C').
% output:
%   Merged-component reduced Amask, C and ACS. After merging, previously
%       each file's ACS has different A. Now they all have the same A as well
%       as STD.

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
    return;
else
    fprintf('%d neurons will be merged into %d new neurons\n\n', sum(MC(:)), size(MC,2));
end

[nr, n2merge] = size(MC);
merged_ROIs = cell(n2merge,1);
ind_del=false(nr,1);

% start merging
for m=1:n2merge
    oldIDs=IDs;
    IDs = find(MC(:, m));  % IDs of neurons within this cluster
%    IDs=setdiff(IDs,oldIDs);
    merged_ROIs{m} = IDs;
    
    % update spatial/temporal components of the merged neuron   
    % data = A(active_pixel, IDs)*C(IDs, :);
    data=[];
    C_start=1;
    for i=1:numel(ACS)
        data = [data ACS(i).Ain(:, IDs)*C(IDs, C_start:(C_start+size(ACS(i).Cin,2)-1))];
        C_start=C_start+size(ACS(i).Cin,2);
    end
    
    data=data./length(IDs);
    [~,I] = max(std(C(IDs, :),0,2)); % choose the most confident(with biggest std) ci.
    ci=C(IDs(I),:);
    for miter=1:10
        ai = data*ci'/(ci*ci');
        ci = ai'*data/(ai'*ai);
    end
    ind_del(IDs(2:end))=true;
    % making ai nicer.
%     temp = ai>quantile(ai, 0.3, 1);
%     ai(~temp(:)) = 0;
   
    Amask(:,IDs(1)) = ai>0;    
    C(IDs(1), :) = ci;
    for i=1:numel(ACS)
        FileA=ACS(i).Ain; FileA(:,IDs(1))=ai;   
        FileSTD=ACS(i).STD; FileSTD(IDs(1))=std(ci);
    end    
end
Amask(:,ind_del)=[];
C(ind_del, :) = [];
for i=1:numel(ACS)
    FileA(:,ind_del)=[]; ACS(i).Ain=FileA;
    FileSTD(ind_del)=[]; ACS(i).STD=FileSTD;
end   

end
