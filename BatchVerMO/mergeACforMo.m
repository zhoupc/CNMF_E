function  [Afinal_alldays,MC,newIDs_alldays,merged_ROIs_alldays] = mergeACforMo(Amask,ACS,merge_thr,M)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   Amask: concatnated Amask from neurons from one or many files.
%   ACS:  structure, having fields of Ain, Cin (concatnated A and C from
%         neurons from all files), and std of Cin(STD).
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%         on correlations of spatial shapes ('A'),  calcium traces ('C').
% output:
%   Afinal: merged As.
%   newIDs: cell array, dim: 1*(number of neurons after merging). Each
%       cell element has the orginal neuron number it has merged from (the
%       nueron number is cumsum across the second dim of A0s).
%   Other outputs are the same as the original quickMerge().

%%%%%%%%Older version
%(%   Merged-component reduced ACS. For those neurons that are merged, they all have the same A as well
%       as STD in each file's ACS. For example, if neuron 1,3,5 are merged,
%       neuron 1 in each file's ACS's A and STD will be the same while neuron 3,5 are
%       deleted in all files's ACS's A and STD. Since C is not used in
%       later steps, C in ACS is not updated.)
%%%%%%%%

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

C=cat(2,ACS.Cin);
STD=std(C,1,2);

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
    if isempty(ind_temp)
        continue
    else
        ind_temp=[i ind_temp];
    end
    mergegroups_intersect = cellfun(@(x) intersect(x,ind_temp),mergegroups,'UniformOutput', false);
    mergegroups_idx = find(~cellfun('isempty',mergegroups_intersect));
    if ~isempty(mergegroups_idx)
        %mergegroups{mergegroups_idx}=union(ind_temp,mergegroups{mergegroups_idx});
        mergegroups{mergegroups_idx(1)}=unique(cat(2,ind_temp,mergegroups{mergegroups_idx}));
        if length(mergegroups_idx)>1
            mergegroups(mergegroups_idx(2:end))=[];
        end
    else        
        mergegroups{end+1}=ind_temp;
    end
end
allneurons=1:size(flag_merge,1);
MC=cellfun(@(x) ismember(allneurons,x),mergegroups,'UniformOutput',false);
MC=cat(1,MC{:});
MC=MC';
numofcells=max(sum(MC,1));
merge_Bool=sum(MC,2)>0;
Aunique_Bool=~merge_Bool;

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
merged_ROIs_alldays=cell(1,numel(M));
newIDs_alldays=cell(1,numel(M));
Afinal_alldays=cell(1,numel(M));

for i=1:numel(M)
    merged_ROIs = cell(n2merge,1);
    newIDs=cell(1,nr);
    ind_del=false(nr,1);
    Afinal=zeros(size(Amask));

    A=cat(2,M{i}{:}); %This day's As.
    Afinal(:,Aunique_Bool)=A(:,Aunique_Bool);
    Aunique_ind=find(Aunique_Bool);
    for ii=1:length(Aunique_ind); newIDs{Aunique_ind(ii)} = Aunique_ind(ii); end
    for m=1:n2merge   %merge A's by their STD deviation.
        IDs = find(MC(:, m));
        merged_ROIs{m} = IDs;
        
        A_temp=A(:,MC(:,m));
        STD_temp=STD(MC(:,m));
        catSTD=diag(STD_temp./sum(STD_temp));
        weightedA=A_temp*catSTD; weightedA=reshape(weightedA,size(A,1),1,[]);
        weightedA=sum(weightedA,3);
        ind_del(IDs(2:end))= true;
        newIDs{IDs(1)} = IDs;
        Afinal(:,IDs(1))=weightedA;
    end
    Afinal(ind_del)=[];
    newIDs(ind_del) = [];
    
    merged_ROIs_alldays{i}=merged_ROIs;
    newIDs_alldays{i}=newIDs;
    Afinal_alldays{i}=Afinal;

end

end

