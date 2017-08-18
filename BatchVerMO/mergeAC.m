function  [Afinal,MC,newIDs,merged_ROIs,close_ind,real_ind] = mergeAC(A,ACS,merge_thr,dmin,d1,d2)
%% merge neurons based on simple spatial and temporal correlation
% input:
%   Amask: concatnated Amask from neurons from one or many files.
%   ACS:  structure, having fields of Ain, Cin (concatnated A and C from
%         neurons from all files), and std of Cin(STD).
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%         on correlations of spatial shapes ('A'),  calcium traces ('C').
%   dmin: min distance for two neurons to be called different.
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

%A=cat(2,ACS.Ain);
%C=cat(2,ACS.Cin_raw);
C=cat(2,ACS.Cin);
STD=std(C,1,2);
qt=quantile(STD,10); %pick the last 10 percent to waste.
real_ind=STD>qt(1);  %%
A=A(:,real_ind);
C=C(real_ind,:);
STD=STD(real_ind);

Amask=A>0;
K = size(C,1);   % number of neurons
%% find neuron pairs to merge
% compute spatial correlation
temp = bsxfun(@times, Amask, 1./sqrt(sum(Amask)));
clear Amask
A_overlap = temp'*temp;

% compute temporal correlation
C_corr = corr(C')-eye(K);

% compute center distance
ctr = round(com(A, d1, d2));
yy = ctr(:,1);
xx = ctr(:,2);
dist_v = sqrt(bsxfun(@minus, xx, xx').^2 + bsxfun(@minus, yy, yy').^2); 


% using merging criterion to detect paired neurons
flag_merge1 = (A_overlap>A_thr)&(C_corr>C_thr);
flag_merge2 = (dist_v<=dmin);
flag_merge=or(flag_merge1,flag_merge2);

MC=merge_detail(flag_merge);
MC1=merge_detail(flag_merge1);
[~,close_ind]=setdiff(MC',MC1','rows','stable');

%%%%
% [l,c] = graph_connected_comp(sparse(flag_merge));     % extract connected components
% MC = bsxfun(@eq, reshape(l, [],1), 1:c);
% MC(:, sum(MC,1)==1) = [];
%%%%

if isempty(MC)
    fprintf('All pairs of neurons are below the merging criterion.\n\n');
    newIDs=[]; merged_ROIs=[];
    return;
else
    fprintf('%d neurons will be merged into %d new neurons\n\n', sum(MC(:)), size(MC,2));
end

[nr, n2merge] = size(MC);
merged_ROIs = cell(n2merge,1);
newIDs=cell(1,nr); %newIDs = num2cell(1:nr);
ind_del=true(nr,1);
Afinal=zeros(size(A));


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
%%%%%%%%%    
%     data=[];
%     for i=1:numel(ACS)
%         FileC=ACS(i).Cin;
%         data = [data A(active_pixel, IDs)*FileC(IDs, :)];
%     end
%     
%     data=data./length(IDs);
%     [~,I] = max(std(C(IDs, :),0,2)); % choose the most confident(with biggest std) ci.
%     ci=C(IDs(I),:);
%     for miter=1:10
%         ai = data*ci'/(ci*ci');
%         ci = ai'*data/(ai'*ai);
%     end
%%%%%%%%
    A_temp=A(active_pixel,MC(:,m));
    STD_temp=STD(MC(:,m));
    catSTD=STD_temp./sum(STD_temp);

    weightedA=A_temp*catSTD;
    
    ind_del(IDs(1))= false;
    newIDs{IDs(1)} = IDs;
    Afinal(active_pixel,IDs(1))=weightedA;
    % making ai nicer.
%     temp = ai>quantile(ai, 0.3, 1);
%     ai(~temp(:)) = 0;
   
    %Amask(:,IDs(1)) = ai>0;
    %C(IDs(1), :) = ci;
%     for i=1:numel(ACS)
%         ACS(i).Ain(active_pixel,IDs(1))=ai;
%         ACS(i).STD(IDs(1))=std(ci);
%         %FileSTD=ACS(i).STD; FileSTD(IDs(1))=std(ci);   ACS(i).STD=FileSTD;
%     end    
end

% for i=1:numel(ACS)
%     ACS(i).Ain(:,ind_del)=[];
%     ACS(i).STD(ind_del)=[];
% %     FileA=ACS(i).Ain;   FileA(:,ind_del)=[]; ACS(i).Ain=FileA;
% %     FileSTD=ACS(i).STD; FileSTD(ind_del)=[]; ACS(i).STD=FileSTD;
% end

newIDs(ind_del) = [];
Afinal(:,ind_del) = [];

% newIDs(ind_del) = [];
% newIDs = find(newIDs);
end
