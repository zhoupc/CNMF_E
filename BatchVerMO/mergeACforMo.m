function  [Afinal_alldays,MC,newIDs,merged_ROIs,close_ind] = mergeACforMo(A,C,merge_thr,M,dmin,d1,d2)
% General Description: Merge neurons based on
%   spatial, temporal correlation, and their distance. Afinal will be
%   different for each day so that they are put into cell array where each day,
%   has its own cell. Yet, each day's A is merged from the exact neurons using the
%   exact weight matrix: max(diff(C,1,2))./STD, see 'C' in input. This
%   function is key in tracking neuron.
%   This function is very similar to mergeAC, but this one does not
%   eliminate neurons with low signal. It also requires C as opposed to
%   ACS. Its emphasis is on the handling of multi-day data in the same way.
% input:
%   A:    concatenated A from neurons from many files registered towards
%           one day. That day should be the most confidently registered A.
%   C:    BigC(in ReadMe) of each file concatenated in time.
%   merge_thr: 1X2 vector, threshold for two metrics {'A', 'C'}. it merge neurons based
%         on correlations of spatial shapes ('A'),  calcium traces ('C').
%   M:    cell array of length number of days. Within each day's cell, it
%           is As registered towards that day. More detailed in ReadMe for MoBatchVer.
%   dmin: min distance for two neurons to be called different.
%   d1 and d2: row and column of the FOV, for calculating center of each
%         neuron.
% output:
%   Afinal_alldays: merged As for each day. Each day's cell array is put
%       into a seperate cell. For the following output variables, it is essentially
%       the same as MergeAC. For ForMo version, each result becomes a cell
%       array, with ex
%   MC:  see 'help merge_detail'.
%   newIDs: cell array, dim: 1*(number of neurons after merging). Each
%       cell element has the orginal neuron number it has merged from (the
%       nueron number is cumsum across the second dim of A0s).
%   merged_ROIs: essentially the same as newIDs. Left over from previous
%       version. Keep if for now.
%   close_ind: ind of neurons for output (that are merged) just because some neurons are
%       close together.

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

STD=std(C,1,2);
STD=max(diff(C,1,2))./STD;

K = size(C,1);   % number of neurons
%% find neuron pairs to merge
% compute spatial correlation
Amask=A>0;
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

clear A_overlap C_corr C flag_merge1 flag_merge2 flag_merge;
display('Deleted some big variables.')


merge_Bool=sum(MC,2)>0;
Aunique_Bool=~merge_Bool; %keep it here for now, future versions may need this.

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
%merged_ROIs_alldays=cell(1,numel(M));
%newIDs_alldays=cell(1,numel(M));
Afinal_alldays=cell(1,numel(M));

for i=1:numel(M)
    Afinal=zeros(size(A));
    ind_del=true(nr,1);
    if i==1
        merged_ROIs = cell(n2merge,1);
        newIDs=cell(1,nr);        
    end

    A=cat(2,M{i}{:}); %This day's As.

    for m=1:n2merge   %merge A's by their STD deviation. 
        display(m)
        IDs = find(MC(:, m));
        ind_del(IDs(1))= false;
        
        A_temp=A(:,MC(:,m));
        STD_temp=STD(MC(:,m));
        catSTD=STD_temp./sum(STD_temp);

        weightedA=A_temp*catSTD;
        
        Afinal(:,IDs(1))=weightedA;
        
        if i==1            
            merged_ROIs{m} = IDs;            
            newIDs{IDs(1)} = IDs;
        end
    end  
    Afinal(:,ind_del)=[];
    Afinal_alldays{i}=Afinal;
    if i==1
        newIDs(ind_del) = [];
    end
end

