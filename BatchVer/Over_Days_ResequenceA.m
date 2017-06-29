function [ns_storage,Ass]=Over_Days_ResequenceA(As,correlation_thresh,max2max2nd,skewnessthresh)
% Input: 
%   As={A1,A2,A3,A4...}
%   Filters for saying they are the same neuron: correlation_thresh, ratio of first max/second max correlation coef, skewnessthresh.
%   Type help findn1n2 for more information for skewnessthresh.
% Output:
%   ns_storage: each column for each day, index for k in A or C such that
%       when you specify A(ns_storage(:,1)) and A(ns_storage(:,2)),
%       each row of the permutated As from time point 1 and point 2 are the index for the same neuron in the original As
%   Ass: in some cases it is impossible to find the same neurons that can be 
%       tracked over many many days, i.e, the ns_storage will be empty. It
%       would be nice to have some resequencing technique that only depends
%       on adjacent two days' A. The implementation here resequence all the
%       A's across days based on the new incoming day's A. Similar neurons
%       will be put in front with those that could not found similar ones
%       tagged behind. Ass is the resequenced As.
% Shijie Gu, techel@live.cn

if or(isempty(As),length(As)==1)
    error('As is empty or has only one A in it. No sense to use this function. Please check data selected.')
elseif length(As)==2;
    ns_storage=findn1n2(As{1},As{2},correlation_thresh,max2max2nd,skewnessthresh);
    ns1=ns_storage(:,1); nn1=1:size(As{1},2); nu1=setdiff(nn1,ns1);
    ns2=ns_storage(:,2); nn2=1:size(As{2},2); nu2=setdiff(nn2,ns2);
    As{1}=[As{1}(:,ns1) As{1}(:,nu1)];
    As{2}=[As{2}(:,ns2) As{2}(:,nu2)];
else
    ns_storage=findn1n2(As{1},As{2},correlation_thresh,max2max2nd,skewnessthresh);
    ns1=ns_storage(:,1); nn1=1:size(As{1},2); nu1=setdiff(nn1,ns1);
    ns2=ns_storage(:,2); nn2=1:size(As{2},2); nu2=setdiff(nn2,ns2);
    As{1}=[As{1}(:,ns1) As{1}(:,nu1)];
    As{2}=[As{2}(:,ns2) As{2}(:,nu2)];
    
    for i=2:length(As)-1
        ns_next=findn1n2(As{i},As{i+1},correlation_thresh,max2max2nd,skewnessthresh);
        [ind_row_storage,ind_row_next] = ismember(ns_storage(:,end),ns_next(:,1));
        
            ind_row_next=ind_row_next(ind_row_next~=0); 
            ns_next=ns_next(:,2);
        ns_next=ns_next(ind_row_next);
        nn1=1:size(ns_storage,1); ns1=nn1(ind_row_storage);
                
        for j=1:i
            nn1=1:size(As{j},2);
            nu1=setdiff(nn1,ns1);
            As{j}=[As{j}(:,ns1) As{j}(:,nu1)];
        end
        ns2=ns_next; nn2=1:size(As{i+1},2); nu2=setdiff(nn2,ns2);
        As{i+1}=[As{i+1}(:,ns2) As{i+1}(:,nu2)];
        ns_storage=ns_storage(ind_row_storage,:);
        ns_storage=[ns_storage ns_next];
    end
end
Ass=As;
end