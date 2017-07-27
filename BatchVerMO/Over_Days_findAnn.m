function [ns_storage]=Over_Days_findAnn(As,correlation_thresh,max2max2nd,skewnessthresh)
% Input: 
%   As={A1,A2,A3,A4...}
%   Filters for saying they are the same neuron: correlation_thresh, ratio of first max/second max correlation coef, skewnessthresh.
%   Type help findn1n2 for more information for skewnessthresh.
% Output:
%   ns_storage: each column for each day, index for k in A or C such that
%       when you specify A(ns_storage(:,1)) and A(ns_storage(:,2)),
%       each row of the permutated A?s from time point 1 and point 2 are the index for the same neuron in the original As

% Shijie Gu, techel@live.cn (ShanghaiTech University, Harvard-MIT)
if or(isempty(As),length(As)==1)
    error('As is empty or has only one A in it. No sense to use this function. Please check data selected.')
elseif length(As)==2;
    ns_storage=findn1n2(As{1},As{2},correlation_thresh,max2max2nd,skewnessthresh);
else
    ns_storage=findn1n2(As{1},As{2},correlation_thresh,max2max2nd,skewnessthresh);
    for i=2:length(As)-1
        ns_next=findn1n2(As{i},As{i+1},correlation_thresh,max2max2nd,skewnessthresh);
        [ind_row_storage,ind_row_next] = ismember(ns_storage(:,end),ns_next(:,1));
        ns_storage=ns_storage(ind_row_storage,:);
            ind_row_next=ind_row_next(ind_row_next~=0); 
            ns_next=ns_next(:,2);
        ns_next=ns_next(ind_row_next);
        ns_storage=[ns_storage ns_next];
    end
end
end