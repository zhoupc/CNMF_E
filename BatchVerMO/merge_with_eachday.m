function MC_new=merge_with_eachday(M,MC)
% This function filters MC so that those MC columns that have neurons merged from
%   each day will stay
% M, the M in mergeACforMo
% MC, the MC you want to filter

eachday_nn=cellfun(@(x) size(x,2),M{1});
MC_new=[];
for m=1:size(MC,2)
    MC_current=MC(:,m);
    MC_current_split=mat2cell(MC_current,eachday_nn,1);
    existence_each_day=cellfun(@(x) any(any(x)),MC_current_split);
    if sum(existence_each_day)==numel(M)
        MC_new(:,end+1)=MC(:,m);
    end
end