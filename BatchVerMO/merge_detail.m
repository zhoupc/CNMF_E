function MC=merge_detail(flag_merge)
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