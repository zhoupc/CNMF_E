function As=Over_Days_ResequenceA(As,correlation_thresh,max2max2nd,skewnessthresh)
% General description: This function puts neurons that are roughly checked(findn1n2)
% as similar over two consecutive days in similar order. This should help
% make cnmfe-BatchVer more robust since the sequence of initiation might
% effect the result of extracted temporal traces. This function intends to
% preserve the sequence of some similar neurons. The idea is illustrated in
% the following example. If in day1 neuron 2 and 3 is tracked similared to
% neuron 4 and 3 in day2. Then the neurons in day2 will be resequenced as 3
% and 4 while other neurons that have nothing to do with day1's neuron do
% not move.
% Input: 
%   As={A1,A2,A3,A4...}
%   Filters for saying they are the same neuron: correlation_thresh, ratio of first max/second max correlation coef, skewnessthresh.
%   Type help findn1n2 for more information for skewnessthresh.
% Output:
%   As: re-sequenced As.
% Shijie Gu, techel@live.cn, ShanghaiTech University; Fee Lab at MIT-BCS.

if or(isempty(As),length(As)==1)
    error('As is empty or has only one A in it. No sense to use this function. Please check data selected.')
else
    for i=1:length(As)-1
        ns_storage=findn1n2(As{i},As{1+1},correlation_thresh,max2max2nd,skewnessthresh);
        if or(isempty(ns_storage),size(ns_storage,1)==1) % no need to re-sequence them.
            continue
        else
            [~,ns_storage_descend_i]=sort(ns_storage(:,1),1,'descend');
            ns_storage=ns_storage(ns_storage_descend_i,:);
            new_location=ns_storage(:,2);
            all_location=1:size(As{i+1},2);
            new_location_old=sort(new_location);
            all_location(new_location_old)=new_location;
            As{i+1}=As{i+1}(:,all_location);
        end
    end 
end