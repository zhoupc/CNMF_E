function threshold_value=threshold_by_diff(P,resolution,number)
% threshold is where consecutive diff smaller than 'resolution' for 'number' pixels
% threshold_value returned in the end is the actual value of P at the
%   threshold.

% Shijie Gu, Fee Lab at McGovern, MIT-BCS; ShanghaiTech University.

if isempty(number)
    number=10;
end

PNR0_sorted=sort(P,'descend');
PNR0_diff=[2; diff(PNR0_sorted)]'; %2 is a random number
[L_low, numRegions_low] = bwlabel(abs(PNR0_diff)<resolution);
sizes_low=[];for i=1:numRegions_low; sizes_low=[sizes_low,sum(L_low==i)]; end
[L_high, numRegions_high] = bwlabel(abs(PNR0_diff)>resolution);
sizes_high=[];for i=1:numRegions_high; sizes_high=[sizes_high,sum(L_high==i)]; end
if find(sizes_low>number)~=1
    threshold=sum(sizes_high(1:find(sizes_low>number)-1))+sum(sizes_low(1:find(sizes_low>number)-1));
    threshold_value=PNR0_sorted(threshold);
else
    threshold_value=PNR0_sorted(1);
    warning('Please lower the ''resolution'' or increase ''number''.')
end

