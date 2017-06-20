function [ci,ind_success] = extract_c(Ysignal,Amask)
% Use neuron k's Amask to find in Ysignal all pixels involved in this neuron.
% Use the mean temperal trace of these pixels extact c
cs=Ysignal(Amask,:);
ci=mean(cs);

if norm(ci)==0  % avoid empty results
    ind_success=false;
else
    ind_success=true;
end

end