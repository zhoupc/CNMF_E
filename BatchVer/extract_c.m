function [ci,ind_success] = extract_c(Ysignal,Amask,A)
% Use neuron k's Amask to find in Ysignal all pixels involved in this neuron.
% Use the mean temperal trace of these pixels extact c
if ~isempty(Amask)    
    cs=Ysignal(Amask,:);
    ci=mean(cs);
    
    if norm(ci)==0  % avoid empty results
        ind_success=false;
    else
        ind_success=true;
    end
elseif ~isempty(A)
    X = A; 
    ci = (X'*X)\(X'*Ysignal);
    ind_success=[];
else
    display('error in extract_c')
end
end