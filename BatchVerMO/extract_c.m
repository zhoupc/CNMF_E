function [ci,ind_success] = extract_c(Ysignal,Amask,A)
% [ci,ind_success] = extract_c(Ysignal,Amask,A)
% General description and Input:
% This function uses one of the two ways to extract temporal traces of a neuron.
% They all need one denoised/filtered/background-subtracted Ysiganl.
%   When you only provide Amask of one neuron into the function, it uses this neuron's Amask to find in Ysignal all pixels involved in this neuron.
%       The mean temperal trace of these pixels will become c.
%   When you provide A, (it will ignore Amask), it will extract one neuron or many neurons at once
%       depending on if you provide Ysignal of one neuron or many neurons. In this condition, it uses
%       regression to get ci.
% Output:
    % (1)Note that this function is on one-neuron scale. Use for loop outside to go
    % through all of the neurons when you use it.
    % (2)No matter true or false is the "ind_success", ci will always have
    % something. "ind_success" is an indicator for further analysis decision 
    % such as, deconcolve or not, etc.
% This function is mainly used in A2C2A, but is also used by itself to
    % extract c such as in the "massive" mode in "demo_endoscope2". 
    
%%Shijie Gu, techel@live.cn (ShanghaiTech University; Fee Lab at MIT-BCS)

if ~isempty(A)
%     X = A; 
%     ci = (X'*X)\(X'*Ysignal);
%     ci(isnan(ci)) = 0;
    ind_nonzero = A>0;
    ai_mask = mean(A(ind_nonzero))*ind_nonzero;
    ci = (A-ai_mask)'*A\((A-ai_mask)'*Ysignal);
    ci(isnan(ci)) = 0;
elseif ~isempty(Amask)    
    cs=Ysignal(Amask,:);
    ci=mean(cs);
    ci(isnan(ci)) = 0; 
else
    error('error in extract_c, please provide A or Amask.')
end
    
if norm(ci)==0  % avoid empty results
    ind_success=false;
else
    ind_success=true;
end

end