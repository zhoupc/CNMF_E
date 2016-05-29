function [ai, ci, ind_success] = extract_ac(HY, Y, ind_ctr, sz) 
%% given a patch of raw & high-pass filtered calcium imaging data, extract 
% spatial and temporal component of one neuron (ai, ci). if succeed, then
% return an indicator ind_succes with value 1; otherwise, 0. 
%% input: 
%       HY:     d X T matrix, filtered patch data 
%       Y:      d X T matrix, raw data 
%       ind_ctr:        scalar, location of the center 
%       sz:         2 X 1 vector, size of the patch 

nr = sz(1); 
nc = sz(2); 
min_corr = 0.85; 

%% find pixels highly correlated with the center 
HY(HY<0) = 0;       % remove some negative signals from nearby neurons 
y0 = HY(ind_ctr, :); 
tmp_corr = reshape(corr(y0', HY'), nr, nc);
data = HY(tmp_corr>min_corr, :); 


%% estimate ci with the mean or rank-1 NMF
ci = mean(data, 1); 
% ci = ci - min(ci); % avoid nonnegative baseline 
% [~, ci] = nnmf(ci, 1); 
if norm(ci)==0
    ai=[]; 
    ind_success=false; 
    return; 
end

%% extract spatial component 
% estiamte the background level using the boundary 
indr = [ones(1, nc), ones(1, nc)*nr, 1:nr, 1:nr];
indc = [1:nc, 1:nc, ones(1, nr)*nc, ones(1,nr)];
ind_bd = sub2ind([nr, nc], indr, indc);     % indices for boundary pixels 
y_bg = median(Y(ind_bd, :), 1);  % take the mean in the boundary as an estimation of the background level 

% sort the data, take the differencing and estiamte ai 
[~, ind_sort] = sort(y_bg, 'ascend'); 
dY = diff(Y(:, ind_sort), 1, 2); 
dY(bsxfun(@lt, abs(dY), std(dY, 0, 2)*2)) = 0; 
dci = diff(ci(ind_sort)); 
dci(dci<std(dci, 0, 2)*2) = 0; 
ai = max(0, dY*dci'/(dci*dci')); 

% post-process ai by bwlabel
temp = full(ai>max(ai)/10);
l = bwlabel(reshape(temp, nr, nc), 4);   % remove disconnected components
temp(l~=l(ind_ctr)) = false;
ai(~temp(:)) = 0;

% estimate ci again 
ind_nonzero = (ai>0); 
ai_mask = mean(ai(ind_nonzero))*ind_nonzero; 
ci = (ai-ai_mask)'*ai\((ai-ai_mask)'*Y); 
ci = ci - median(ci); 
ci(ci<-std(diff(ci))) = 0; 
% return results 
if norm(ai)==0
    ind_success= false; 
else
    ind_success=true; 
end