function [ai, ind_success] = extract_a(ci, Y_box, HY_box, Ybg_Sn_box, ind_ctr, sz)
% given the temporal component, background-subtracted and denoised calcium imaging data, extract
% spatial component of one neuron ai. Method is regression. Y=ac+Ybg+noise. Afterwards, 
% we threshold the spatial shape and remove those too small or emoty results. those If succeed, then
% return an indicator ind_succes with value 1; otherwise, 0.
% some inputs explained:
%       ind_ctr:    scalar, location of the center
%       sz:         2 X 1 vector, size of the patch

% Shijie Gu, techel@live.cn. modified from extract_ac() by pcZ.

%% preparations 
nr = sz(1);
nc = sz(2);
min_pixels = 5;
y_bg = mean(Ybg_Sn_box);

%% estimate ai 
T = length(ci); 
X = [ones(T,1), y_bg', ci']; 
temp = (X'*X)\(X'*Y'); 
ai = max(0, temp(3,:)'); 

%% threshold the spatial shape and remove those too small or emoty results.
% remove outliers not continuous with the main.
temp =  full(ai>quantile(ai(:), 0.5)); 
l = bwlabel(reshape(temp, nr, nc), 4); 
temp(l~=l(ind_ctr)) = false; 
ai(~temp(:)) = 0; 

if sum(ai(:)>0) < min_pixels %the ROI is too small
    ind_success=false;
    return;
end

if or(norm(ai)==0,any(isnan(ai)))
    ind_success= false;
else
    ind_success=true;
end