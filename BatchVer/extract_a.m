function [ai, ai_raw, ind_success] = extract_a(ci, Y_box, HY_box, Amask, ind_ctr, sz, sn, options)
% given the temporal component, background-subtracted and denoised calcium imaging data, extract
% spatial component of one neuron ai. Method is regression. 
    % Y=ac+Ybg+noise for non-denoised Y. 
    % Y=ac for denoised Y.
    % Use sn==1 or [] to tell function which one suits current supply and demand.
    % Input:
    % for both, ind_ctr, sz, are needed for quality control.
    %(1) Y=ac+Ybg+noise for non-denoised Y will need everything.
    %(2) Y=ac for denoised Y. can have Amask=[] and Y=[].
% Afterwards, 
% we threshold the spatial shape and remove those too small or empty results. For those succeeded, then
% return an indicator ind_succes with value 1; otherwise, 0.
% some inputs explained:
%       ind_ctr:    scalar, location of the center
%       sz:         2 X 1 vector, size of the patch

%%Shijie Gu, techel@live.cn. modified from extract_ac() by Pengcheng Zhou.

%% preparations 
if isempty(ci)
    error('No ci is provided.')
end
nr = sz(1);
nc = sz(2);
if isempty(options)
    min_pixels=15;
else
    min_pixels = options.min_pixels;
end

HY=HY_box;

%% estimate ai    
if sn==1
    X=ci';
    temp = (X'*X)\(X'*HY');
    ai = max(0, temp');
else
    Y=Y_box;
    T = length(ci); 
    y_bg = median(HY(~Amask, :), 1); 
    X = [ones(T,1), y_bg', ci']; 
    temp = (X'*X)\(X'*Y'); 
    ai = max(0, temp(3,:)'); 
end

%%
ai(isnan(ai))=0;
ai_raw=ai;
%% threshold the spatial shape and remove those too small or empty results.
% remove outliers not continuous with the main.
temp =  full(ai>quantile(ai(:), 0.5)); 
l = bwlabel(reshape(temp, nr, nc), 4); 
temp(l~=l(ind_ctr)) = false; 
ai(~temp(:)) = 0; 

if sum(ai(:)>0) < min_pixels %the ROI is too small
    ind_success=false;
elseif norm(ai)==0
    ind_success= false;
else
    ind_success=true;
end