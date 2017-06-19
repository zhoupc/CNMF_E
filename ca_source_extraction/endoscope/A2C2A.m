function A2C2A(File, A, Acenter, Amask, i, options)
%% a greedy method for detecting ROIs and initializing CNMF. in each iteration,
% Modified from "greedyROI_endoscope"
% Main dependence is "extract_a" modified from "extract_ac" by Pengcheng.

% For each neuron, from Amask, get C, use the regression method in "extract_ac" to find A. 
% Peel A*C off Ysignal.
% Iterate for all the neurons indicated by Amask.

% Shijie Gu, techel@live.cn
%% Input:
%   Y:  d X T matrx, imaging data
%   K:  scalar, maximum number of neurons to be detected.
%   options: struct data of paramters/options
%       d1:     number of rows
%       d2:     number of columns
%       gSiz:   maximum size of a neuron
%       nb:     number of background
%       min_corr: minimum threshold of correlation for segementing neurons
%   sn:     d X 1 vector, noise level of each pixel
%   debug_on: options for showing procedure of detecting neurons
%   save_avi: save the video of initialization procedure. string: save
%   video; true: just play it; false: interactive mode. (the name of this
%      argument is very misleading after several updates of the code. sorry)

%% Output:
%`      results: struct variable with fields {'Ain', 'Cin', 'Sin', 'kernel_pars'}
%           Ain:  d X K' matrix, estimated spatial component
%           Cin:  K'X T matrix, estimated temporal component
%           Sin:  K' X T matrix, inferred spike counts within each frame
%           kernel_pars: K'X1 cell, parameters for the convolution kernel
%           of each neuron
%       center: K' X 2, coordinate of each neuron's center
%       Cn:  d1*d2, correlation image
%       save_avi:  options for saving avi.

%% parameters and preparations
d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSiz = options.gSiz;    % average size of neurons
min_pixel = options.min_pixel;  % minimum number of pixels to be a neuron
A=A;                    %concatnated A's from all files
Amask=Amask;
Acenter=Acenter;        %concatnated neuron center from all files

Ysignal=File(i).Ysignal;
Y=File(i).Y;

Ysignal(isnan(Ysignal)) = 0;    % remove nan values
Ysignal = double(Ysignal);
T = size(Ysignal, 2);

deconv_options_0= options.deconv_options;
deconv_flag = options.deconv_flag;

K = size(Amask,2);
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
Sin = zeros(K, T);      % spike counts
Cin_raw = zeros(K, T);
kernel_pars = cell(K,1);    % parameters for the convolution kernels of all neurons
STD=zeros(K,1);

%% start initialization
for mcell = 1:K;
    r=Acenter(mcell,1);
    c=Acenter(mcell,2);
    
    % Use center point data to roughly check whether this is a good starting point
    ind_p=sub2ind([d1 d2], r, c);
    y0 = Ysignal(ind_p, :);    
    y0_std = std(diff(y0));
    if max(diff(y0))< 3*y0_std % signal is weak
        Ain(:,K)=A(:,K);
        Cin_raw(K,:)=y0;
        Cin(K,:)=y0;
        continue;
    end
    
    % extract ci    
    [ci_raw,ind_success_ci] = extact_c(Ysignal,Amask);    
    Cin_raw(K,:)=ci_raw;
    
    if ~ind_success_ci
        Ain(:,K)=A(:,K);
        Cin(K,:)=ci_raw;
        continue;
    end
    
    if and(ind_success_ci,deconv_flag)
        % deconv the temporal trace
        [ci, si, deconv_options] = deconvolveCa(ci_raw, deconv_options_0, 'sn', 1);  % sn is 1 because i normalized c_raw already
        % save this initialization
        Cin(k, :) = ci;       
        Sin(k, :) = si;
        kernel_pars{k} = reshape(deconv_options.pars, 1, []);
    end
    
    ci=Cin(K,:);    
    STD(K)=std(ci);
    
    % extract ai        
    % select its neighbours for estimation of ai, the box size is
    %[2*gSiz+1, 2*gSiz+1]
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind); % size of the box
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    HY_box = Ysignal(ind_nhood, :);      % extract temporal component from HY_box
    Y_box = Y(ind_nhood, :);

    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % index of the center

    %% extract ai
    sz = [nr, nc];
    [ai, ind_success_ai] = extract_a(ci, Y_box, HY_box, Ybg_Sn_box, ind_ctr, sz);    
    if ind_success_ai==false
        Ain(:,K)=A(:,K);
    else
        Ain(:,K)=ai;
    end
    
    %% further processing ci if there is a healthy signal    
        % update the raw data
        Ysignal(ind_nhood, :) = HY_box - ai*ci;
        Y(ind_nhood, :)=Y_box - ai*ci;
   
end

File(i).Ain = sparse(Ain(:, 1:K));
File(i).Cin = Cin(1:K, :);
File(i).Cin_raw = Cin_raw(1:K, :);
File(i).STD=STD;
if deconv_flag
    File(i).Sin = Sin(1:K, :);
    File(i).kernel_pars = cell2mat(kernel_pars(1:K));
end
end