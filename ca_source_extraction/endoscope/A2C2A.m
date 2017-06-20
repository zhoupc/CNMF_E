function ACS=A2C2A(File, A0s, j, options, sn)
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
%   sn: For background and denoised
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

if ~isempty(j)
    A=A0s{j};                    %concatnated A's from all files
elseif isempty(j)
    A=A0s;
end
    
Acenter = com(A, d1, d2); %nr x 2 matrix, with the center of mass coordinates
Amask=A>0;

if isempty(File)
    File=[];
else
    File=File;
end
   
if isempty(sn) % means it is not denoised, used in normal cnmf-e
    global Y
    % create a spatial filter for removing background
    gSig = options.gSig;    % width of the gaussian kernel approximating one neuron
    psf = fspecial('gaussian', round(gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;
    
    % filter the data
    HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate');
    HY = reshape(HY, d1*d2, []);
    HY = bsxfun(@minus, HY, median(HY, 2));
    Ysignal=HY;
    
elseif sn==1 
    Ysignal=File.Ysignal;
    Y=File.Y;
else
    display('This function only takes input of sn being 1 or []')
    return
end

Ysignal(isnan(Ysignal)) = 0;    % remove nan values
Ysignal = double(Ysignal);
T = size(Ysignal, 2);

deconv_options_0= options.deconv_options;
deconv_flag = options.deconv_flag;

K = size(Amask,2);
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
%Sin = zeros(K, T);      % spike counts
%Cin_raw = zeros(K, T);
%kernel_pars = cell(K,1);    % parameters for the convolution kernels of all neurons
STD=zeros(K,1);
del_ind=false(1,K);
%% start initialization
for k = 1:K;
    r=Acenter(k,1);
    c=Acenter(k,2);
    
    % Use center point data to roughly check whether this is a good starting point
    ind_p=sub2ind([d1 d2], r, c);
    y0 = Ysignal(ind_p, :);    
    y0_std = std(diff(y0));
    if max(diff(y0))< 3*y0_std % signal is weak
        Ain(:,k)=A(:,k);
%        Cin_raw(k,:)=y0;
        Cin(k,:)=y0;
        continue;
    end
    
    % extract ci    
    [ci_raw,ind_success_ci] = extract_c(Ysignal,Amask(:,k));    
%    Cin_raw(k,:)=ci_raw;
    
    if ~ind_success_ci
        Ain(:,k)=A(:,k);
        Cin(k,:)=ci_raw;
        continue;
    end
    
    if and(and(ind_success_ci,deconv_flag),sn==1)
        % deconv the temporal trace
        [ci, ~, ~] = deconvolveCa(ci_raw, deconv_options_0, 'sn', sn);  % sn is 1 if Ysignal has no noise.
        % save this initialization
        Cin(k, :) = ci;       
%        Sin(k, :) = si; can get it later in iteration
%        kernel_pars{k} = reshape(deconv_options.pars, 1, []);
    else
        Cin(k, :)=ci_raw;
        ci=ci_raw;
    end      
    STD(k)=std(ci);
    
    % extract ai        
    % select its neighbours for estimation of ai, the box size is
    %[2*gSiz+1, 2*gSiz+1]
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);               % size of the box
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    HY_box = Ysignal(ind_nhood, :);      % extract temporal component from HY_box
    Y_box = Y(ind_nhood, :);
    Amask=Amask(ind_nhood,k);
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % index of the center

    sz = [nr, nc];
    [ai, ind_success_ai] = extract_a(ci, Y_box, HY_box, Amask, ind_ctr, sz, sn);  
    if and(ind_success_ai==false,sn==1)
        Ain(:,k)=A(:,k);
    elseif and(ind_success_ai==false,sn==[])
        del_ind(k)=true;
    else
        Ain(:,k)=ai;
    end
          
    % update the raw data
    Y(ind_nhood, :)=Y_box - ai*ci;
    if sn==1
        Ysignal(ind_nhood, :) = HY_box - ai*ci;            
    else
        rsub = max(1, -2*gSiz+r):min(d1, 2*gSiz+r);
        csub = max(1, -2*gSiz+c):min(d2, 2*gSiz+c);
        [cind, rind] = meshgrid(csub, rsub);
        ind_nhood_HY = sub2ind([d1, d2], rind(:), cind(:));
        [nr2, nc2] = size(cind);
        
        Hai = imfilter(reshape(Ain(ind_nhood_HY, k), nr2, nc2), psf, 'replicate');
        HY_box = HY(ind_nhood_HY, :) - Hai(:)*ci;
            % HY_box = bsxfun(@minus, HY_box, median(HY_box, 2));
        Ysignal(ind_nhood_HY, :) = HY_box;
    end
    
end


if sn==1
    ACS.Ain = [ACS.Ain sparse(Ain(:, ~del_ind))];
    ACS.Cin = [ACS.Cin; Cin(~del_ind, :)];
    ACS.STD=[ACS.STD STD];
else
    ACS.Ain = sparse(Ain(:, ~del_ind));
    ACS.Cin = Cin(~del_ind, :);
end
% if deconv_flag
%     File(i).Sin = Sin(~del_ind, :);
%     File(i).kernel_pars = cell2mat(kernel_pars(~del_ind));
end