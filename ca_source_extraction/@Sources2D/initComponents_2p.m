function center = initComponents_2p(obj, Y, K, debug_on, save_avi)
%% a greedy method for detecting ROIs and initializing CNMF. in each iteration,
% it searches the one with large (peak-median)/noise level and large local
% correlation
%% Input:
%   Y:  d X T matrx, imaging data
%   K:  scalar, maximum number of neurons to be detected.
%   debug_on: options for showing procedure of detecting neurons
%% Output:
%       Ain:  d X K' matrix, estimated spatial component
%       Cin:  K'X T matrix, estimated temporal component
%       bin:  d X nb matrix/vector, spatial components of the background
%       Cin:  nb X T matrix/vector, temporal components of the background
%       center: K' X 2, coordinate of each neuron's center
%       res:  d X T, residual after initializing Ain, Cin, bin, fin

%% Author: Pengcheng Zhou, Carnegie Mellon University.
% the method is an modification of greedyROI method used in Neuron paper of Eftychios
% Pnevmatikakis et.al. https://github.com/epnev/ca_source_extraction/blob/master/utilities/greedyROI2d.m
%% In each iteration of initializing neurons, it searchs the one with maximum
% value of (max-median)/noise * Cn, which selects pixels with high SNR and
% local correlation.

%% parameters
if exist('K', 'var')
    K = 200; 
end
Y = obj.reshape(Y, 1); 
Y_std = get_noise_fft(Y);
options = obj.options;

if ~exist('debug_on', 'var'); debug_on = false; end
if ~exist('save_avi', 'var'); save_avi=false; end
d1 = options.d1;
d2 = options.d2;
gSig = options.gSig;
gSiz = options.gSiz;
if and(isempty(gSiz), isempty(gSig)); gSig = 3; gSiz = 10; end
if isempty(gSiz); gSiz=3*gSig; end
if isempty(gSig); gSig=gSiz/3; end
if isfield(options, 'min_corr')
    min_corr = options.min_corr;    % minimum local correaltion value to start one neuron
else
    min_corr = 0.3;
end
if isfield(options, 'deconv_flag')
    deconv_flag = options.deconv_flag;
    deconv_options_0= options.deconv_options;
else
    deconv_flag = false; 
end 
if isfield(options, 'min_pnr')
    min_pnr = options.min_pnr;
else
    min_pnr = 2;
end
nb = options.nb;        % number of the background
pSiz = 1;       % after selecting one pixel, take the mean of square box
%near the pixel as temporal activity. the box size is (2*pSiz+1)
psf = ones(gSig)/(gSig^2);

maxIter = 5;            % iterations for refining results
sz = 4;            %distance of neighbouring pixels for computing local correlation

if ~ismatrix(Y); Y = reshape(Y, d1*d2, []); end;
[~, T] = size(Y);       % number of frames
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
Cin_raw = zeros(K, T);      % temporal components
Sin = zeros(K, T);      % temporal components
center = zeros(K, 2);   % center of the initialized components
kernel_pars = cell(K,1);    % parameters for the convolution kernels of all neurons

%% compute correlation image and (max-median)/std ratio
ind_frame = round(linspace(1, T, min(T, 1000)));    % select few frames for speed issues
%tmp_noise = randn(d1*d2, length(ind_frame)); 
C1 = correlation_image(full(Y(:, ind_frame)), sz, d1, d2);
Cb =  zeros(size(C1)); %correlation_image(full(Y(:, ind_frame(1:3:end)))+tmp_noise(:, 1:3:end), [gSiz, gSiz+1], d1, d2);  %backgroung correlatin 
Cn = C1-Cb; %here Cb is the background correlation. for 2photon imaging results. It might be useful when the background signal is large 
Y_median = median(Y(:, ind_frame), 2);
Y = bsxfun(@minus, Y, Y_median);
% Y_std = sqrt(mean(Y.*Y, 2));

%% find local maximum
k = 0;      %number of found components
min_pixel = floor(gSig^2/2);  % minimum number of peaks to be a neuron
peak_ratio = full(max(Y, [], 2))./Y_std; %(max-median)/std
peak_ratio(isinf(peak_ratio)) = 0;  % avoid constant values
peak_ratio(peak_ratio<min_pnr) = 0; 

% save_avi = false;   %save procedures for demo
if debug_on
    figure('position', [100, 100, 800, 650]); %#ok<*UNRCH>
      figure('position', [100, 100, 1200, 800], 'color', [1,1,1]*0.9); %#ok<*UNRCH>
    set(gcf, 'defaultAxesFontSize', 20); 
    ax_cn = axes('position', [0.04, 0.5, 0.3, 0.4]); 
    ax_cn_varying = axes('position', [0.36, 0.5, 0.3, 0.4]); 
    ax_cn_box = axes('position', [0.68, 0.54, 0.24, 0.32]); 
    ax_raw = axes('position', [0.05, 0.25, 0.92, 0.2]); 
    ax_trace = axes('position', [0.05, 0.01, 0.92, 0.2]); 
    axes(ax_cn); 
    imagesc(Cn);  hold on;  
    axis equal off tight; 
    title('correlation image'); 
    axes(ax_cn_varying);
    imagesc(Cn, [0,1]);
    axis equal off tight; hold on;
    if save_avi
        avi_file = VideoWriter('greedyROI_example.avi');
        avi_file.open();
    end
end

max_thresh = min_pnr * (min_corr);
while k<K
    %% find the pixel with the maximum ratio
    [max_v, ind_p] = max(peak_ratio.*(Cn(:)));
    peak_ratio(ind_p) = 0;  % no longer visit this pixel any more
    if max_v<max_thresh; break; end
    if Cn(ind_p)<min_corr; continue; end % ignore this local maximum due to small local correlation
    if  max_v/(min_corr)< min_pnr;     continue;    end
    
    [r, c] = ind2sub([d1,d2], ind_p);
    
    % select its neighbours for computing correlation
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    nr = length(rsub);  %size of the neighboring matrix
    nc = length(csub);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    Y_box = Y(ind_nhood, :);
    
    % draw a small area near the peak and extract the mean activities
    r0 = rsub(1); c0 = csub(1);
    rsub = (max(1, -pSiz+r):min(d1, pSiz+r)) - r0+1;
    csub = (max(1, -pSiz+c):min(d2, pSiz+c)) -c0+1;
    [cind, rind] = meshgrid(csub, rsub);
    ind_peak = sub2ind([nr, nc], rind(:), cind(:));
    y0 = mean(Y_box(ind_peak, :), 1);
%     y0(y0<0) = 0;
    
    % compute the correlation between the peak and its neighbours
    temp = reshape(corr(y0', Y_box'), nr, nc);
    active_pixel = full(temp>min_corr/2);
    l = bwlabel(active_pixel, 8);   % remove disconnected components
    active_pixel(l~=mode(l(ind_peak))) = false;
    tmp_v = sum(active_pixel(:));    %number of pixels with above-threshold correlation
    if debug_on
        axes(ax_cn_varying); cla;
        imagesc(reshape(Cn(:), d1, d2), [0, Cn(ind_p)]); 
        title(sprintf('neuron %d', k+1));
        axis equal off tight; hold on;
        plot(c,r, 'om');
        axes(ax_cn_box);
        imagesc(temp, [min_corr, 1]);
        axis equal off tight;
        title('local corr.');
        axes(ax_raw); cla; hold on; 
        plot(y0/max(y0)); title('activity in the center');
        axis off tight; 
        if ~save_avi; pause; end
    end
    if tmp_v<min_pixel;         continue;  end % neuron is too small
    
    %% save neuron
    %   nonzero area
    data = Y_box(active_pixel(:), :);
    ind_active = ind_nhood(active_pixel(:));  %indices of active pixels within the whole frame
    peak_ratio(ind_nhood(ind_peak)) = 0;    % the small area near the peak is not able to initialize neuron anymore
    
    % do a rank-1 matrix factorization in this small area
    [ai, ci_raw] = finetune2d(data, y0, maxIter);
    if norm(ai)==0;        continue;   end
    k = k+1;
    if deconv_flag
        % deconv the temporal trace
        [ci, si, deconv_options] = deconvolveCa(ci_raw, deconv_options_0);  % sn is 1 because i normalized c_raw already
        % save this initialization
        Ain(ind_active, k) = ai;
        Cin(k, :) = ci;
        Sin(k, :) = si;
        Cin_raw(k, :) = ci_raw;
        %                 kernel_pars(k, :) = kernel.pars;
        kernel_pars{k} = reshape(deconv_options.pars, 1, []);
    else
        ci = ci_raw;
        Ain(ind_active, k) = ai;
        Cin(k, :) = ci_raw;
        Cin_raw(k, :) = ci_raw;
    end
    ci = reshape(ci, 1, []); 
    Y(ind_active, :) = data-ai*ci;
    center(k, :) = [r, c];
    
    if debug_on
        axes(ax_cn);
        plot(c, r, '.r');
        axes(ax_cn_varying);
        plot(c,r, 'or');
        axes(ax_cn_box);
        temp = zeros(nr, nc); temp(active_pixel) = ai;
        imagesc(temp);
        axis equal off tight;
        title('spatial component');
        axes(ax_trace);  cla; hold on;
        if deconv_flag
                    plot(ci_raw, 'linewidth', 2); title('temporal component'); axis tight;
            plot(ci, 'r', 'linewidth', 1);  axis tight;
            legend('raw trace', 'denoised trace');
        else
            plot(ci_raw, 'r', 'linewidth', 1); title('temporal component'); axis tight;
        end
        title('temporal component');
        if exist('avi_file', 'var')
            temp = getframe(gcf); 
            temp.cdata = imresize(temp.cdata, [800, 640]); 
            avi_file.writeVideo(temp); 
        elseif ~save_avi
            temp = input('type s to stop the debug mode:  ', 's');
            if strcmpi(temp, 's')
                save_avi = true;
            end
        else
            drawnow(); 
        end
    end
    
    if mod(k, 10)==0
        fprintf('%d/%d neurons have been detected\n', k, K);
    end
    
    if k==K;   break; end
    
    %% udpate peak_ratio and correlation image
    tmp_old = peak_ratio(ind_active);
    tmp_new = max(Y(ind_active, :), [], 2)./Y_std(ind_active);
    temp = zeros(nr, nc);
    temp(active_pixel) = max(0, tmp_old-tmp_new); % after each iteration, the peak ratio can not be increased
    peak_ratio(ind_nhood) = max(0, peak_ratio(ind_nhood) - reshape(imfilter(temp, psf), [], 1)); % update peak_ratio, results are smoothed
    Cn(ind_nhood) = correlation_image(full(Y(ind_nhood, ind_frame)), sz, nr, nc)-reshape(Cb(ind_nhood), nr, nc);  % update local correlation
end

center = center(1:k, :);
obj.A = sparse(Ain(:, 1:k));
obj.C = Cin(1:k, :);
obj.C_raw = Cin_raw(1:k, :);
obj.S = Sin(1:k, :);
obj.P.kernel_pars = kernel_pars(1:k); 

if exist('avi_file', 'var'); avi_file.close(); end
res = bsxfun(@plus, Y, Y_median);

% %% initialize background
% tsub = max(1, round(T/1000));
% [bin, f] = nnmf(max(res(:, 1:tsub:T), 0), nb);
% fin = imresize(f, [nb, T]);
% fin = HALS_temporal(max(res, 0), bin, fin, maxIter);
% obj.b = HALS_spatial(max(res, 0), bin, fin, [], maxIter);
% obj.f = fin; 
end

function [ai, ci] = finetune2d(data, ci, nIter)
%do matrix factorization given the model data = ai*ci, where ai>=0
%
%Input:
%   data:   d x T matrix, small patch containing one neuron
%   ci:     initial value for trace
%   nIter  number of coordinate descent steps
%
%Output:
%   ai  M x N matrix, result of the fine-tuned neuron shape
%   ci  1 x T matrix, result of the neuron
%% copied from greedyROI.m

if ~exist('nIter', 'var'), nIter = 1; end
data(data<0)= 0;
%do block coordinate descent
for iter = 1:nIter,
    %update basis
    ai = max(0, (data*ci')/(ci*ci'));
    norm_ai = norm(ai, 2);
    if norm_ai==0; break;     end
    ai = ai/norm_ai;
    ci =  (ai'*data);
    %     ci(ci<0) = 0;
end
% [b, sn] = estimate_baseline_noise(ci); 
% ci = ci - b; 
% ai = ai*sn; 
% ci = ci/sn; 
end