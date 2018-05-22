function [results, center, Cn, PNR, save_avi] = greedyROI_endoscope(Y, K, options,debug_on, save_avi)
%% a greedy method for detecting ROIs and initializing CNMF. in each iteration,
% it searches the one with large (peak-median)/noise level and large local
% correlation. It's the same with greedyROI_corr.m, but with some features
% specialized for endoscope data
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

%% Author: Pengcheng Zhou, Carnegie Mellon University. zhoupc1988@gmail.com
% the method is an modification of greedyROI method used in Neuron paper of Eftychios
% Pnevmatikakis et.al. https://github.com/epnev/ca_source_extraction/blob/master/utilities/greedyROI2d.m
% In each iteration of peeling off neurons, it searchs the one with maximum
% value of (max-median)/noise * Cn, which achieves a balance of SNR and
% local correlation.


%% use correlation to initialize NMF
%% parameters
d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSig = options.gSig;    % width of the gaussian kernel approximating one neuron
gSiz = options.gSiz;    % average size of neurons
ssub = options.ssub;
tsub = options.tsub;
if (ssub~=1) || (tsub~=1)
    d1_raw = d1;
    d2_raw = d2;
    T_raw = size(Y, ndims(Y));
    Y = dsData(reshape(double(Y), d1_raw, d2_raw, T_raw), options);
    [d1, d2, ~] = size(Y);
    options.d1 = d1;
    options.d2 = d2;
    gSig = gSig/ssub;
    gSiz = round(gSiz/ssub);
    options.min_pixel = options.min_pixel/(ssub^2);
    options.bd = round(options.bd/ssub);
    
end
min_corr = options.min_corr;    %minimum local correlations for determining seed pixels
min_pnr = options.min_pnr;               % peak to noise ratio for determining seed pixels
min_v_search = min_corr*min_pnr;
seed_method = options.seed_method; % methods for selecting seed pixels
if strcmpi(seed_method, 'manual')
    min_corr = min_corr/2;
    min_pnr = min_pnr/2;
end
% kernel_0 = options.kernel;
deconv_options_0= options.deconv_options;
min_pixel = options.min_pixel;  % minimum number of pixels to be a neuron
deconv_flag = options.deconv_flag;
% smin = options.smin;
% boudnary to avoid for detecting seed pixels
try
    bd = options.bd;
catch
    bd = round(gSiz/2);
end
sig = 3;    % thresholding noise by sig*std()

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false;
elseif ischar(save_avi)
    avi_name = save_avi;
    debug_on = true;
elseif save_avi
    debug_on = true; % turn on debug mode
else
    save_avi = false; %don't save initialization procedure
end
% debug mode and exporting results
if ~exist('debug_on', 'var')
    debug_on = false;
end

if ~ismatrix(Y); Y = reshape(Y, d1*d2, []); end;  % convert the 3D movie to a matrix
Y(isnan(Y)) = 0;    % remove nan values
Y = double(Y);
T = size(Y, 2);

%% preprocessing data
% create a spatial filter for removing background
if gSig>0
    if options.center_psf
        psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
        ind_nonzero = (psf(:)>=max(psf(:,1)));
        psf = psf-mean(psf(ind_nonzero));
        psf(~ind_nonzero) = 0;
    else
        psf = fspecial('gaussian', round(gSiz), gSig);
    end
else
    psf = [];
end

% filter the data
if isempty(psf)
    % no filtering
    HY = Y;
else
    HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate');
end

HY = reshape(HY, d1*d2, []);
% HY_med = median(HY, 2);
% HY_max = max(HY, [], 2)-HY_med;    % maximum projection
HY = bsxfun(@minus, HY, median(HY, 2));
HY0=HY; % saved as intact filtered result.
HY_max = max(HY, [], 2);
Ysig = GetSn(HY);
PNR = reshape(HY_max./Ysig, d1, d2);
PNR0 = PNR;
PNR(PNR<min_pnr) = 0;

% estimate noise level and thrshold diff(HY)
% dHY = diff(HY(:, 1:nf:end), 1, 2);  %
% Ysig = std(dHY(:, 1:5:end), 0, 2);
% dHY(bsxfun(@lt, dHY, Ysig*sig)) =0;    % all negative and noisy spikes are removed
HY_thr = HY;
HY_thr(bsxfun(@lt, HY_thr, Ysig*sig)) = 0;

% compute loal correlation
Cn = correlation_image(HY_thr, [1,2], d1,d2);
Cn0 = Cn;   % backup
Cn(isnan(Cn)) = 0;
% Cn = Cn + randn(size(Cn))*(1e-100);

% screen seeding pixels as center of the neuron
v_search = Cn.*PNR;
v_search(or(Cn<min_corr, PNR<min_pnr)) = 0;
ind_search = false(d1*d2,1);  % showing whether this pixel has been searched before
ind_search(v_search==0) = true; % ignore pixels with small correlations or low peak-noise-ratio

% ignore boundaries pixels when determinging seed pixels
if length(bd) ==1
    bd = ones(1,4)*bd;
end
ind_bd = false(size(v_search));
ind_bd(1:bd(1), :) = true;
ind_bd((end-bd(2)+1):end, :) = true;
ind_bd(:, 1:bd(3)) = true;
ind_bd(:, (end-bd(4)+1):end) = true;

% show local correlation
if debug_on
    figure('position', [100, 100, 1200, 800], 'color', [1,1,1]*0.9); %#ok<*UNRCH>
    set(gcf, 'defaultAxesFontSize', 20);
    ax_cn = axes('position', [0.04, 0.5, 0.3, 0.4]);
    ax_pnr_cn = axes('position', [0.36, 0.5, 0.3, 0.4]);
    ax_cn_box = axes('position', [0.68, 0.54, 0.24, 0.32]);
    ax_trace = axes('position', [0.05, 0.05, 0.92, 0.4]);
    axes(ax_cn);
    imagesc(Cn0);
    %     imagesc(Cn.*PNR, quantile(Cn(:).*PNR(:), [0.5, 0.99]));
    axis equal off; hold on;
    axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
    %     title('Cn * PNR');
    title('Cn');
    if exist('avi_name', 'var')
        avi_file = VideoWriter(avi_name);
        avi_file.FrameRate = 1;
        avi_file.open();
    elseif save_avi
        avi_file = VideoWriter('initialization.avi');
        avi_file.FrameRate = 1;
        avi_file.open();
    end
    
end

%% start initialization
if ~exist('K', 'var')||isempty(K)
    K = floor(sum(v_search(:)>0)/10);
else
    K = min(floor(sum(v_search(:)>0)/10), K);
end
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
Sin = zeros(K, T);    % spike counts
Cin_raw = zeros(K, T);
kernel_pars = cell(K,1);    % parameters for the convolution kernels of all neurons
center = zeros(K, 2);   % center of the initialized components

%% do initialization in a greedy way
searching_flag = true;
k = 0;      %number of found components
[ii, jj] = meshgrid(1:d2, 1:d1);
pixel_v = (ii*10+jj)*(1e-10);
while searching_flag
    %% find local maximum as initialization point
    %find all local maximum as initialization point
    
    tmp_d = max(3, round(gSiz/4));
    v_search = medfilt2(v_search,3*[1, 1])+pixel_v; % add an extra value to avoid repeated seed pixels within one ROI.
    v_search(ind_search) = 0;
    v_max = ordfilt2(v_search, tmp_d^2, true(tmp_d));
    % set boundary to be 0
    v_search(ind_bd) = 0;
    
    if strcmpi(seed_method, 'manual') %manually select seed pixels
        tmp_fig = figure('position', [200, 200, 1024, 412]);
        subplot(121); cla;
        imagesc(PNR);  hold on;
        title('PNR');
        plot(center(1:k, 2), center(1:k, 1), '*r');
        axis equal off;
        axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
        colorbar;
        
        subplot(122);
        imagesc(Cn, [min_corr, 1]); %, [0, max(max(min_v_search(:)*0.99), min_v_search)]);
        colorbar;
        hold on;
        axis equal;
        axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
        drawnow;
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        title('click neuron centers for initialziation');
        xlabel('click invalid pixels to stop', 'color', 'r');
        ind_localmax = zeros(K,1);
        for tmp_k=1:K
            figure(tmp_fig);
            [tmp_x, tmp_y] = ginput(1);
            tmp_x = round(tmp_x); tmp_y = round(tmp_y);
            if isempty(tmp_x)||or(tmp_x<1, tmp_x>d2) || or(tmp_y<1, tmp_y>d1) ||(v_search(tmp_y, tmp_x)==0)
                break;
            end
            plot(tmp_x, tmp_y, '*r', 'linewidth', 2);
            drawnow();
            ind_localmax(tmp_k) = sub2ind([d1,d2], tmp_y, tmp_x);
        end
        close(tmp_fig);
        
        ind_localmax = ind_localmax(1:(tmp_k-1));
        if isempty(ind_localmax)
            break;
        end
    else
        % automatically select seed pixels
        ind_search(v_search<min_v_search) = true;    % avoid generating new seed pixels after initialization
        ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0));
        if(isempty(ind_localmax)); break; end
    end
    [~, ind_sort] = sort(v_search(ind_localmax), 'descend');
    ind_localmax = ind_localmax(ind_sort);
    [r_peak, c_peak] = ind2sub([d1,d2],ind_localmax);
    
    %% try initialization over all local maximums
    for mcell = 1:length(ind_localmax)
        % find the starting point
        ind_p = ind_localmax(mcell);
        %         max_v = max_vs(mcell);
        max_v = v_search(ind_p);
        if mcell==1
            img_clim = [0, max_v];
        end
        ind_search(ind_p) = true; % indicating that this pixel has been searched.
        if max_v<min_v_search % all pixels have been tried for initialization
            continue;
        end
        [r, c]  = ind2sub([d1, d2], ind_p);
        
        % roughly check whether this is a good starting point
        y0 = HY(ind_p, :);
        y0_std = std(diff(y0));
        %         y0(y0<median(y0)) = 0;
        %         if (k>=1) && any(corr(Cin(1:k, :)', y0')>0.9) %already found similar temporal traces
        %             continue;
        %         end
        if max(diff(y0))< 3*y0_std % signal is weak
            continue;
        end
        
        % select its neighbours for estimation of ai and ci, the box size is
        %[2*gSiz+1, 2*gSiz+1]
        rsub = max(1, -gSiz+r):min(d1, gSiz+r);
        csub = max(1, -gSiz+c):min(d2, gSiz+c);
        [cind, rind] = meshgrid(csub, rsub);
        [nr, nc] = size(cind);
        ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
        HY_box = HY(ind_nhood, :);      % extract temporal component from HY_box
        Y_box = Y(ind_nhood, :);    % extract spatial component from Y_box
        ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
        
        % neighbouring pixels to update after initialization of one
        % neuron
        rsub = max(1, -2*gSiz+r):min(d1, 2*gSiz+r);
        csub = max(1, -2*gSiz+c):min(d2, 2*gSiz+c);
        [cind, rind] = meshgrid(csub, rsub);
        ind_nhood_HY = sub2ind([d1, d2], rind(:), cind(:));
        [nr2, nc2] = size(cind);
        
        %% show temporal trace in the center
        if debug_on
            axes(ax_pnr_cn); cla;
            imagesc(reshape(v_search, d1, d2), img_clim); % [0, max_v]);
            title(sprintf('neuron %d', k+1));
            axis equal off; hold on;
            axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
            plot(c_peak(mcell:end), r_peak(mcell:end), '.r');
            plot(c,r, 'or', 'markerfacecolor', 'r', 'markersize', 10);
            axes(ax_cn_box);
            imagesc(reshape(Cn(ind_nhood), nr, nc), [0, 1]);
            axis equal off tight;
            title('correlation image');
            axes(ax_trace); cla;
            plot(HY_box(ind_ctr, :)); title('activity in the center'); axis tight;
            if ~save_avi; pause; end
            if exist('avi_file', 'var')
                frame = getframe(gcf);
                frame.cdata = imresize(frame.cdata, [800, 1200]);
                avi_file.writeVideo(frame);
            end
        end
        
        %% extract ai, ci
        sz = [nr, nc];
        if options.center_psf
            [ai, ci_raw, ind_success] =  extract_ac(HY_box, Y_box, ind_ctr, sz, options.spatial_constraints);
        else
            [ai, ci_raw, ind_success] =  extract_ac(HY_box, Y_box, ind_ctr, sz, options.spatial_constraints);
        end
        if or(any(isnan(ai)), any(isnan(ci_raw))); ind_success=false; end
        if sum(ai)<=min_pixel; ind_success = false; end
        %         if max(ci_raw)<min_pnr;
        %             ind_success=false;
        %         end
        if sum(ai(:)>0)<min_pixel; ind_success=false; end
        if ind_success
            k = k+1;
            
            if deconv_flag
                % deconv the temporal trace
                [ci, si, deconv_options] = deconvolveCa(ci_raw, deconv_options_0);  % sn is 1 because i normalized c_raw already
                % save this initialization
                Ain(ind_nhood, k) = ai;
                Cin(k, :) = ci;
                Sin(k, :) = si;
                Cin_raw(k, :) = ci_raw-deconv_options.b;
                %                 kernel_pars(k, :) = kernel.pars;
                kernel_pars{k} = reshape(deconv_options.pars, 1, []);
            else
                ci = ci_raw;
                Ain(ind_nhood, k) = ai;
                Cin(k, :) = ci_raw;
                Cin_raw(k, :) = ci_raw;
            end
            ci = reshape(ci, 1,[]);
            center(k, :) = [r, c];
            
            % avoid searching nearby pixels
            ind_search(ind_nhood(ai>max(ai)*0.5)) = true;
            
            % update the raw data
            Y(ind_nhood, :) = Y_box - ai*ci;
            % update filtered data
            if isempty(psf)
                Hai = reshape(Ain(ind_nhood_HY, k), nr2, nc2);
            else
                Hai = imfilter(reshape(Ain(ind_nhood_HY, k), nr2, nc2), psf, 'replicate');
            end
            HY_box = HY(ind_nhood_HY, :) - Hai(:)*ci;
            %             HY_box = bsxfun(@minus, HY_box, median(HY_box, 2));
            HY(ind_nhood_HY, :) = HY_box;
            
            % update the maximum projection of HY
            Ysig_box = Ysig(ind_nhood_HY);
            temp = max(HY_box, [], 2);
            tmp_PNR = temp./Ysig_box;
            tmp_PNR(or(isnan(tmp_PNR), tmp_PNR<min_pnr)) = 0;
            PNR(ind_nhood_HY) = tmp_PNR;
            
            HY_box_thr = HY_box;  %thresholded version of HY
            HY_box_thr(bsxfun(@lt, HY_box, Ysig_box*sig)) = 0;
            
            % update correlation image
            tmp_Cn = correlation_image(HY_box_thr, [1,2], nr2, nc2);
            tmp_Cn(or(isnan(tmp_Cn), tmp_Cn<min_corr)) = 0;
            Cn(ind_nhood_HY) = tmp_Cn;
            
            % update search value
            v_search(ind_nhood_HY) = Cn(ind_nhood_HY).*PNR(ind_nhood_HY);
            v_search(ind_bd) = 0;
            v_search(ind_search) = 0;
        else
            continue;
        end
        
        %% display results
        if debug_on
            axes(ax_cn);
            plot(c, r, '.r');
            axes(ax_pnr_cn);
            plot(c,r, 'or');
            axes(ax_cn_box);
            imagesc(reshape(ai, nr, nc));
            axis equal off tight;
            title('spatial component');
            axes(ax_trace);  cla; hold on;
            plot(ci_raw, 'linewidth', 2); title('temporal component'); axis tight;
            if deconv_flag
                plot(ci, 'r', 'linewidth', 2);  axis tight;
                legend('raw trace', 'denoised trace');
            end
            if exist('avi_file', 'var')
                frame = getframe(gcf);
                frame.cdata = imresize(frame.cdata, [800, 1200]);
                avi_file.writeVideo(frame);
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
            fprintf('%d neurons have been detected\n', k);
        end
        
        if k==K
            searching_flag = false;
            break;
        end
    end
end
center = center(1:k, :);
results.Ain = Ain(:, 1:k);
results.Cin = Cin(1:k, :);
results.Cin_raw = Cin_raw(1:k, :);
if deconv_flag
    results.Sin = Sin(1:k, :);
    results.kernel_pars = reshape(cell2mat(kernel_pars(1:k)), k, []);
end
% Cin(Cin<0) = 0;
Cn = Cn0;
PNR = PNR0;
if ssub~=1
    Ain = reshape(full(results.Ain), d1,d2, []);
    if ~isempty(Ain)
        results.Ain = sparse(reshape(imresize(Ain, [d1_raw, d2_raw]), d1_raw*d2_raw, []));
    else
        results.Ain = sparse(zeros(d1_raw*d2_raw, 0));
    end
    Cn =imresize(Cn, [d1_raw, d2_raw]);
    PNR = imresize(PNR, [d1_raw, d2_raw]);
    center = center*ssub-1; 
end
if tsub~=1
    if k==0
        results.Cin = zeros(0, T_raw);
        results.Cin_raw = zeros(0, T_raw);
        results.Sin = zeros(0, T_raw);
    else
        results.Cin = imresize(results.Cin, [k, T_raw], 'box');
        results.Cin_raw = imresize(results.Cin_raw, [k, T_raw], 'box');
        if deconv_flag
            results.Sin = imresize(results.Sin, [k, T_raw]);
        end
    end
end
if exist('avi_file', 'var')
    close(gcf);
    if avi_file.Duration==0
        warning('off', 'MATLAB:audiovideo:VideoWriter:noFramesWritten')
        avi_file.close();
        delete(avi_name);
    else
        avi_file.close();
    end
end
end
