function [Ain, Cin, center, Cn, PNR] = greedyROI_endoscope(Y, K, options,debug_on, save_avi)
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
%   save_avi: save the video of initialization procedure

%% Output:
%       Ain:  d X K' matrix, estimated spatial component
%       Cin:  K'X T matrix, estimated temporal component
%       center: K' X 2, coordinate of each neuron's center
%       Cn:  d1*d2, correlation image

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
min_corr = options.min_corr;    %minimum local correlations for determining seed pixels
min_pnr = 15;               % peak to noise ratio for determining seed pixels
% boudnary to avoid for detecting seed pixels
try
    bd = options.bd;
catch
    bd = gSig*2;
end
% min_pixel = 5;  % minimum number of pixels to be a neuron
sig = 10;    % thresholding noise by sig*std()

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false;
elseif ischar(save_avi)
    avi_name = save_avi;
    save_avi = true;
    debug_on = true;
elseif save_avi
    avi_name = 'greedyROI_example.avi';
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
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

% filter the data
HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate');
HY = reshape(HY, d1*d2, []);
% HY_med = median(HY, 2);
% HY_max = max(HY, [], 2)-HY_med;    % maximum projection
HY = bsxfun(@minus, HY, median(HY, 2)); 
HY_max = max(HY, [], 2); 
Ysig = get_noise_fft(HY, options);
PNR = reshape(HY_max./Ysig, d1, d2);
PNR(PNR<min_pnr) = 0;
PNR0 = PNR;

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
Cn = Cn + rand(size(Cn))*(1e-6);
clear dHY;

% screen seeding pixels as center of the neuron
v_search = Cn.*PNR;
v_search(or(Cn<min_corr, PNR<min_pnr)) = 0;
ind_search = false(d1*d2,1);  % showing whether this pixel has been searched before
ind_search(v_search==0) = true; % ignore pixels with small correlations or low peak-noise-ratio

% show local correlation
if debug_on
    figure('position', [100, 100, 1290, 646]); %#ok<*UNRCH>
    subplot(231);
        imagesc(Cn); colorbar;
%     imagesc(Cn.*PNR, quantile(Cn(:).*PNR(:), [0.5, 0.99])); colorbar;
    axis equal off tight; hold on;
%     title('Cn * PNR');
    title('Cn'); 
    if save_avi
        avi_file = VideoWriter(avi_name);
        avi_file.open();
    end
end

%% start initialization
if ~exist('K', 'var')||isempty(K); K = sum(v_search(:)>0);
else
    K = min(sum(v_search(:)>0), K);
end
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
center = zeros(K, 2);   % center of the initialized components

%% do initialization in a greedy way
searching_flag = true;
k = 0;      %number of found components
% set boundary to be 0
ind_bd = false(size(v_search));
ind_bd(1:bd, :) = true;
ind_bd((end-bd+1):end, :) = true;
ind_bd(:, 1:bd) = true;
ind_bd(:, (end-bd):end) = true;
while searching_flag
    %% find local maximum as initialization point
    %find all local maximum as initialization point
    tmp_d = round(gSig)+1;
    v_search(ind_search) = 0;
    v_max = ordfilt2(v_search, tmp_d^2, true(tmp_d));
    % set boundary to be 0
    v_search(ind_bd) = 0;
    
    ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0));
    if(isempty(ind_localmax)); break; end
    [~, ind_sort] = sort(v_search(ind_localmax), 'descend');
    ind_localmax = ind_localmax(ind_sort);
    [r_peak, c_peak] = ind2sub([d1,d2],ind_localmax);
    
    %% try initialization over all local maximums
    for mcell = 1:length(ind_localmax);
        % find the starting point
        ind_p = ind_localmax(mcell);
        %         max_v = max_vs(mcell);
        max_v = v_search(ind_p);
        ind_search(ind_p) = true; % indicating that this pixel has been searched.
        if max_v==0; % all pixels have been tried for initialization
            continue;
        end;
        [r, c]  = ind2sub([d1, d2], ind_p);
        
        % roughly check whether this is a good starting point
        y0 = HY(ind_p, :);
        y0_std = std(diff(y0));
        y0(y0<median(y0)) = 0;
        if (k>=1) && any(corr(Cin(1:k, :)', y0')>0.9) %already found similar temporal traces
            continue;
        end
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
        
        %% show temporal trace in the center
        if debug_on
            subplot(232); cla;
            imagesc(reshape(v_search, d1, d2), [0, max_v]); colorbar;
            title(sprintf('neuron %d', k+1));
            axis equal off tight; hold on;
            plot(c_peak(mcell:end), r_peak(mcell:end), '.r');
            plot(c,r, 'og');
            subplot(233);
            imagesc(reshape(Cn(ind_nhood), nr, nc));
            axis equal off tight;
            title('correlation image');
            subplot(2,3,4:6); cla;
            plot(HY_box(ind_ctr, :)); title('activity in the center'); axis tight;
            if ~save_avi; pause; end
        end
        
        %% extract ai, ci
        sz = [nr, nc];
        [ai, ci, ind_success] =  extract_ac(HY_box, Y_box, ind_ctr, sz);
        if or(any(isnan(ai)), any(isnan(ci))); ind_success=false; end
        if max(ci)/get_noise_fft(ci)<min_pnr; ind_success=false; end
        if ind_success
            % save this initialization
            k = k+1;
            Ain(ind_nhood, k) = ai;
            Cin(k, :) = ci;
            center(k, :) = [r, c];
            
            % avoid searching nearby pixels that are highly correlated with this one
            ind_search(ind_nhood(ai>max(ai)*options.merge_thr)) = true;
            
            % update the raw data
            Y(ind_nhood, :) = Y_box - ai*ci;
            % update filtered data
            Hai = imfilter(reshape(ai, nr, nc), psf, 'replicate');  % filter ai
            HY_box = HY_box - Hai(:)*ci;       
%             HY_box_med = median(HY_box, 2);
            HY_box = bsxfun(@minus, HY_box, median(HY_box, 2)); 
            HY(ind_nhood, :) = HY_box;
            % update the maximum projection of HY
            %             dHY_box = diff(HY(ind_nhood, 1:nf:end), 1, 2);
            Ysig_box = get_noise_fft(HY_box, options);
            %             dHY_box(bsxfun(@lt, dHY_box, Ysig_box*sig)) = 0;
%             temp = max(HY_box, [], 2) - HY_box_med;
            temp = max(HY_box, [], 2);
            tmp_PNR = temp./Ysig_box;
            tmp_PNR(or(isnan(tmp_PNR), tmp_PNR<min_pnr)) = 0;
            PNR(ind_nhood) = tmp_PNR;
            HY_box_thr = HY_box;
%             HY_box_thr(bsxfun(@lt, HY_box, HY_box_med+Ysig_box*sig)) = 0;
            HY_box_thr(bsxfun(@lt, HY_box, Ysig_box*sig)) = 0;
            
            % update correlation image
            tmp_Cn = correlation_image(HY_box_thr, [1,2], nr, nc);
            tmp_Cn(or(isnan(tmp_Cn), tmp_Cn<min_corr)) = 0;
            Cn(ind_nhood) = tmp_Cn;
            
            % update search value
            v_search(ind_nhood) = Cn(ind_nhood).*PNR(ind_nhood);
            v_search(ind_bd) = 0;
            v_search(ind_search) = 0;
        else
            continue;
        end
        
        %% display results
        if debug_on
            subplot(231);
            plot(c, r, '.r');
            subplot(232);
            plot(c,r, 'or');
            subplot(233);
            imagesc(reshape(ai, nr, nc));
            axis equal off tight;
            title('spatial component');
            subplot(2,3,4:6); cla;
            plot(ci); title('temporal component'); axis tight;
            if save_avi;
                frame = getframe(gcf);
                frame.cdata = imresize(frame.cdata, [646, 1290]);
                avi_file.writeVideo(frame);
            else
                temp = input('type s to stop the debug mode:  ', 's');
                if strcmpi(temp, 's')
                    debug_on = false;
                end
            end
        end
        
        
        if mod(k, 10)==0
            fprintf('%d/%d neurons have been detected\n', k, K);
        end
        
        if k==K;
            searching_flag = false;
            break;
        end
    end
end
center = center(1:k, :);
Ain = sparse(Ain(:, 1:k));
Cin = Cin(1:k, :);
% Cin(Cin<0) = 0;
Cn = Cn0;
PNR = PNR0;
if save_avi; avi_file.close(); end
end