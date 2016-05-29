function [Ain, Cin, center, Cn] = greedyROI_endoscope(Y, K, options,debug_on, save_avi)
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

%% Author: Pengcheng Zhou, Carnegie Mellon University.
% the method is an modification of greedyROI method used in Neuron paper of Eftychios
% Pnevmatikakis et.al. https://github.com/epnev/ca_source_extraction/blob/master/utilities/greedyROI2d.m
%% In each iteration of peeling off neurons, it searchs the one with maximum
% value of (max-median)/noise * Cn, which achieves a balance of SNR and
% local correlation.


%% use correlation to initialize NMF
%% parameters
d1 = options.d1;
d2 = options.d2;
gSig = options.gSig;
gSiz = options.gSiz;
min_corr = options.min_corr;
% boudnary to be removed 
try 
    bd = options.bd; 
catch 
    bd = gSig*2;
end
% min_pixel = 5;  % minimum number of pixels to be a neuron
sig = 2;    % thresholding noise by sig*std()
max_sig_ratio = 10; 

% maxIter = 5;
if ~exist('debug_on', 'var'); debug_on = false; end
if ~exist('save_avi', 'var'); save_avi = false; end
if ~ismatrix(Y); Y = reshape(Y, d1*d2, []); end;

%% preprocessing data
% spatially filter data
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate');
HY = reshape(HY, d1*d2, []);
T = size(HY, 2);
HY_max = max(HY, [], 2)-mean(HY, 2);    % maximum projection

% estimate noise level and thrshold diff(HY)
dHY = diff(HY, 1, 2);
Ysig = std(dHY(:, 1:10:end), 0, 2);
dHY(bsxfun(@lt, dHY, Ysig*sig)) = 0;    % all negative and noisy spikes are removed
HY_max(HY_max<= Ysig*max_sig_ratio) = 0;    %

% compute loal correlation
Cn = correlation_image(full(dHY), [1,3], d1,d2);
Cn0 = Cn;   % backup
Cn(isnan(Cn)) = 0;
Cn = Cn + rand(size(Cn))*(1e-6);
clear dHY;

% screen searching pixels as center of the neuron
v_search = Cn.*reshape(HY_max./Ysig, d1, d2);
v_search(or(Cn<min_corr, v_search<0)) = 0;
ind_search = false(d1*d2,1);  % showing whether this pixel has been searched before 

% show local correlation
if debug_on
    figure('position', [100, 100, 1290, 646]); %#ok<*UNRCH>
    subplot(231);
    imagesc(Cn); colorbar;
    axis equal off tight; hold on;
    title('correlation image');
    if save_avi
        avi_file = VideoWriter('greedyROI_example.avi');
        avi_file.open();
    end
end

%% start initialization
if ~exist('K', 'var')||isempty(K); K = round(d1*d2/25); end
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
center = zeros(K, 2);   % center of the initialized components

%% do initialization in a greedy way
searching_flag = true;
    k = 0;      %number of found components
while searching_flag   
    %% find local maximum as initialization point     
    %find all local maximum as initialization point
    tmp_d = round(gSiz/2) + 1;
    v_search(ind_search) = 0; 
    v_max = ordfilt2(v_search, tmp_d^2, true(tmp_d));
    % set boundary to be 0
    v_search(1:bd, :) = 0;
    v_search((end-bd+1):end, :) = 0;
    v_search(:, 1:bd) = 0;
    v_search(:, (end-bd):end) = 0;
    
    ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0));
    if(isempty(ind_localmax)); break; end 
    [max_vs, ind_sort] = sort(v_search(ind_localmax), 'descend');
    ind_localmax = ind_localmax(ind_sort);
    [r_peak, c_peak] = ind2sub([d1,d2],ind_localmax);
    
    %% try initialization over all local maximums 
    for mcell = 1:length(ind_localmax);
        % find the starting point
        ind_p = ind_localmax(mcell);
        max_v = max_vs(mcell);
        if max_v==0; % all pixels has been tried for initialization 
            searching_flag = false;
            break;
        end;  
        ind_search(ind_p) = true; % indicating that this pixel has been searched. 
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
            HY(ind_nhood, :) = HY_box - Hai(:)*ci;
            % update the maximum projection of HY
            dHY_box = diff(HY(ind_nhood, :), 1, 2);
            Ysig_box = std(dHY_box(:, 1:10:end), 0, 2);
            dHY_box(bsxfun(@lt, dHY_box, Ysig_box)) = 0;
            temp  = HY(ind_nhood, :);
            temp = max(temp, [], 2) - mean(temp, 2); 
            temp(temp<Ysig_box*max_sig_ratio) = 0;
            HY_max(ind_nhood, :) = temp;
            
            
            % update correlation image
            tmp_Cn = correlation_image(dHY_box, [1,3], nr, nc);
            tmp_Cn(or(isnan(tmp_Cn), tmp_Cn<min_corr)) = 0;
            Cn(ind_nhood) = tmp_Cn;
            
            % update search value
            v_search(ind_nhood) = Cn(ind_nhood).*HY_max(ind_nhood)./Ysig_box;
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
            plot(ci); title('temporal component');
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
if save_avi; avi_file.close(); end
