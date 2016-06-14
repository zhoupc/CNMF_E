classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data.
    % Author: Pengcheng Zhou, zhoupc1988@gmail.com with modifications from
    % Eftychios Pnevmatikakis
    
    %% properties
    properties
        A;          % spatial components of neurons
        C;          % temporal components of neurons
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        S;          % spike counts
        Cn;         % correlation image
        Coor;       % neuron contours
        Df;         % background for each component to normalize the filtered raw data
        C_df;       % temporal components of neurons and background normalized by Df
        S_df;       % spike counts of neurons normalized by Df
        options;    % options for model fitting
        P;          % some estimated parameters
        Fs = nan;    % frame rate
        indicator = 'GCaMP6f'; 
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = Sources2D(varargin)
            obj.options = CNMFSetParms();
            obj.P = struct('p', 2);
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:});
            end
        end
        
        %% update parameters
        function updateParams(obj, varargin)
            obj.options = CNMFSetParms(obj.options, varargin{:});
        end
        
        %% data preprocessing
        function Y = preprocess(obj,Y,p)
            [obj.P,Y] = preprocess_data(Y,p,obj.options);
        end
        
        %% fast initialization
        function [center] = initComponents(obj, Y, K, tau)
            if nargin<4 ;    tau = [];             end
            [obj.A, obj.C, obj.b, obj.f, center] = initialize_components(Y, K, tau, obj.options);
        end
        
        %% fast initialization for microendoscopic data
        [center, Cn, pnr] = initComponents_endoscope(obj, Y, K, patch_sz, debug_on, save_avi);
        
        %% update spatial components
        function updateSpatial(obj, Y)
            [obj.A, obj.b, obj.C] = update_spatial_components(Y, ...
                obj.C, obj.f, obj.A, obj.P, obj.options);
        end
        
        %% udpate spatial components without background
        function updateSpatial_nb(obj, Y)
            [obj.A, obj.C] = update_spatial_components_nb(Y, ...
                obj.C, obj.A, obj.P, obj.options);
        end
        %% update temporal components
        function updateTemporal(obj, Y)
            [obj.C, obj.f, obj.P, obj.S] = update_temporal_components(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        %% update temporal components without background
        function updateTemporal_nb(obj, Y)
            [obj.C, obj.P, obj.S] = update_temporal_components_nb(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        %% merge found components
        function [nr, merged_ROIs] = merge(obj, Y)
            [obj.A, obj.C, nr, merged_ROIs, obj.P, obj.S] = merge_components(...
                Y,obj.A, [], obj.C, [], obj.P,obj.S, obj.options);
        end
        
        %% compute the residual
        function [Y_res] = residual(obj, Yr)
            Y_res = Yr - obj.A*obj.C - obj.b*obj.f;
        end
        
        %% take the snapshot of current results
        function [A, C,  b, f, P, S] = snapshot(obj)
            A = obj.A;
            C = obj.C;
            b = obj.b;
            f = obj.f;
            P = obj.P;
            try
                S = obj.S;
            catch
                S = [];
            end
        end
        
        %% extract DF/F signal after performing NMF
        function [C_df, Df, S_df] = extractDF_F(obj, Y, i)
            if ~exist('i', 'var')
                i = size(obj.A, 2) + 1;
            end
            
            [obj.C_df, obj.Df, obj.S_df] = extract_DF_F(Y, [obj.A, obj.b],...
                [obj.C; obj.f], obj.S, i);
            
            C_df =  obj.C_df;
            Df = obj.Df;
            S_df = obj.S_df;
            
        end
        
        %% order_ROIs
        function [srt] = orderROIs(obj, srt)
            %% order neurons
            % srt: sorting order
            
            if nargin<2; srt=[]; end
            [obj.A, obj.C, obj.S, obj.P, srt] = order_ROIs(obj.A, obj.C,...
                obj.S, obj.P, srt);
        end
        
        %% view contours
        function [Coor, json_file] = viewContours(obj, Cn, contour_threshold, display, ind)
            if or(isempty(Cn), ~exist('Cn', 'var') )
                Cn = reshape(obj.P.sn, obj.options.d1, obj.options.d2);
            end
            if nargin<4; display=0; end
            if nargin<5; ind=1:size(obj.A, 2); end
            [obj.Coor, json_file] = plot_contours(obj.A(:, ind), Cn, ...
                contour_threshold, display);
            Coor = obj.Coor;
        end
        
        %% plot components
        function plotComponents(obj, Y, Cn)
            if ~exist('Cn', 'var')
                Cn = [];
            end
            view_components(Y, obj.A, obj.C, obj.b, obj.f, Cn, obj.options);
        end
        
        %% plot components GUI
        function plotComponentsGUI(obj, Y, Cn)
            if ~exist('Cn', 'var')
                Cn = [];
            end
            plot_components_GUI(Y,obj.A,obj.C,obj.b,obj.f,Cn,obj.options)
        end
        
        %% make movie
        function makePatchVideo(obj, Y)
            make_patch_video(obj.A, obj.C, obj.b, obj.f, Y, obj.Coor,...
                obj.options);
        end
        
        %% new methods added by PC, Since 02/05/2016
        %% down sample data and initialize it
        [obj_ds, Y_ds] = downSample(obj, Y)
        
        %% up-sample the results
        upSample(obj, obj_ds, T);
        
        %% copy the objects
        function obj_new = copy(obj)
            % Instantiate new object of the same class.
            obj_new = feval(class(obj));
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                obj_new.(p{i}) = obj.(p{i});
            end
        end
        
        %% quick merge neurons based on spatial and temporal correlation
        merged_ROIs = quickMerge(obj, temporal_component)
        
        
        %% quick view
        viewNeurons(obj, ind, C2, folder_nm);
        displayNeurons(obj, ind, C2, folder_nm);
        
        %% delete neurons
        function delete(obj, ind)
            obj.A(:, ind) = [];
            obj.C(ind, :) = [];
            try  %#ok<TRYNC>
                obj.S(ind, :) = [];
                obj.P.gn(ind) = [];
                obj.P.b(ind) = [];
                obj.P.c1(ind) = [];
                obj.P.neuron_sn(ind) = [];
                obj.Coor(ind) = [];
            end
        end
        
        %% estimate neuron centers
        function [center] = estCenter(obj)
            center = com(obj.A, obj.options.d1, obj.options.d2);
        end
        
        %% update A & C using HALS
        function [obj, IDs] = HALS_AC(obj, Y)
            %update A,C,b,f with HALS
            Y = obj.reshape(Y, 1);
            [obj.A, obj.C, obj.b, obj.f, IDs] = HALS_2d(Y, obj.A, obj.C, obj.b,...
                obj.f, obj.options);
        end
        
        %% update backgrounds
        [B, F] = updateBG(obj, Y, nb, model)
        
        %% reshape data
        function Y = reshape(obj, Y, dim)
            % reshape the imaging data into diffrent dimensions
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            if dim==1; Y=reshape(Y, d1*d2, []);  %each frame is a vector
            else  Y = reshape(full(Y), d1, d2, []);    %each frame is an image
            end
        end
        
        %% deconvolve all temporal components
        function C0 = deconvTemporal(obj)
            C0 = obj.C;
            [obj.C, obj.P, obj.S] = deconv_temporal(obj.C, obj.P, obj.options);
        end
        
        %% update background
        function [Y, ind_bg, bg] = linearBG(obj, Y)
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            % model the background as linear function
            if ndims(Y)==3; Y = neuron.reshape(Y,1); end
            T = size(Y,2);
            % find background area
            ind_frame = round(linspace(1, T, min(T, 500)));
            tmp_C1 = correlation_image(Y(:, ind_frame), [1, 2], d1, d2);
            tmp_C2 = correlation_image(Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
            tmp_Cn = tmp_C1(:)-tmp_C2(:);
            ind = (tmp_Cn>0);
            ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), 0.1));
            % get the mean activity within the selected area
            Ybg = mean(Y(ind_bg(:), :), 1);
            f1 = (Ybg-mean(Ybg))/std(Ybg);
            
            % regress DY over df to remove the trend
            df = diff(f1,1);
            b1 = diff(Y, 1,2)*df'/(df*df');
            
            Yres = Y-b1*f1;
            b0 = median(Yres, 2);
            %             b0 = quantile(Yres, 0.05,2);
            Y = bsxfun(@minus, Yres, b0);
            bg.b = [b1, b0];
            bg.f = [f1; ones(1,T)];
        end
        
        %% select background pixels
        function [f1, ind_bg] = findBG(obj, Y, q)
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            if nargin<3;    q = 0.1; end;   %quantiles for selecting background pixels
            % model the background as linear function
            if ndims(Y)==3; Y = neuron.reshape(Y,1); end
            T = size(Y,2);
            % find background area
            ind_frame = round(linspace(1, T, min(T, 500)));
            tmp_C1 = correlation_image(Y(:, ind_frame), [1, 2], d1, d2);
            tmp_C2 = correlation_image(Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
            tmp_Cn = tmp_C1(:)-tmp_C2(:);
            ind = (tmp_Cn>0);
            ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), q));
            % get the mean activity within the selected area
            Ybg = mean(Y(ind_bg(:), :), 1);
            f1 = (Ybg-mean(Ybg))/std(Ybg);
        end
        %% play movie
        function playMovie(obj, Y, min_max, t_pause)
            % play movies
            if ismatrix(Y); Y=obj.reshape(Y, 2); end
            [~, ~, T] = size(Y);
            if (nargin<3) || (isempty(min_max));
                temp = Y(:, :, randi(T, min(100, T), 1));
                min_max = quantile(temp(:), [0.2, 0.9999]);
                min_max(1) = max(min_max(1), 0);
            end
            if nargin<4; t_pause=0.01; end
            figure;
            for t=1:size(Y,3)
                imagesc(Y(:, :, t), min_max);
                axis equal; axis off;
                title(sprintf('Frame %d', t));
                pause(t_pause);
            end
        end
        
        %% trim spatial components
        function [ind_nonzero] = trimSpatial(obj, ratio)
            % remove small nonzero pixels
            if nargin<2;    ratio = 50; end;
            tmp_A = obj.A;
            ind_nonzero = bsxfun(@gt, tmp_A, max(tmp_A, [], 1)/ratio);
            obj.A(~ind_nonzero) = 0;
        end
        
        %% solve A & C with regression
        function [ind_delete] = regressAC(obj, Y)
            if ~ismatrix(Y); Y=obj.reshape(Y,2); end
            maxIter = 1;
            A2 = obj.A;  % initial spatial components
            for miter=1:maxIter
                ind_nonzero = obj.trimSpatial(50); % remove tiny nonzero pixels
                active_pixel = determine_search_location(A2, 'dilate', obj.options);
                K = size(A2, 2); %number of neurons
                
                tmp_C = (A2'*A2)\(A2'*Y);        % update temporal components
                tmp_C(tmp_C<0) = 0;
                % %                 tmp_C(bsxfun(@gt, quantile(tmp_C, 0.95, 2), tmp_C)) = 0;
                %                 tmp_C(bsxfun(@gt, mean(tmp_C, 2)+std(tmp_C, 0.0, 2), tmp_C)) = 0;
                A2 = (Y*tmp_C')/(tmp_C*tmp_C');           % update spatial components
                A2(or(A2<0, ~active_pixel)) = 0;
                
                for m=1:K
                    % remove disconnected pixels
                    img = obj.reshape(A2(:,m), 2);
                    lb = bwlabel(img>max(img(:))/50, 4);
                    img(lb~=mode(lb(ind_nonzero(:, m)))) = 0 ; %find pixels connected to the original components
                    A2(:,m) = img(:);
                end
                
                center_ratio = sum(A2.*ind_nonzero, 1)./sum(A2, 1); %
                % captured too much pixels, ignore newly captured pixels
                ind_too_much = (center_ratio<0.4);
                A2(:, ind_too_much) = A2(:, ind_too_much) .*ind_nonzero(:, ind_too_much);
                % captured too much noiser pixels or all pxiels are zero, ignore this neuron
                ind_delete = or(center_ratio<0.01, sum(A2,1)==0);
                A2(:, ind_delete) = [];
                obj.A = A2;
            end
            temp = (A2'*A2)\(A2'*Y);
            temp(temp<0) = 0;
            obj.C = temp;
        end
        
        %% regress A given C
        function [A2, ind_neuron] = regressA(obj, Y, C)
            if ~ismatrix(Y); Y=obj.reshape(Y,2); end
            A2 = obj.A;  % initial spatial components
            ind_nonzero = obj.trimSpatial(50); % remove tiny nonzero pixels
            active_pixel = determine_search_location(A2, 'dilate', obj.options);
            K = size(A2, 2); %number of neurons
            
            A2 = (Y*C')/(C*C');           % update spatial components
            A2(or(A2<0, ~active_pixel)) = 0;
            
            for m=1:K
                % remove disconnected pixels
                img = obj.reshape(A2(:,m), 2);
                lb = bwlabel(img>max(img(:))/50, 4);
                img(lb~=mode(lb(ind_nonzero(:, m)))) = 0 ; %find pixels connected to the original components
                A2(:,m) = img(:);
            end
            
            center_ratio = sum(A2.*ind_nonzero, 1)./sum(A2, 1); %
            % captured too much pixels, ignore newly captured pixels
            ind_too_much = (center_ratio<0.3);
            A2(:, ind_too_much) = A2(:, ind_too_much) .*ind_nonzero(:, ind_too_much);
            % captured too much noiser pixels or all pxiels are zero, ignore this neuron
            ind_delete = or(center_ratio<0., sum(A2,1)==0);
            A2(:, ind_delete) = [];
            %             C(ind_delete) = [];
            ind_neuron = (1:K);
            ind_neuron(ind_delete) =[];
        end
        
        %% view results
        function runMovie(obj, Y, min_max, save_avi, avi_name, S)
            ctr = obj.estCenter();
            if ~exist('save_avi', 'var')||isempty(save_avi); save_avi=false; end
            if ~exist('avi_name', 'var'); avi_name = []; end
            if ~exist('S', 'var');  S = []; end
            run_movie(Y, obj.A, obj.C, obj.Cn, min_max, obj.Coor, ctr, 5, 1, save_avi, avi_name, S)
        end
        
        %% function
        function image(obj, a, min_max)
            if isvector(a); a = obj.reshape(a,2); end
            if nargin<3; imagesc(a); else imagesc(a, min_max); end
        end
        
        %% normalize
        function normalize(obj)
            norm_A = max(obj.A, [], 1);
            tmp_A = bsxfun(@times, obj.A, 1./norm_A);
            tmp_C = bsxfun(@times, obj.C, norm_A');
            obj.A = tmp_A; obj.C = tmp_C;
        end
        
        %% load data
        function [Y, neuron] = load_data(obj, nam, sframe, num2read)
            ssub = obj.options.ssub;    % spatial downsampling factor
            tsub = obj.options.tsub;    % temporal downsampling factor
            d1 = obj.options.d1;        % image height
            d2 = obj.options.d2;        % image width
            Tbatch = round((2^28)/(d1*d2)/tsub)*tsub; % maximum memory usage is 2GB
            [~, ~, file_type] = fileparts(nam);
            if strcmpi(file_type, '.mat')
                data = matfile(nam);
                Ysiz = data.Ysiz;
                numFrame = Ysiz(3);
                img = data.Y(:, :, 1);
            elseif strcmpi(file_type, '.tif') || strcmpi(file_type, '.tiff')
                numFrame = length(imfinfo(nam));
                img = imread(nam);
            end
            num2read = min(num2read, numFrame-sframe+1); % frames to read
            
            if Tbatch>=num2read
                % load all data because the file is too small
                if strcmpi(file_type, '.mat')
                    Yraw = data.Y;
                elseif strcmpi(file_type, '.tif') || strcmpi(file_type, '.tiff')
                    Yraw = bigread2(nam, sframe, num2read);
                else
                    fprintf('\nThe input file format is not supported yet\n\n');
                    return;
                end
                [neuron, Y] = obj.downSample(double(Yraw));
            else
                % load data in batch model
                [d1s, d2s] = size(imresize(img, 1/ssub));  % size of spatial downsampled data
                Ts = floor(num2read/tsub);  % frame number after downsampling
                Y = zeros(d1s, d2s, Ts);    % downsampled data size
                lastframe = sframe + Ts*tsub -1;  % index of the last frame for loading
                frame0 = sframe;
                while sframe<= lastframe
                    tmp_num2read =  min(lastframe-sframe+1, Tbatch);
                    if strcmpi(file_type, '.mat')
                        fprintf('load data from frame %d to frame %d of %d total frames\n', sframe, sframe+tmp_num2read-1, lastframe);
                        Yraw = data.Y(:, :, sframe:(sframe+tmp_num2read-1));
                    elseif strcmpi(file_type, '.tif') || strcmpi(file_type, '.tiff')
                        Yraw = bigread2(nam, sframe, tmp_num2read);
                    else
                        fprintf('\nThe input file format is not supported yet\n\n');
                        return;
                    end
                    [neuron, temp] = obj.downSample(double(Yraw));
                    Y(:, :, (sframe-frame0)/tsub + (1:size(temp, 3))) = temp;
                    sframe = sframe + size(Yraw,3);
                end
            end
        end
        
        %% merge neurons
        function img = overlapA(obj, ind, ratio)
            %merge all neurons' spatial components into one singal image
            if nargin<2 || isempty(ind)
                AA = obj.A;
            else
                AA = obj.A(:, ind);
            end
            if nargin<3
                ratio = 0.1;
            end
            AA = bsxfun(@times, AA, 1./sum(AA,1));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            [d, K] = size(AA);
            
            col = randi(6, 1, K);
            img = zeros(d, 3);
            for m=1:3
                img(:, m) = sum(bsxfun(@times, AA, mod(col, 2)), 2);
                col = round(col/2);
            end
            img = obj.reshape(img, 2);
            img = img/max(img(:))*(2^16);
            img = uint16(img);
        end
        
        %% play video
        function playAC(obj, avi_file, cell_id)
            if nargin<3
                cell_id = 1:size(obj.C, 1);
            end
            [K, T] = size(obj.C(cell_id, :));
            % draw random color for each neuron
            tmp = mod((1:K)', 6)+1;
            col = zeros(K, 3);
            for m=1:3
                col(:, m) = mod(tmp, 2);
                tmp = round(tmp/2);
            end
            figure;
            % play
            if nargin>1
                fp = VideoWriter(avi_file);
                fp.open();
            end
            
            cmax = max(reshape(obj.A*obj.C(:, 1:100:end), 1, []));
            for m=1:T
                img = obj.A(:, cell_id)*bsxfun(@times, obj.C(cell_id,m), col);
                img = obj.reshape(img, 2)/cmax*500;
                imagesc(uint8(img));
                axis equal off tight;
                title(sprintf('Time %.2f seconds', m/obj.Fs));
                
                pause(.1);
                if nargin>1
                    frame = getframe(gcf);
                    fp.writeVideo(frame);
                end
            end
            if nargin>1; fp.close(); end
        end
        
        %% find neurons from the residual
        % you can do it in either manual or automatic way
        function [center, Cn, pnr] = pickNeurons(obj, Y, patch_par, seed_method)
            if ~exist('patch_par', 'var')||isempty(patch_par)
                seed_method = [3,3];
            end
            if ~exist('seed_method', 'var')||isempty(seed_method)
                seed_method = 'auto';
            end
            neuron = obj.copy();
            neuron.options.seed_method = seed_method;
            [center, Cn, pnr] = neuron.initComponents_endoscope(Y, [], patch_par, false, false);
            obj.A = [obj.A, neuron.A];
            obj.C = [obj.C; neuron.C];
        end
        
        %% post process spatial component
        function A_ = post_process_spatial(obj, A_)
            if ~exist('A_', 'var');
                A_ = obj.A;
            end
            A_ = threshold_components(A_, obj.options);
            obj.A = A_;
        end
        
        %% estimate local background
        function Ybg = localBG(obj, Ybg, ssub, rr, IND, method)
            if ~exist('rr', 'var')||isempty(rr); rr=obj.options.gSiz; end
            if ~exist('ssub', 'var')||isempty(ssub); ssub = 1; end
            if ~exist('IND', 'var'); IND = []; end
            if ~exist('method', 'var'); method = 'regression'; end
            Ybg = lle(obj.reshape(Ybg, 2), ssub, rr, IND, method);
            Ybg = obj.reshape(Ybg, 1);
        end
        
        %% save results
        function save_results(obj, file_nm, Ybg) %#ok<INUSD>
            warning('off', 'all');
            neuron_results = struct(obj);  %#ok<NASGU>
            if exist('Ybg', 'var')
                save(file_nm, 'neuron_results', 'Ybg');
            else
                save(file_nm, 'neuron_results');
            end
            warning('on', 'all');
            fprintf('results has been saved into file %s\n', file_nm);
        end
        
        %% event detection
        function E = event_detection(obj, sig, w)
            % detect events by thresholding S with sig*noise
            % can get at most one spike
            % sig: threshold of the minimum amplitude of the events
            
            if ~exist('sig', 'var')|| isempty(sig)
                sig=5;
            end
            
            if ~exist('w', 'var')||isempty(w)
                w = obj.Fs;
            end
            E =obj.C;    % event detection
            Emin = ordfilt2(E, 1, ones(1, w));
            Emax = ordfilt2(E, w, ones(1, w));
            E(E~=Emax) = 0;  % only select local maximums
            for m=1:size(E,1)
                E(m, E(m, :)-Emin(m, :)< obj.P.neuron_sn{m}*sig) = 0; % remove small transients
            end
        end
        
        % compute correlation image and peak to noise ratio for endoscopic
        % data. unlike the correlation image for two-photon data,
        % correlation image of the microendoscopic data needs to be
        % spatially filtered first. otherwise neurons are significantly
        % overlapped.
        function [Cn, PNR] = correlation_pnr(obj, Y)
            [Cn, PNR] = correlation_image_endoscope(Y, obj.options);
            %             obj.Cn = Cn;
            %             obj.PNR = PNR;
        end
    end
    
end
