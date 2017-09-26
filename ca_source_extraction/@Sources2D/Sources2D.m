classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data. It
    % Both CNMF and CNMF-E model can be used.
    
    % Author: Pengcheng Zhou, Columbia University, 2017
    % zhoupc1988@gmail.com
    
    %% properties
    properties
        % spatial
        A;          % spatial components of neurons
        % temporal
        C;          % temporal components of neurons
        C_raw;      % raw traces of temporal components
        S;          % spike counts
        kernel;     % calcium dynamics. this one is less used these days.
        % background
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        W;          % a sparse weight matrix matrix
        b0;         % constant baselines for each pixel
        % optiosn
        options;    % options for model fitting
        P;          % some estimated parameters or parameters relating to data
        % data info
        Fs = nan;    % frame rate
        file = '';
        frame_range;  % frame range of the data
        %quality control
        ids;        % unique identifier for each neuron
        tags;       % tags bad neurons with multiple criterion using a 16 bits number
        % ai indicates the bit value of the ith digits from
        % right to the left
        % a1 = 0 means that neuron has too few nonzero pixels
        % a2 = 0 means that neuron has silent calcium transients
        % a3 = 0 means that neuron has unrealistic neuron
        % shapes
        % a4 = 0 indicates that user doesn't want this neuron
        % a5 = 0 indicates that the neuron has been merged with
        % others
        %others
        Cn;         % correlation image
        PNR;        % peak-to-noise ratio image.
        Coor;       % neuron contours
        Df;         % background for each component to normalize the filtered raw data
        C_df;       % temporal components of neurons and background normalized by Df
        S_df;       % spike counts of neurons normalized by Df
        
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = Sources2D(varargin)
            obj.options = CNMFSetParms();
            obj.P =struct('mat_file', [], 'mat_data', [], 'indicator', '', 'k_options', 0, ...
                'k_snapshot', 0, 'k_del', 0, 'k_merge', 0, 'k_trim', 0, 'sn', [], ...
                'kernel_pars',[], 'k_ids', 0);
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:});
            end
            obj.kernel = create_kernel('exp2');
        end
        
        %------------------------------------------------------------------OPTIONS------------%
        %% update parameters
        function updateParams(obj, varargin)
            obj.options = CNMFSetParms(obj.options, varargin{:});
        end
        
        %------------------------------------------------------------------PREPROCESSING-------%
        %% data preprocessing
        function Y = preprocess(obj,Y,p)
            [obj.P,Y] = preprocess_data(Y,p,obj.options);
        end
        
        %------------------------------------------------------------------LOAD DATA-----------%
        %% select data
        function nam = select_data(obj, nam)
            %% select file
            if ~exist('nam', 'var') || isempty(nam) || ~exist(nam, 'file')
                % choose files using graphical interfaces
                try
                    load .dir.mat;          %load previous path
                catch
                    dir_nm = [cd(), filesep];   %use the current path
                end
                [file_nm, dir_nm] = uigetfile(fullfile(dir_nm, '*.tif;*.mat;*.h5;*.avi'));
                nam = [dir_nm, file_nm];        % full name of the data file
            else
                % file has been selected
                nam =  get_fullname(nam);
                [dir_nm, ~, ~] = fileparts(nam);
            end
            
            
            if dir_nm~=0
                save .dir.mat dir_nm;
            else
                fprintf('no file was selected. STOP!\n');
                return;
            end
            obj.file = nam;
        end
        
        function Ypatch = load_patch_data(obj, patch_pos, frame_range)
            %% load data within selected patch position and selected frames
            mat_data = obj.P.mat_data;
            
            dims = mat_data.dims;
            d1 = dims(1);
            d2 = dims(2);
            T = dims(3);
            if ~exist('patch_pos', 'var') || isempty(patch_pos)
                patch_pos = [1, d1, 1, d2];
            end
            if ~exist('frame_range', 'var') || isempty(frame_range)
                frame_range = [1, T];
            end
            Ypatch = get_patch_data(mat_data, patch_pos, frame_range, false);
        end
        %% distribute data and be ready to run source extraction
        function getReady(obj, pars_envs)
            % data and its extension
            nam = obj.file;
            
            % parameters for scaling things
            if ~exist('pars_envs', 'var') || isempty(pars_envs)
                memory_size_to_use = 16.0;  %GB
                memory_size_per_patch = 1.0;  % GB;
                patch_dims = [64, 64];
            else
                memory_size_to_use = pars_envs.memory_size_to_use;
                memory_size_per_patch = pars_envs.memory_size_per_patch;
                patch_dims = pars_envs.patch_dims;
            end
            
            % overlapping area
            w_overlap = obj.options.ring_radius;
            
            % distribute data
            [data, dims, obj.P.folder_analysis] = distribute_data(nam, patch_dims, w_overlap, memory_size_per_patch, memory_size_to_use);
            obj.P.mat_data = data;
            obj.P.mat_file = data.Properties.Source;
            obj.updateParams('d1', dims(1), 'd2', dims(2));
            obj.P.numFrames = dims(3);
            
            % estimate the noise level for each pixel
            try
                obj.P.sn = obj.P.mat_data.sn;
            catch
                obj.P.sn = obj.estimate_noise();
                mat_data = obj.P.mat_data;
                mat_data.Properties.Writable = true;
                mat_data.sn = obj.P.sn;
                mat_data.Properties.Writable = false;
            end
        end
        
        %% estimate noise
        function sn = estimate_noise(obj, frame_range, method)
            mat_data = obj.P.mat_data;
            dims = mat_data.dims;
            T = dims(3);
            if ~exist('frame_range', 'var') || isempty(frame_range)
                frame_range = [1, T];
            end
            if ~exist('method', 'var') || isempty(method)
                method = 'psd';
            end
            if strcmpi(method, 'hist')
                % fit the histogram of with a parametric gaussian shape.
                % it works badly when the
                % most frames are involved in calcium transients.
                foo = @(Y) GetSn_hist(Y, false);
            elseif strcmpi(method, 'std')
                % simly use the standard deviation. works badly when
                % neurons have large calcium transients
                foo = @(Y) std(Y, 0, 2);
            else
                % default options is using the power spectrum method.
                % works badly when noises have temporal correlations
                foo = @GetSn;
            end
            block_idx_r = mat_data.block_idx_r;
            block_idx_c = mat_data.block_idx_c;
            nr_block = length(block_idx_r) - 1;
            nc_block = length(block_idx_c) - 1;
            sn = cell(nr_block, nc_block);
            fprintf('estimating the nosie level for every pixel.......\n');
            tic;
            parfor mblock=1:(nr_block*nc_block)
                [m, n] = ind2sub([nr_block, nc_block], mblock);
                r0 = block_idx_r(m); %#ok<PFBNS>
                r1 = block_idx_r(m+1);
                c0 = block_idx_c(n); %#ok<PFBNS>
                c1 = block_idx_c(n+1);
                Ypatch = get_patch_data(mat_data, [r0, r1, c0, c1], frame_range, false);
                Ypatch = double(reshape(Ypatch, [], T));
                tmp_sn = reshape(foo(Ypatch), r1-r0+1, c1-c0+1);
                if m~=nr_block
                    tmp_sn(end-1, :) = [];
                end
                if n~=nc_block
                    tmp_sn(:, end-1) = [];
                end
                sn{mblock} = tmp_sn;
            end
            fprintf('Time cost for estimating the nosie levels:  %.3f \n\n', toc);
            sn = cell2mat(sn);
        end
        
        %% pick parameters
        [gSig, gSiz, ring_radius, min_corr, min_pnr] = set_parameters(obj)
        %------------------------------------------------------------------INITIALIZATION-----------%
        %% fast initialization, for 2P data, CNMF
        function [center] = initComponents(obj, Y, K, tau)
            if nargin<4 ;    tau = [];             end
            [obj.A, obj.C, obj.b, obj.f, center] = initialize_components(Y, K, tau, obj.options);
        end
        
        %% initialize neurons using patch method
        % for 1P and 2P data, CNMF and CNMF-E
        [center, Cn, PNR] = initComponents_parallel(obj, K, frame_range, save_avi, use_parallel)
        
        %% fast initialization for microendoscopic data
        [center, Cn, pnr] = initComponents_endoscope(obj, Y, K, patch_sz, debug_on, save_avi);
        
        [center] = initComponents_2p(obj,Y, K, options, sn, debug_on, save_avi);
        
        %% pick neurons from the residual
        % for 1P and 2P data, CNMF and CNMF-E
        [center, Cn, PNR] = initComponents_residual_parallel(obj, K, save_avi, use_parallel)
        
        %------------------------------------------------------------------UPDATE MODEL VARIABLES---%
        %% update background components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_background_parallel(obj, use_parallel)
        
        %% update spatial components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_spatial_parallel(obj, use_parallel)
        
        %% update spatial components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_temporal_parallel(obj, use_parallel)
        
        %% update spatial components, 2P data, CNMF
        function updateSpatial(obj, Y)
            [obj.A, obj.b, obj.C] = update_spatial_components(Y, ...
                obj.C, obj.f, obj.A, obj.P, obj.options);
        end
        
        %% udpate spatial components without background, no background components
        function updateSpatial_nb(obj, Y)
            [obj.A, obj.C] = update_spatial_components_nb(Y, ...
                obj.C, obj.A, obj.P, obj.options);
        end
        
        %% udpate spatial components using NNLS and thresholding
        function updateSpatial_nnls(obj, Y)
            [obj.A, obj.C] = update_spatial_nnls(Y, ...
                obj.C, obj.A, obj.P, obj.options);
        end
        
        %% update temporal components for endoscope data
        updateSpatial_endoscope(obj, Y, numIter, method, IND_thresh);
        
        %% update temporal components
        function updateTemporal(obj, Y)
            [obj.C, obj.f, obj.P, obj.S] = update_temporal_components(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        %% udpate temporal components with fast deconvolution
        updateTemporal_endoscope(obj, Y, allow_deletion)
        updateTemporal_endoscope_parallel(obj, Y, smin)
        
        %% update temporal components without background
        function updateTemporal_nb(obj, Y)
            [obj.C, obj.P, obj.S] = update_temporal_components_nb(...
                Y, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        %------------------------------------------------------------------MERGE NEURONS ---------%
        %% merge found components
        function [nr, merged_ROIs] = merge(obj, Y)
            [obj.A, obj.C, nr, merged_ROIs, obj.P, obj.S] = merge_components(...
                Y,obj.A, [], obj.C, [], obj.P,obj.S, obj.options);
        end
        
        %% quick merge neurons based on spatial and temporal correlation
        [merged_ROIs, newIDs] = quickMerge(obj, merge_thr)
        
        %% quick merge neurons based on spatial and temporal correlation
        [merged_ROIs, newIDs] = MergeNeighbors(obj, dmin, method)
        
        %------------------------------------------------------------------EXPORT---%
        function save_neurons(obj)
            folder_nm = [obj.P.log_folder, 'neurons'];
            obj.viewNeurons([], obj.C_raw, folder_nm);
        end
        
        show_demixed_video(obj,save_avi, kt, center_ac, range_ac, range_Y, multi_factor)        %% compute the residual
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
        
        %% extract DF/F signal for microendoscopic data
        function [C_df,C_raw_df, Df] = extract_DF_F_endoscope(obj, Ybg)
            options_ = obj.options;
            A_ = obj.A;
            C_ = obj.C;
            [K,T] = size(obj.C);  % number of frames
            Ybg = bsxfun(@times, A_, 1./sum(A_.^2, 1))' * Ybg;   % estimate the background for each neurons
            if isempty(options_.df_window) || (options_.df_window > T)
                % use quantiles of the whole recording session as the
                % baseline
                if options_.df_prctile == 50
                    Df = median(Ybg,2);
                else
                    Df = prctile(Ybg,options_.df_prctile,2);
                end
                C_df = spdiags(Df,0,K,K)\C_;
                C_raw_df = spdiags(Df,0,K,K)\obj.C_raw;
            else
                % estimate the baseline for each short period
                if options_.df_prctile == 50
                    Df = medfilt1(Ybg,options_.df_window,[],2,'truncate');
                else
                    Df = zeros(size(Ybg));
                    for i = 1:size(Df,1)
                        df_temp = running_percentile(Ybg(i,:), options_.df_window, options_.df_prctile);
                        Df(i,:) = df_temp(:)';
                    end
                end
                C_df = C_./Df;
                C_raw_df = obj.C_raw./Df;
            end
        end
        
        %% order_ROIs
        function [srt] = orderROIs(obj, srt)
            %% order neurons
            % srt: sorting order
            nA = sqrt(sum(obj.A.^2));
            nr = length(nA);
            if nargin<2
                srt='srt';
            end
            K = size(obj.C, 1);
            
            if ischar(srt)
                if strcmpi(srt, 'decay_time')
                    % time constant
                    if K<0
                        disp('Are you kidding? You extracted 0 neurons!');
                        return;
                    else
                        taud = zeros(K, 1);
                        for m=1:K
                            temp = ar2exp(obj.P.kernel_pars(m));
                            taud(m) = temp(1);
                        end
                        [~, srt] = sort(taud);
                    end
                elseif strcmp(srt, 'mean')
                    [~, srt] = sort(mean(obj.C,2)'.*sum(obj.A), 'descend');
                elseif strcmp(srt, 'circularity')
                    % order neurons based on its circularity
                    tmp_circularity = zeros(K,1);
                    for m=1:K
                        [w, r] = nnmf(obj.reshape(obj.A(:, m),2), 1);
                        ky = sum(w>max(w)*0.3);
                        kx = sum(r>max(r)*0.3);
                        tmp_circularity(m) = abs((kx-ky+0.5)/((kx+ky)^2));
                    end
                    [~, srt] = sort(tmp_circularity, 'ascend');
                elseif strcmpi(srt, 'snr')
                    snrs = var(obj.C, 0, 2)./var(obj.C_raw-obj.C, 0, 2);
                    [~, srt] = sort(snrs, 'descend');
                end
            end
            obj.A = obj.A(:, srt);
            obj.C = obj.C(srt, :);
            obj.ids = obj.ids(srt);
            obj.tags = obj.tags(srt);
            
            try
                obj.C_raw = obj.C_raw(srt,:);
                obj.S = obj.S(srt,:);
                obj.P.kernel_pars = obj.P.kernel_pars(srt, :);
                obj.P.neuron_sn = obj.P.neuron_sn(srt);
            end
        end
        
        %% view contours
        function [Coor, json_file] = viewContours(obj, Cn, contour_threshold, display, ind, ln_wd)
            if or(isempty(Cn), ~exist('Cn', 'var') )
                Cn = reshape(obj.P.sn, obj.options.d1, obj.options.d2);
            end
            if (~exist('display', 'var') || isempty(display)); display=0; end
            if (~exist('ind', 'var') || isempty(ind)); ind=1:size(obj.A, 2); end
            if (~exist('ln_wd', 'var') || isempty(ln_wd)); ln_wd = 1; end
            [obj.Coor, json_file] = plot_contours(obj.A(:, ind), Cn, ...
                contour_threshold, display, [], [], ln_wd);
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
        
        %% quick view
        viewNeurons(obj, ind, C2, folder_nm);
        displayNeurons(obj, ind, C2, folder_nm);
        
        %% delete neurons
        function delete(obj, ind)
            % write the deletion into the log file
            if ~exist('ind', 'var') || isempty(ind)
                return;
            end
            try
                % folders and files for saving the results
                log_file =  obj.P.log_file;
                flog = fopen(log_file, 'a');
                log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
                ids_del = obj.ids(ind);
                if isempty(ids_del)
                    fprintf('no neurons to be deleted\n');
                    fclose(flog);
                    return;
                end
                tmp_str = get_date();
                tmp_str=strrep(tmp_str, '-', '_');
                eval(sprintf('log_data.ids_del_%s = ids_del;', tmp_str));
                
                fprintf(flog, '[%s]\b', get_minute());
                fprintf('Deleted %d neurons: ', length(ids_del));
                fprintf(flog, 'Deleted %d neurons: \n', length(ids_del));
                for m=1:length(ids_del)
                    fprintf('%2d, ', ids_del(m));
                    fprintf(flog, '%2d, ', ids_del(m));
                end
                fprintf(flog, '\n');
                fprintf('\nThe IDS of these neurons were saved as intermediate_results.spatial_%s\n\n', tmp_str);
                fprintf(flog, '\tThe IDS of these neurons were saved as intermediate_results.spatial_%s\n\n', tmp_str);
                fclose(flog);
            end
            
            obj.A(:, ind) = [];
            obj.C(ind, :) = [];
            if ~isempty(obj.S)
                try obj.S(ind, :) = []; catch; end
            end
            if ~isempty(obj.C_raw)
                try obj.C_raw(ind, :) = []; catch;  end
            end
            if isfield(obj.P, 'kernel_pars')&&(  ~isempty(obj.P.kernel_pars))
                try obj.P.kernel_pars(ind, :) = []; catch; end
            end
            try  obj.ids(ind) = [];   catch;   end
            try obj.tags(ind) =[]; catch; end
            
            % save the log
        end
        
        %% merge neurons
        
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
            if dim==1
                Y=reshape(Y, d1*d2, []);  %each frame is a vector
            else
                Y = reshape(full(Y), d1, d2, []);    %each frame is an image
            end
        end
        
        %% deconvolve all temporal components
        C_ = deconvTemporal(obj, use_parallel)
        
        %% play movie
        function playMovie(obj, Y, min_max, col_map, avi_nm, t_pause)
            Y = double(Y);
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            % play movies
            figure('papersize', [d2,d1]/max(d1,d2)*5);
            width = d2/max(d1,d2)*500;
            height =d1/max(d1,d2)*500;
            set(gcf, 'position', [500, 200, width, height]);
            axes('position', [0,0, 1, 1]);
            if ~exist('col_map', 'var') || isempty(col_map)
                col_map = jet;
            end
            if exist('avi_nm', 'var') && ischar(avi_nm)
                avi_file = VideoWriter(avi_nm);
                if ~isnan(obj.Fs)
                    avi_file.FrameRate = obj.Fs;
                end
                avi_file.open();
                avi_flag = true;
            else
                avi_flag = false;
            end
            if ismatrix(Y); Y=obj.reshape(Y, 2); end
            [~, ~, T] = size(Y);
            if (nargin<3) || (isempty(min_max))
                temp = Y(:, :, randi(T, min(100, T), 1));
                min_max = quantile(temp(:), [0.2, 0.9999]);
                min_max(1) = max(min_max(1), 0);
                min_max(2) = max(min_max(2), min_max(1)+0.1);
            end
            if ~exist('t_pause', 'var'); t_pause=0.01; end
            for t=1:size(Y,3)
                imagesc(Y(:, :, t), min_max); colormap(col_map);
                axis equal; axis off tight;
                if isnan(obj.Fs)
                    title(sprintf('Frame %d', t));
                else
                    text(1, 10, sprintf('Time = %.2f', t/obj.Fs), 'fontsize', 15, 'color', 'w');
                end
                pause(t_pause);
                if avi_flag
                    temp = getframe(gcf);
                    temp.cdata = imresize(temp.cdata, [width, height]);
                    avi_file.writeVideo(temp);
                end
            end
            if avi_flag
                avi_file.close();
            end
        end
        
        %% export AVI
        function exportAVI(obj, Y, min_max, avi_nm, col_map)
            % export matrix data to movie
            %min_max: 1*2 vector, scale
            %avi_nm: string, file name
            %col_map: colormap
            
            T = size(Y, ndims(Y));
            Y = Y(:);
            if ~exist('col_map', 'var') || isempty(col_map)
                col_map = jet;
            end
            if ~exist('avi_nm', 'var') || isempty(avi_nm)
                avi_nm = 'a_movie_with_no_name.avi';
            end
            if ~exist('min_max', 'var') || isempty(min_max)
                min_max = [min(Y(1:10:end)), max(Y(1:10:end))];
            end
            
            Y = uint8(64*(Y-min_max(1))/diff(min_max));
            Y(Y<1) = 1;
            Y(Y>64) = 64;
            col_map = uint8(255*col_map);
            Y = reshape(col_map(Y, :), obj.options.d1, [], T, 3);
            Y = permute(Y, [1,2,4,3]);
            
            avi_file = VideoWriter(avi_nm);
            avi_file.open();
            for m=1:T
                %                 temp.cdata = squeeze(Y(:, :, :, m));
                %                 temp.colormap = [];
                avi_file.writeVideo(squeeze(Y(:, :, :, m)));
            end
            avi_file.close();
        end
        
        %% trim spatial components
        function [ind_small] = trimSpatial(obj, thr, sz)
            % remove small nonzero pixels
            if nargin<2;    thr = 0.01; end
            if nargin<3;    sz = 5; end
            
            se = strel('square', sz);
            ind_small = false(size(obj.A, 2), 1);
            for m=1:size(obj.A,2)
                ai = obj.A(:,m);
                ai_open = imopen(obj.reshape(ai,2), se);
                
                temp = full(ai_open>max(ai)*thr);
                l = bwlabel(obj.reshape(temp,2), 4);   % remove disconnected components
                [~, ind_max] = max(ai_open(:));
                
                ai(l(:)~=l(ind_max)) = 0;
                if sum(ai(:)>0) < obj.options.min_pixel %the ROI is too small
                    ind_small(m) = true;
                end
                obj.A(:, m) = ai(:);
            end
            ind_small = find(ind_small);
            obj.delete(ind_small);
        end
        
        %% keep spatial shapes compact
        function compactSpatial(obj)
            ind_del = false(1, size(obj.A,2));
            for m=1:size(obj.A, 2)
                ai = obj.reshape(obj.A(:, m), 2);
                ai = spatial_constraints(ai);
                if sum(ai(:))>=obj.options.min_pixel
                    obj.A(:, m) = ai(:);
                else
                    ind_del(m) = true;
                end
            end
            obj.delete(ind_del);
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
            run_movie(Y, obj.A, obj.C_raw, obj.Cn, min_max, obj.Coor, ctr, 5, 1, save_avi, avi_name, S)
        end
        
        %% function
        function image(obj, a, min_max)
            if isvector(a); a = obj.reshape(a,2); end
            if nargin<3
                imagesc(a);
            else
                imagesc(a, min_max);
            end
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
            elseif strcmpi(file_type, '.h5') || strcmpi(file_type, '.hdf5')
                temp = h5info(nam);
                dataset_nam = ['/', temp.Datasets.Name];
                dataset_info = h5info(nam, dataset_nam);
                dims = dataset_info.Dataspace.Size;
                ndims = length(dims);
                numFrame = dims(end);
                img = squeeze(h5read(nam, dataset_nam, ones(1, ndims), [1,d1, d2, 1, 1]));
            end
            num2read = min(num2read, numFrame-sframe+1); % frames to read
            
            if Tbatch>=num2read
                % load all data because the file is too small
                if strcmpi(file_type, '.mat')
                    Yraw = data.Y(:, :, (1:num2read)+sframe-1);
                elseif strcmpi(file_type, '.tif') || strcmpi(file_type, '.tiff')
                    Yraw = bigread2(nam, sframe, num2read);
                elseif strcmpi(file_type, '.h5') || strcmpi(file_type, 'hdf5')
                    Yraw = squeeze(h5read(nam, dataset_nam, ones(1, ndims), [1, d1, d2, 1, num2read]));
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
                        fprintf('load data from frame %d to frame %d of %d total frames\n', sframe, sframe+tmp_num2read-1, lastframe);
                        Yraw = bigread2(nam, sframe, tmp_num2read);
                    elseif strcmpi(file_type, '.h5') || strcmpi(file_type, 'hdf5')
                        fprintf('load data from frame %d to frame %d of %d total frames\n', sframe, sframe+tmp_num2read-1, lastframe);
                        Yraw = squeeze(h5read(nam, dataset_nam, [1,1,1,1,sframe], [1, d1, d2, 1, tmp_num2read]));
                    else
                        fprintf('\nThe input file format is not supported yet\n\n');
                        return;
                    end
                    [neuron, temp] = obj.downSample(double(Yraw));
                    Y(:, :, (sframe-frame0)/tsub + (1:size(temp, 3))) = temp;
                    sframe = sframe + size(Yraw,3);
                end
            end
            neuron.options.min_pixel = ceil(obj.options.min_pixel/(ssub^2));
        end
        
        %% estimate noise
        function sn = estNoise(obj, Y)
            fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
            Y = obj.reshape(Y, 1);
            sn = get_noise_fft(Y,obj.options);
            obj.P.sn = sn(:);
            fprintf('  done \n');
        end
        
        %% reconstruct background
        function Ybg = reconstruct_background(obj)
            %%reconstruct background using the saved data
            % input:
            %   use_parallel: boolean, do initialization in patch mode or not.
            %       default(true); we recommend you to set it false only when you want to debug the code.
            
            %% Author: Pengcheng Zhou, Columbia University, 2017
            %% email: zhoupc1988@gmail.com
            
            %% process parameters
            
            try
                % map data
                mat_data = obj.P.mat_data;
                
                % dimension of data
                dims = mat_data.dims;
                d1 = dims(1);
                d2 = dims(2);
                T = dims(3);
                obj.options.d1 = d1;
                obj.options.d2 = d2;
                
                % parameters for patching information
                patch_pos = mat_data.patch_pos;
                block_pos = mat_data.block_pos;
                
                % number of patches
                [nr_patch, nc_patch] = size(patch_pos);
            catch
                error('No data file selected');
            end
            
            % frames to be loaded for initialization
            T = diff(obj.frame_range) + 1;
            
            % threshold for detecting large residuals
%             thresh_outlier = obj.options.thresh_outlier;
            
            % options
            %% start updating the background
            bg_model = obj.options.background_model;
            Ybg = zeros(d1, d2, T);
            for mpatch=1:(nr_patch*nc_patch)
                tmp_patch = patch_pos{mpatch};
                
                if strcmpi(bg_model, 'ring')
                    W_ring = obj.W{mpatch}; 
                    b0_ring = obj.b0{mpatch}; 
                    % load data 
                    Ypatch = get_patch_data(mat_data, tmp_patch, obj.frame_range, true);
                    Ypatch = reshape(Ypatch, [], T); 
                    tmp_block = block_pos{mpatch};
                    
                    % find the neurons that are within the block
                    mask = zeros(d1, d2);
                    mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
                    ind = (reshape(mask(:), 1, [])* obj.A>0);
                    A_patch = obj.A(logical(mask), ind);
                    C_patch = obj.C(ind, :);
                    
                    % reconstruct background 
                    Ymean = mean(Ypatch,2);
                    Cmean = mean(C_patch , 2);
                    Ypatch = bsxfun(@minus, double(Ypatch), Ymean);
                    C_patch = bsxfun(@minus, C_patch, Cmean);
                    Bf = W_ring*(double(Ypatch) - A_patch*C_patch);
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, Bf, b0_ring), diff(tmp_patch(1:2))+1, [], T); 
                elseif strcmpi(bg_model, 'nmf')
                    b_nmf = obj.b{mpatch};
                    f_nmf = obj.f{mpatch};
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(b_nmf*f_nmf, diff(tmp_patch(1:2))+1, [], T);
                else
                    b_svd = obj.b{mpatch};
                    f_svd = obj.f{mpatch};
                    b0_svd = obj.b0{mpatch};
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, b_svd*f_svd, b0_svd), diff(tmp_patch(1:2))+1, [], T);
                end
                
            end
            
        end
        
        
        
        %% concateneate b, f, b0, W
        %         function Ybg = concatenate_var(obj)
        %             %% concatenate a signle variable for divided version of b, f, b0, W
        %             % input:
        %             %    nam: {'b', 'f', 'b0'}
        %             %% Author: Pengcheng Zhou, Columbia University, 2017
        %             %% email: zhoupc1988@gmail.com
        %
        %             %% process parameters
        %
        %             try
        %                 % map data
        %                 mat_data = obj.P.mat_data;
        %                 mat_file = mat_data.Properties.Source;
        %                 tmp_dir = fileparts(mat_file);
        %
        %                 % dimension of data
        %                 dims = mat_data.dims;
        %                 d1 = dims(1);
        %                 d2 = dims(2);
        %                 T = dims(3);
        %                 obj.options.d1 = d1;
        %                 obj.options.d2 = d2;
        %
        %                 % parameters for patching information
        %                 patch_pos = mat_data.patch_pos;
        %                 block_pos = mat_data.block_pos;
        %
        %                 % number of patches
        %                 [nr_patch, nc_patch] = size(patch_pos);
        %             catch
        %                 error('No data file selected');
        %             end
        %
        %             % frames to be loaded for initialization
        %             frame_range = obj.frame_range;
        %             T = diff(frame_range) + 1;
        %
        %             % threshold for detecting large residuals
        %             thresh_outlier = obj.options.thresh_outlier;
        %
        %             % options
        %             %% start updating the background
        %             bg_model = obj.options.background_model;
        %             Ybg = zeros(d1, d2, T);
        %             for mpatch=1:(nr_patch*nc_patch)
        %
        %
        %             end
        %
        %         end
        
        %% merge neurons
        function [img, col0, AA] = overlapA(obj, ind, ratio)
            %merge all neurons' spatial components into one singal image
            if nargin<2 || isempty(ind)
                AA = obj.A;
            else
                AA = obj.A(:, ind);
            end
            if nargin<3
                ratio = 0.3;
            end
            
            %             v_max = max(max(AA,1));
            v_max = 1;
            AA = bsxfun(@times, AA, v_max./max(AA,[],1));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            [d, K] = size(AA);
            
            if K==2
                col = [4, 2];
            elseif K==3
                col = [4, 2, 1];
            else
                col = randi(6, 1, K);
            end
            col0 = col;
            img = zeros(d, 3);
            for m=3:-1:1
                img(:, m) = sum(bsxfun(@times, AA, mod(col, 2)), 2);
                col = floor(col/2);
            end
            img = obj.reshape(img, 2);
            img = img/max(img(:))*(2^16);
            img = uint16(img);
        end
        
        %% play video
        function playAC(obj, avi_file, cell_id, indt)
            if nargin<3 || isempty(cell_id)
                cell_id = 1:size(obj.C, 1);
            end
            [K, T] = size(obj.C(cell_id, :));
            if ~exist('indt', 'var')
                indt = [1, T];
            end
            % draw random color for each neuron
            tmp = mod((1:K)', 6)+1;
            col = zeros(K, 3);
            for m=1:3
                col(:, m) = mod(tmp, 2);
                tmp = floor(tmp/2);
            end
            figure;
            height = obj.options.d1*512/obj.options.d2;
            set(gcf, 'position', [675, 524, 512, height]);
            axes('position', [0,0,1,1]);
            % play
            if nargin>1
                fp = VideoWriter(avi_file);
                fp.FrameRate = obj.Fs;
                fp.open();
            end
            
            cmax = max(reshape(obj.A*obj.C, 1, []));
            for m=indt(1):indt(2)
                img = obj.A(:, cell_id)*bsxfun(@times, obj.C(cell_id,m), col);
                img = obj.reshape(img, 2)/cmax*1000;
                imagesc(uint8(img));
                axis equal off tight;
                text(5, 10, sprintf('Time: %.2f seconds',(m-indt(1))/obj.Fs), 'color', 'w',...
                    'fontsize', 16);
                pause(.01);
                if exist('fp', 'var')
                    frame = getframe(gcf);
                    frame.cdata = imresize(frame.cdata, [512, height]);
                    fp.writeVideo(frame);
                end
            end
            if nargin>1; fp.close(); end
        end
        
        %% find neurons from the residual
        % you can do it in either manual or automatic way
        function [center, Cn, pnr] = pickNeurons(obj, Y, patch_par, seed_method, debug_on)
            if ~exist('patch_par', 'var')||isempty(patch_par)
                patch_par = [3,3];
            end
            if ~exist('seed_method', 'var')||isempty(seed_method)
                seed_method = 'auto';
            end
            if ~exist('debug_on', 'var')||isempty(debug_on)
                debug_on = false;
            end
            neuron = obj.copy();
            neuron.options.seed_method = seed_method;
            neuron.options.gSig = 1;
            neuron.options.center_psf = 0;
            [center, Cn, pnr] = neuron.initComponents_endoscope(Y, [], patch_par, debug_on, false);
            obj.A = [obj.A, neuron.A];
            obj.C = [obj.C; neuron.C];
            obj.C_raw = [obj.C_raw; neuron.C_raw];
            if obj.options.deconv_flag
                obj.S = [obj.S; neuron.S];
                obj.P.kernel_pars = [obj.P.kernel_pars; neuron.P.kernel_pars];
            end
        end
        
        %% post process spatial component
        A_ = post_process_spatial(obj, A_)
        
        %% merge neruons
        [merged_ROIs, newIDs, obj_bk] = merge_high_corr(obj, show_merge, merge_thr);
        [merged_ROIs, newIDs, obj_bk] = merge_close_neighbors(obj, show_merge, merge_thr);
        [merged_ROIs, newIDs, obj_bk] = merge_neurons_dist_corr(obj, show_merge, merge_thr, dmin, method_dist, max_decay_diff);
        
        %% tag neurons with quality control
        function tags_ = tag_neurons_parallel(obj)
            A_ = obj.A;
            %             C_ = obj.C_;
            S_ = obj.S;
            K = size(A_, 2);
            tags_ = zeros(K, 1, 'like', uint16(0));
            min_pixel = obj.options.min_pixel;
            
            % check the number of nonzero pixels
            nz_pixels = full(sum(A_>0, 1));
            tags_ = tags_ + uint16(nz_pixels'<min_pixel);
            
            % check the number of calcium transients after the first frame
            if obj.options.deconv_flag
                nz_spikes = full(sum(S_(:,2:end), 2));
                tags_ = tags_ + uint16(nz_spikes<1)*2;
            end
            
            obj.tags = tags_;
        end
        %% estimate local background
        function [Ybg, results] = localBG(obj, Ybg, ssub, rr, IND, sn, thresh)
            if ~exist('rr', 'var')||isempty(rr); rr=obj.options.gSiz; end
            if ~exist('ssub', 'var')||isempty(ssub); ssub = 1; end
            if ~exist('IND', 'var') ||isempty(IND); IND = []; end
            if ~exist('sn', 'var')||isempty(sn)
                if isfield(obj.P, 'sn')
                    sn = obj.reshape(obj.P.sn, 2);
                else
                    sn = [];
                end
            else
                sn = obj.reshape(sn, 2);
            end
            
            if ~exist('thresh', 'var')||isempty(thresh); thresh = []; end
            [Ybg, results] = local_background(obj.reshape(Ybg, 2), ssub, rr, IND, sn, thresh);
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
        
        %% reconstruct background signal given the weights
        function Ybg = reconstructBG(obj, Y, weights)
            if ~exist('weights', 'var')||isempty(weights)
                try
                    weights = obj.P.weights;
                catch
                    Ybg = [];
                    disp('no regression weights given');
                end
            end
            Y = obj.reshape(Y-obj.A*obj.C, 2);
            [d1,d2, ~] = size(Y);
            dims = weights.dims;
            b0 = mean(Y,3);
            Y = bsxfun(@minus, Y, b0);
            Y = imresize(Y, dims);
            Y = reshape(Y, [], size(Y,3));
            Ybg = zeros(size(Y));
            parfor m=1:size(Y,1)
                w = weights.weights{m};
                Ybg(m,:) = w(2,:)*Y(w(1,:),:);
            end
            Ybg = bsxfun(@plus, imresize(Ybg, [d1,d2]), b0);
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
        
        %% save results
        function file_path = save_workspace(obj)
            fprintf('------------- SAVE THE WHOLE WORKSPACE ----------\n\n');

            obj.compress_results();
            file_path = [obj.P.log_folder,  strrep(get_date(), ' ', '_'), '.mat'];
            evalin('base', sprintf('save(''%s''); ', file_path));
            
            try
                fp = fopen(obj.P.log_file, 'a');
                fprintf(fp, '\n--------%s--------\n[%s]\bSave the current workspace into file \n\t%s\n\n', get_date(), get_minute(), file_path);
                fprintf('The current workspace has been saved into file \n\t%s\n\n', file_path);
                fp.close();
            end
            
        end
        
        %% compress A, S and W
        function compress_results(obj)
            obj.A = sparse(obj.A);
            obj.S = sparse(obj.S);
            if ~iscell(obj.W)
                obj.W = sparse(obj.W);
            end
        end
        
        %% compute correlation image and peak to noise ratio for endoscopic
        % data. unlike the correlation image for two-photon data,
        % correlation image of the microendoscopic data needs to be
        % spatially filtered first. otherwise neurons are significantly
        % overlapped.
        function [Cn, PNR] = correlation_pnr(obj, Y)
            [Cn, PNR] = correlation_image_endoscope(Y, obj.options);
        end
        
        %% compute correlation image and peak to noise ratio in parallel
        [Cn, PNR] = correlation_pnr_parallel(obj, frame_range)
        
        %% convert struct data to Sources2D object
        function obj = struct2obj(obj, var_struct)
            temp = fieldnames(var_struct);
            for m=1:length(temp)
                try %#ok<*TRYNC>
                    eval(sprintf('obj.%s=var_struct.%s;', temp{m}, temp{m}));
                end
            end
        end
        
        %% convert Sources2D object to a struct variable
        function neuron = obj2struct(obj, ind)
            if ~exist('ind', 'var') || isempty(ind)
                neuron.A = sparse(obj.A);
                neuron.C = obj.C;
                neuron.C_raw = obj.C_raw;
                neuron.S = sparse(obj.S);
                neuron.options = obj.options;
                neuron.P = obj.P;
                neuron.b = obj.b;
                neuron.f = obj.f;
                neuron.W = obj.W;
                neuron.b0 = obj.b0;
                neuron.Fs = obj.Fs;
                neuron.frame_range = obj.frame_range;
                neuron.kernel = obj.kernel;
                neuron.file = obj.kernel;
                neuron.Cn = obj.Cn;
                neuron.ids = obj.ids;
                neuron.tags = obj.tags;
            else
                neuron.A = sparse(obj.A(:, ind));
                neuron.C = obj.C(ind,:);
                neuron.C_raw = obj.C_raw(ind,:);
                neuron.S = sparse(obj.S(ind,:));
                neuron.ids = obj.ids(ind);
                neuron.tags = obj.tags(ind);
            end
            
        end
        
        %% show contours of the all neurons
        function Coor = show_contours(obj, thr, ind, img, with_label)
            %% show neuron contours
            %% inputs:
            %   thr: threshold for the compactness of the neuron
            %   ind: indices of the neurons to be shown
            %   img: the background image. by default, it uses the
            %   correlation image
            %   with_label: include the label of not
            %% outputs:
            %   Coor: contours of all neurosn
            K = size(obj.A, 2);
            if ~exist('ind', 'var') || isempty(ind)
                ind = (1:K);
            end
            if ~exist('thr', 'var') || isempty(thr)
                thr = 0.9;
            end
            if ~exist('img', 'var') || isempty(img)
                img = obj.Cn;
            end
            
            if ~exist('with_label', 'var') || isempty(with_label)
                with_label = false;
            end
            
            if isempty(obj.Coor) || (length(obj.Coor)~=K)
                Coor = obj.get_contours(thr);
                obj.Coor = Coor;
            end
            figure;
            plot_contours(obj.A(:, ind), img, thr,with_label, [], obj.Coor(ind), 2);
            colormap gray;
            file_path = [obj.P.log_folder,  'contours_%dneurons', strrep(get_date(), ' ', '_'), '.pdf'];
            saveas(gcf, file_path);
        end
        
        %% get contours of the all neurons
        function Coor = get_contours(obj, thr, ind)
            A_ = obj.A;
            if exist('ind', 'var')
                A_ = A_(:, ind);
            end
            if ~exist('thr', 'var') || isempty(thr)
                thr = 0.9;
            end
            num_neuron = size(A_,2);
            if num_neuron==0
                Coor ={};
                return;
            else
                Coor = cell(num_neuron,1);
            end
            for m=1:num_neuron
                % smooth the image with median filter
                A_temp = medfilt2(obj.reshape(full(A_(:, m)),2), [3, 3]);
                % find the threshold for detecting nonzero pixels
                
                A_temp = A_temp(:);
                [temp,ind] = sort(A_temp(:).^2,'ascend');
                temp =  cumsum(temp);
                ff = find(temp > (1-thr)*temp(end),1,'first');
                if ~isempty(ff)
                    pvpairs = { 'LevelList' , [0,0]+A_temp(ind(ff)), 'ZData', obj.reshape(A_temp,2)};
                    h = matlab.graphics.chart.primitive.Contour(pvpairs{:});
                    temp = h.ContourMatrix;
                    %temp = medfilt1(temp')';
                    temp = medfilt1(temp(:, 2:end)')';
                    Coor{m} = [temp, temp(:, 1)]; 
                end
                
            end
        end
        
        %% manually draw ROI and show the mean fluorescence traces within the ROI
        function y = drawROI(obj, Y, img, type)
            Y = obj.reshape(Y,1);
            if ~exist('img', 'var') || isempty(img)
                img = mean(Y, 2);
            end
            if ~exist('type', 'var') || isempty(img)
                type = 'ROI';
            end
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            figure;
            while true
                if d1>d2
                    subplot(131);
                else
                    subplot(311);
                end
                obj.image(img);
                if strcmpi(type, 'roi')
                    temp = imfreehand();
                    ind = temp.createMask();
                else
                    [c, r] = ginput(1);
                    ind = sub2ind([d1,d2], round(r), round(c));
                end
                y = mean(Y(ind(:), :), 1);
                if d1>d2
                    subplot(1,3,2:3);
                else
                    subplot(3,1,2:3);
                end
                plot(y);
            end
        end
        %% determine nonzero pixels
        
        %         function Coor = get_contours(obj, thr, ind)
        %             A_ = obj.A;
        %             if exist('ind', 'var')
        %                 A_ = A_(:, ind);
        %             end
        %             if ~exist('thr', 'var') || isempty(thr)
        %                 thr = 0.995;
        %             end
        %             num_neuron = size(A_,2);
        %             if num_neuron==0
        %                 Coor ={};
        %                 return;
        %             else
        %                 Coor = cell(num_neuron,1);
        %             end
        %             for m=1:num_neuron
        %                 % smooth the image with median filter
        %                 img = medfilt2(obj.reshape(full(A_(:, m)),2), [3, 3]);
        %                 % find the threshold for detecting nonzero pixels
        %                 temp = sort(img(img>1e-9));
        %                 if ~any(temp)
        %                     Coor{m} = [];
        %                     continue;
        %                 end
        %                 temp_sum = cumsum(temp);
        %                 ind = find(temp_sum>=temp_sum(end)*(1-thr),1);
        %                 v_thr = temp(ind);
        %
        %                 % find the connected components
        %                 [~, ind_max] = max(img(:));
        %                 temp = bwlabel(img>v_thr);
        %                 img = double(temp==temp(ind_max));
        %                 v_nonzero = imfilter(img, [0,-1/4,0;-1/4,1,-1/4; 0,-1/4,0]);
        %                 vv = v_nonzero(v_nonzero>1e-9)';
        %                 [y, x] = find(v_nonzero>1e-9);
        %                 xmx = bsxfun(@minus, x, x');
        %                 ymy = bsxfun(@minus, y, y');
        %                 dist_pair = xmx.^2 + ymy.^2;
        %                 dist_pair(diag(true(length(x),1))) = inf;
        %                 seq = ones(length(x)+1,1);
        %                 for mm=1:length(x)-1
        %                     [v_min, seq(mm+1)] = min(dist_pair(seq(mm), :)+vv);
        %                     dist_pair(:,seq(mm)) = inf;
        %                     if v_min>3
        %                         seq(mm+1) = 1;
        %                         break;
        %                     end
        %                 end
        %                 Coor{m} = [smooth(x(seq), 2)'; smooth(y(seq),2)'];
        %             end
        %
        %         end
        
        %         %% update background
        %         function [Y, ind_bg, bg] = linearBG(obj, Y)
        %             d1 = obj.options.d1;
        %             d2 = obj.options.d2;
        %             % model the background as linear function
        %             if ndims(Y)==3; Y = neuron.reshape(Y,1); end
        %             T = size(Y,2);
        %             % find background area
        %             ind_frame = round(linspace(1, T, min(T, 500)));
        %             tmp_C1 = correlation_image(Y(:, ind_frame), [1, 2], d1, d2);
        %             tmp_C2 = correlation_image(Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
        %             tmp_Cn = tmp_C1(:)-tmp_C2(:);
        %             ind = (tmp_Cn>0);
        %             ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), 0.1));
        %             % get the mean activity within the selected area
        %             Ybg = mean(Y(ind_bg(:), :), 1);
        %             f1 = (Ybg-mean(Ybg))/std(Ybg);
        %
        %             % regress DY over df to remove the trend
        %             df = diff(f1,1);
        %             b1 = diff(Y, 1,2)*df'/(df*df');
        %
        %             Yres = Y-b1*f1;
        %             b0 = median(Yres, 2);
        %             %             b0 = quantile(Yres, 0.05,2);
        %             Y = bsxfun(@minus, Yres, b0);
        %             bg.b = [b1, b0];
        %             bg.f = [f1; ones(1,T)];
        %         end
        %
        %         %% select background pixels
        %         function [f1, ind_bg] = findBG(obj, Y, q)
        %             d1 = obj.options.d1;
        %             d2 = obj.options.d2;
        %             if nargin<3;    q = 0.1; end;   %quantiles for selecting background pixels
        %             % model the background as linear function
        %             if ndims(Y)==3; Y = neuron.reshape(Y,1); end
        %             T = size(Y,2);
        %             % find background area
        %             ind_frame = round(linspace(1, T, min(T, 500)));
        %             tmp_C1 = correlation_image(Y(:, ind_frame), [1, 2], d1, d2);
        %             tmp_C2 = correlation_image(Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
        %             tmp_Cn = tmp_C1(:)-tmp_C2(:);
        %             ind = (tmp_Cn>0);
        %             ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), q));
        %             % get the mean activity within the selected area
        %             Ybg = mean(Y(ind_bg(:), :), 1);
        %             f1 = (Ybg-mean(Ybg))/std(Ybg);
        %         end
    end
    
end
