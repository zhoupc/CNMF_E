classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data. It
    % Both CNMF and CNMF-E model can be used.
    
    % Author: Pengcheng Zhou, Columbia University, 2017
    % zhoupc1988@gmail.com
    
    %% properties
    properties
        % spatial
        A;          % spatial components of neurons
        A_prev;     % previous estimation of A
        % temporal
        C;          % temporal components of neurons
        C_prev;     % previous estimation of C
        C_raw;      % raw traces of temporal components
        S;          % spike counts
        kernel;     % calcium dynamics. this one is less used these days.
        % background
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        W;          % a sparse weight matrix matrix
        b0;         % constant baselines for each pixel
        b0_new;  % correct the changes in b0 when we update A & C 
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
        % a1 = 1 means that neuron has too few nonzero pixels
        % a2 = 1 means that neuron has silent calcium transients
        % a3 = 1 indicates that the residual after being deconvolved has
        % much smaller std than the raw data, which is uaually the case
        % when temporal traces are random noise
        % a4 = 1 indicates that the CNMF-E fails in deconvolving temporal
        % traces. this usually happens when neurons are false positives.
        
        % a4 = 1 indicates that the neuron has small PNR
        %others
        Cn;         % correlation image
        PNR;        % peak-to-noise ratio image.
        Coor;       % neuron contours
        neurons_per_patch; 
        Df;         % background for each component to normalize the filtered raw data
        C_df;       % temporal components of neurons and background normalized by Df
        S_df;       % spike counts of neurons normalized by Df
        batches = cell(0);  % results for each small batch data
        file_id = [];    % file id for each batch.
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
                try
                    save .dir.mat dir_nm;
                catch
                    fprintf('You do not have the access to write file to the current folder');
                end
            else
                fprintf('no file was selected. STOP!\n');
                return;
            end
            obj.file = nam;
        end
        
        %% select multiple datasets
        function nams = select_multiple_files(obj, nams)
            if ~isempty(nams)
                nams = unique(nams);
                ind_keep = true(size(nams));
                for m=1:length(nams)
                    if ~exist(nams{m}, 'file')
                        ind_keep = false;
                    else
                        nams{m} = get_fullname(nams{m});
                    end
                end
                
                nams = nams(ind_keep);
                fprintf('\n----------------------------------------\n');
                for m=1:length(nams)
                    fprintf('%d:\t%s\n', m, nams{m});
                end
                fprintf('\n----------------------------------------\n');
                obj.file = nams;
                return;
                
            end
            fprintf('\t-------------------------- GUIDE --------------------------\n');
            fprintf('\ttype -i to remove the selected i-th file\n');
            fprintf('\ttype 0 to stop the data selection\n');
            fprintf('\ttype anything else to select more data\n');
            fprintf('\t--------------------------  END  --------------------------\n\n');
            if ~exist('nams', 'var') || isempty(nams)
                nams = cell(0);
            end
            while true
                if isempty(nams)
                    fprintf('You haven''t selected any file yet\n');
                else
                    ind_keep = true(size(nams));
                    for m=1:length(nams)
                        if ~exist(nams{m}, 'file')
                            ind_keep = false;
                        end
                    end
                    
                    nams = nams(ind_keep);
                    fprintf('\n----------------------------------------\n');
                    for m=1:length(nams)
                        fprintf('%d:\t%s\n', m, nams{m});
                    end
                    fprintf('\n----------------------------------------\n');
                    
                end
                
                % make choise
                
                your_ans = input('make  your choice: ');
                if your_ans==0
                    break;
                elseif your_ans<0 && round(abs(your_ans))<length(nams)
                    nams(round(abs(your_ans))) = [];
                else
                    nams{end+1} = obj.select_data([]);  %#ok<AGROW>
                end
            end
            obj.file = unique(nams);
        end
        
        
        %% load the patched file into the memory 
        function map_data_to_memory(obj)
           mat_file = obj.P.mat_file; 
           if exist(mat_file, 'file')
               data_nam = sprintf('mat_data_%d', string2hash(mat_file));
               if isempty(evalin('base', sprintf('whos(''%s'')', data_nam)))
                   % the data has not been mapped yet 
                   evalin('base', sprintf('%s=load(''%s''); ', data_nam, mat_file));
               end
           end
        end
        
        %% load data within selected patch position and selected frames
        function Ypatch = load_patch_data(obj, patch_pos, frame_range)
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
            if iscell(obj.file)
                nam = obj.file{1};
            else
                nam = obj.file;
            end
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
            
            %% map the data to the memory
            temp = cast(0, data.dtype); 
            temp = whos('temp'); 
            data_space = prod(data.dims) * temp.bytes / (2^30); 
            if data_space < memory_size_to_use
                obj.map_data_to_memory();
            end
        end
        
        %% distribute data and be ready to run batch mode source extraction
        function getReady_batch(obj, pars_envs, log_folder)
            %% parameters for scaling things
            if ~exist('pars_envs', 'var') || isempty(pars_envs)
                pars_envs = struct('memory_size_to_use', 4, ...   % GB, memory space you allow to use in MATLAB
                    'memory_size_per_patch', 0.3, ...   % GB, space for loading data within one patch
                    'patch_dims', [64, 64], ...  %GB, patch size
                    'batch_frames', []);       % number of frames per batch
            end
            batch_frames = pars_envs.batch_frames;
            if ~exist('log_folder', 'var')||isempty(log_folder)||(~exist(log_folder, 'dir'))
                obj.P.log_folder = [cd(), filesep];
            else
                obj.P.log_folder = log_folder;
            end
            %% distribute all files
            nams = unique(obj.file);
            obj.file = nams;
            % pre-allocate spaces for saving results of multiple patches
            batches_ = cell(100, 1);
            file_id_ = zeros(100,1, 'like', uint8(1));
            numFrames = zeros(100,1);
            
            k_batch = 0;
            for m=1:length(nams)
                neuron = obj.copy();
                neuron.file = nams{m};
                neuron.getReady(pars_envs);
                if m==1
                    obj.options = neuron.options;
                end
                T = neuron.P.numFrames;
                if isempty(batch_frames)
                    batch_frames = T;
                end
                nbatches = round(T/batch_frames);
                if nbatches<=1
                    frame_range_t = [1, T];
                else
                    frame_range_t = round(linspace(1, T, nbatches+1));
                end
                nbatches = length(frame_range_t) -1;
                
                for n=1:nbatches
                    k_batch = k_batch +1;
                    tmp_neuron = neuron.copy();
                    tmp_neuron.frame_range = frame_range_t(n:(n+1))-[0, 1*(n~=nbatches)];
                    batch_i.neuron = tmp_neuron;
                    batch_i.shifts = [];  % shifts allowed
                    file_id_(k_batch) = m;
                    batches_{k_batch} = batch_i;
                    numFrames(k_batch) = diff(tmp_neuron.frame_range)+1;
                end
            end
            
            obj.batches = batches_(1:k_batch);
            obj.file_id = file_id_(1:k_batch);
            obj.P.numFrames = numFrames(1:k_batch);
        end
        
        %% estimate noise
        function sn = estimate_noise(obj, frame_range, method)
            mat_data = obj.P.mat_data;
            dims = mat_data.dims;
            T = dims(3);
            if ~exist('frame_range', 'var') || isempty(frame_range)
                frame_range = [1, min(T, 3000)];
            end
            T = diff(frame_range)+1;
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
            for mblock=1:(nr_block*nc_block)
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
        [center, Cn, PNR] = initComponents_parallel(obj, K, frame_range, save_avi, use_parallel, use_prev)
        
        %% initialize neurons for multiple batches
        % for 1P and 2P data, CNMF and CNMF-E
        [center, Cn, PNR] = initComponents_batch(obj, K, frame_range, save_avi, use_parallel)
        
        %% fast initialization for microendoscopic data
        [center, Cn, pnr] = initComponents_endoscope(obj, Y, K, patch_sz, debug_on, save_avi);
        
        [center] = initComponents_2p(obj,Y, K, options, sn, debug_on, save_avi);
        
        %% pick neurons from the residual
        % for 1P and 2P data, CNMF and CNMF-E
        [center, Cn, PNR] = initComponents_residual_parallel(obj, K, save_avi, use_parallel, min_corr, min_pnr, seed_method)
        
        %------------------------------------------------------------------UPDATE MODEL VARIABLES---%
        %% update background components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_background_parallel(obj, use_parallel)
        
        %% update  background for all batches
        update_background_batch(obj, use_parallel)
        
        %% update spatial components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_spatial_parallel(obj, use_parallel, update_sn)
        
        %% update spatial components in batch mode
        update_spatial_batch(obj, use_parallel);
        
        %% update temporal components in parallel
        % for 1P and 2P data, CNMF and CNMF-E
        update_temporal_parallel(obj, use_parallel, use_c_hat)
        
        %% update temporal components in batch mode
        % for 1P and 2P data, CNMF and CNMF-E
        update_temporal_batch(obj, use_parallel)
        
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
        function save_neurons(obj, folder_nm, with_craw)
            if ~exist('folder_nm', 'var') || isempty(folder_nm)
                folder_nm = [obj.P.log_folder, 'neurons'];
            end
            if ~exist('with_craw', 'var')||isempty(with_craw)
                with_craw = true;
            end
            if exist(folder_nm, 'dir')
                temp = cd();
                cd(folder_nm);
                delete *;
                cd(temp);
            else
                mkdir(folder_nm);
            end
            if with_craw
                obj.viewNeurons([], obj.C_raw, folder_nm);
            else
                obj.viewNeurons([], [], folder_nm);
            end
        end
        
        %% export the result as a video
        avi_filename = show_demixed_video(obj,save_avi, kt, frame_range, center_ac, range_ac, range_Y, multi_factor, use_craw)        %% compute the residual
        %         function [Y_res] = residual(obj, Yr)
        %             Y_res = Yr - obj.A*obj.C - obj.b*obj.f;
        %         end
        
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
                    if obj.options.deconv_flag
                        temp = mean(obj.C,2)'.*sum(obj.A);
                    else
                        temp = mean(obj.C,2)'.*sum(obj.A)./obj.P.neuron.sn';
                    end
                    [~, srt] = sort(temp, 'descend');
                elseif strcmp(srt, 'sparsity_spatial')
                    temp = sqrt(sum(obj.A.^2, 1))./sum(abs(obj.A), 1);
                    [~, srt] = sort(temp);
                elseif strcmp(srt, 'sparsity_temporal')
                    temp = sqrt(sum(obj.C_raw.^2, 2))./sum(abs(obj.C_raw), 2);
                    [~, srt] = sort(temp, 'descend');
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
                elseif strcmpi(srt, 'pnr')
                    pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);
                    [~, srt] = sort(pnrs, 'descend');
                else %if strcmpi(srt, 'snr')
                    snrs = var(obj.C, 0, 2)./var(obj.C_raw-obj.C, 0, 2);
                    [~, srt] = sort(snrs, 'descend');
                end
            end
            obj.A = obj.A(:, srt);
            obj.C = obj.C(srt, :);
            
            try
                obj.C_raw = obj.C_raw(srt,:);
                obj.S = obj.S(srt,:);
                obj.P.kernel_pars = obj.P.kernel_pars(srt, :);
                obj.P.neuron_sn = obj.P.neuron_sn(srt);
                
                obj.ids = obj.ids(srt);
                obj.tags = obj.tags(srt);
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
        
        %% concatenate results from multi patches
        function concatenate_temporal_batch(obj)
            T_k= obj.P.numFrames;
            T = sum(T_k);
            K = size(obj.A,2);
            C_ = zeros(K, T);
            C_raw_ = zeros(K, T);
            deconv_flag = obj.options.deconv_flag;
            if deconv_flag
                S_ = zeros(K, T);
            end
            t = 0;
            nbatches = length(obj.batches);
            for mbatch=1:nbatches
                neuron_k = obj.batches{mbatch}.neuron;
                try
                    C_(:, t+(1:T_k(mbatch))) = neuron_k.C;
                catch
                    pause;
                end
                C_raw_(:, t+(1:T_k(mbatch))) = neuron_k.C_raw;
                if deconv_flag
                    S_(:, t+(1:T_k(mbatch))) = neuron_k.S;
                end
                t = t+T_k(mbatch);
            end
            obj.C = C_;
            obj.C_raw = C_raw_;
            if deconv_flag
                obj.S = S_;
            end
        end
        %% quick view
        ind_del = viewNeurons(obj, ind, C2, folder_nm);
        displayNeurons(obj, ind, C2, folder_nm);
        
        %% function remove false positives
        function ids = remove_false_positives(obj, show_delete)
            if ~exist('show_delete', 'var')
                show_delete = false;
            end
            tags_ = obj.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
            ids = find(tags_);
            if isempty(ids)
                fprintf('all components are good \n');
            else
                if show_delete
                    obj.viewNeurons(ids, obj.C_raw);
                else
                    obj.delete(ids);
                end
            end
        end
        
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

                fprintf(flog, '[%s]\b', get_minute());
                fprintf('Deleted %d neurons: ', length(ids_del));
                fprintf(flog, 'Deleted %d neurons: \n', length(ids_del));
                for m=1:length(ids_del)
                    fprintf('%2d, ', ids_del(m));
                    fprintf(flog, '%2d, ', ids_del(m));
                end
                fprintf(flog, '\n');
                fprintf('\nThe IDS of these neurons were saved as intermediate_results.ids_del_%s\n\n', tmp_str);
                if obj.options.save_intermediate
                    tmp_str = get_date();
                    tmp_str=strrep(tmp_str, '-', '_');
                    eval(sprintf('log_data.ids_del_%s = ids_del;', tmp_str));
                    
                    fprintf(flog, '\tThe IDS of these neurons were saved as intermediate_results.ids_del_%s\n\n', tmp_str);
                end
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
        C_ = deconvTemporal(obj, use_parallel, method_noise)
        
        %% decorrelate all tmeporal components 
        C_ = decorrTemporal(obj, wd); 
        
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
                    temp.cdata = imresize(temp.cdata, [height, width]);
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
%             ind_small = find(ind_small);
%             obj.delete(ind_small);
        end
        
        %% keep spatial shapes compact
        function compactSpatial(obj)
            for m=1:size(obj.A, 2)
                ai = obj.reshape(obj.A(:, m), 2);
                ai = circular_constraints(ai);
                obj.A(:, m) = ai(:);
            end
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
        
        %% reconstruct baseline
        function b0_ = reconstruct_b0(obj)
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
                % number of patches
                [nr_patch, nc_patch] = size(patch_pos);
            catch
                error('No data file selected');
                b0_= [];
                return;
            end
            
            b0_ = zeros(d1, d2);
            for mpatch=1:(nr_patch*nc_patch)
                b0_patch = obj.b0{mpatch};
                tmp_patch = patch_pos{mpatch};
                r0 = tmp_patch(1);
                r1 = tmp_patch(2);
                c0 = tmp_patch(3);
                c1 = tmp_patch(4);
                b0_(r0:r1, c0:c1) = reshape(b0_patch, r1-r0+1, c1-c0+1);
            end
        end
        
        %% concatenate W
        function W = concatenate_W(obj)
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
                W= [];
                return;
            end
            d = d1*d2;
            temp = get_nhood(obj.options.ring_radius);
            num_neighbors = length(temp);
            ii_all = zeros(d*num_neighbors, 1);
            jj_all = zeros(d*num_neighbors, 1);
            vv_all = zeros(d*num_neighbors, 1);
            k = 0;
            for mpatch=1:(nr_patch*nc_patch)
                W_patch = obj.W{mpatch};
                tmp_patch = patch_pos{mpatch};
                tmp_block = block_pos{mpatch};
                nr_patch = diff(tmp_patch(1:2))+1;
                nc_patch = diff(tmp_patch(3:4))+1;
                nr_block = diff(tmp_block(1:2))+1;
                nc_block = diff(tmp_block(3:4))+1;
                [ii, jj, vv] = find(W_patch);
                [ii1, ii2] = ind2sub([nr_patch, nc_patch], ii);
                ii1 = ii1 + tmp_patch(1)-1;
                ii2 = ii2 + tmp_patch(3)-1;
                
                [jj1, jj2] = ind2sub([nr_block, nc_block], jj);
                jj1 = jj1 + tmp_block(1)-1;
                jj2 = jj2 + tmp_block(3)-1;
                
                ii = sub2ind([d1, d2], ii1, ii2);
                jj = sub2ind([d1, d2], jj1, jj2);
                ii_all(k+(1:length(ii))) = ii;
                jj_all(k+(1:length(jj))) = jj;
                vv_all(k+(1:length(vv))) = vv;
                
                k = k + length(ii);
            end
            W = sparse(ii_all(1:k), jj_all(1:k), vv_all(1:k), d, d);
        end
        %% reconstruct background
        function Ybg = reconstruct_background(obj, frame_range)
            %%reconstruct background using the saved data
            % input:
            %   frame_range:  [frame_start, frame_end], the range of frames to be loaded
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
            
            if ~exist('frame_range', 'var')||isempty(frame_range)
                frame_range = obj.frame_range;
            end
            % frames to be loaded for initialization
            T = diff(frame_range) + 1;
            
            bg_model = obj.options.background_model;
            bg_ssub = obj.options.bg_ssub;
            % reconstruct the constant baseline
            if strcmpi(bg_model, 'ring')
                b0_ = obj.reconstruct_b0();
                b0_new_ = obj.reshape(obj.b0_new, 2); 
            end
            
            %% start updating the background
            Ybg = zeros(d1, d2, T);
            for mpatch=1:(nr_patch*nc_patch)
                tmp_patch = patch_pos{mpatch};
                if strcmpi(bg_model, 'ring')
                    W_ring = obj.W{mpatch};
                    %                     b0_ring = obj.b0{mpatch};
                    % load data
                    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
                    [nr_block, nc_block, ~] = size(Ypatch);
                    Ypatch = reshape(Ypatch, [], T);
                    tmp_block = block_pos{mpatch};
                    tmp_patch = patch_pos{mpatch};
                    b0_ring = b0_(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4));
                    b0_ring = reshape(b0_ring, [], 1);
                    
                    b0_patch = reshape(b0_new_(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)), [], 1);
                    
                    % find the neurons that are within the block
                    mask = zeros(d1, d2);
                    mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
                    ind = (reshape(mask(:), 1, [])* obj.A_prev>0);
                    
                    A_patch = obj.A_prev(logical(mask), ind);
                    C_patch = obj.C_prev(ind, frame_range(1):frame_range(2));
                    
                    % reconstruct background
                    %                     Cmean = mean(C_patch , 2);
                    Ypatch = bsxfun(@minus, double(Ypatch), b0_ring);
%                     b0_ring = b0_(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4));
%                     b0_ring = reshape(b0_ring, [], 1);
%                     
                    if bg_ssub==1
                        Bf = W_ring*(double(Ypatch) - A_patch*C_patch);
                        Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, Bf, b0_patch), diff(tmp_patch(1:2))+1, [], T);
                    else
                        [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                        temp = reshape(double(Ypatch)-A_patch*C_patch, nr_block, nc_block, []);
                        temp = imresize(temp, 1./bg_ssub, 'nearest');
                        Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                        Bf = imresize(Bf, [nr_block, nc_block], 'nearest');
                        Bf = Bf((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1, :);
                        Bf = reshape(Bf, [], T);
                        Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, Bf, b0_patch), diff(tmp_patch(1:2))+1, [], T);
                    end
                elseif strcmpi(bg_model, 'nmf')
                    b_nmf = obj.b{mpatch};
                    f_nmf = obj.f{mpatch};
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(b_nmf*f_nmf(:, frame_range(1):frame_range(2)), diff(tmp_patch(1:2))+1, [], T);
                else
                    b_svd = obj.b{mpatch};
                    f_svd = obj.f{mpatch};
                    b0_svd = obj.b0{mpatch};
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, b_svd*f_svd(:, frame_range(1):frame_range(2)), b0_svd), diff(tmp_patch(1:2))+1, [], T);
                end
                
            end
            
        end
        
        % compute RSS
        function [RSS_total, RSS] = compute_RSS(obj, frame_range)
            %% compute RSS
            % input:
            %   frame_range:  [frame_start, frame_end], the range of frames to be loaded
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
            
            if ~exist('frame_range', 'var')||isempty(frame_range)
                frame_range = obj.frame_range;
            end
            % frames to be loaded for initialization
            T = diff(frame_range) + 1;
            
            bg_model = obj.options.background_model;
            bg_ssub = obj.options.bg_ssub;
            % reconstruct the constant baseline
            if strcmpi(bg_model, 'ring')
                b0_ = obj.reconstruct_b0();
                b0_new_ = obj.reshape(obj.b0_new, 2); 
            end
            
            %% start reconstructing the background
            RSS = cell(nr_patch, nc_patch);
            W_all = obj.W;
            b0_all = cell(nr_patch, nc_patch);  % including b0 within block
            b0_all_new = b0_all; 
            A_all = cell(nr_patch, nc_patch);
            C_all = cell(nr_patch, nc_patch);
            A_all_prev = cell(nr_patch, nc_patch); 
            C_all_prev = cell(nr_patch, nc_patch); 
            b_ = obj.b;
            f_ = obj.f;
            for mpatch=1:nr_patch*nc_patch
                if strcmpi(bg_model, 'ring')
                    tmp_block = block_pos{mpatch};
                else
                    tmp_block = patch_pos{mpatch};
                end
                b0_all{mpatch} = b0_(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4));
                b0_all_new{mpatch} = b0_new_(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4));
                % find the neurons that are within the block
                mask = zeros(d1, d2);
                mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
                ind = (reshape(mask(:), 1, [])* obj.A>0);
                A_all{mpatch} = obj.A(logical(mask), ind);
                C_all{mpatch} = obj.C(ind, frame_range(1):frame_range(2));
                if strcmpi(bg_model, 'ring')
                    ind = (reshape(mask(:), 1, [])* obj.A_prev>0);
                    A_all_prev{mpatch} = obj.A_prev(logical(mask), ind);
                    C_all_prev{mpatch} = obj.C_prev(ind, frame_range(1):frame_range(2));
                end
            end
            
            %%
           parfor mpatch=1:(nr_patch*nc_patch)
                tmp_block = block_pos{mpatch};
                tmp_patch = patch_pos{mpatch};
                
                if strcmpi(bg_model, 'ring')
                    W_ring = W_all{mpatch};
                    mask = zeros(d1, d2);
                    mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
                    mask(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)) = 2;
                    ind_patch = (mask(mask>0)==2);
                    
                    %                     b0_ring = obj.b0{mpatch};
                    % load data
                    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
                    [nr_block, nc_block, ~] = size(Ypatch);
                    Ypatch = reshape(Ypatch, [], T);
                    
%                     tmp_block = block_pos{mpatch};
%                     tmp_patch = patch_pos{mpatch};
                    b0_ring = b0_all{mpatch};
                    b0_new_patch = b0_all_new{mpatch};
                    b0_ring = reshape(b0_ring, [], 1);
%                     b0_new_patch = reshape(b0_new_patch, [], 1); 
                    
                    % find the neurons that are within the block
                    A_patch = A_all{mpatch};
                    C_patch = C_all{mpatch};
                    
                    A_patch_prev = A_all_prev{mpatch};
                    C_patch_prev = C_all_prev{mpatch};
                    
                    % compute Y-A*C
                    YmAC = double(Ypatch(ind_patch, :)) -A_patch(ind_patch, :)*C_patch;
                    b0_new_patch = reshape(b0_new_patch(ind_patch), [], 1); 
                    
                    % reconstruct background
                    Ypatch = bsxfun(@minus, double(Ypatch), b0_ring);
%                     b0_ring = b0_all_patch{mpatch};
%                     b0_ring = reshape(b0_ring, [], 1);
                    if bg_ssub==1                        
                        Bf = W_ring*(double(Ypatch) - A_patch_prev*C_patch_prev);
                    else
                        [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                        temp = reshape(double(Ypatch)-A_patch_prev*C_patch_prev, nr_block, nc_block, []);
                        temp = imresize(temp, 1./bg_ssub, 'nearest');
                        Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                        Bf = imresize(Bf, [nr_block, nc_block], 'nearest');
                        Bf = Bf((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1, :);
                        Bf = reshape(Bf, [], T);
                    end
%                     b0_ring = mean(YmAC, 2); 
                    Ybg_patch = bsxfun(@plus, Bf, b0_new_patch); %, diff(tmp_patch(1:2))+1, [], T);                  
                elseif strcmpi(bg_model, 'nmf')
                    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
                    YmAC = double(Ypatch) - A_all{mpatch}*C_all{mpatch};
                    b_nmf = b_{mpatch};
                    f_nmf = f_{mpatch};
                    Ybg_patch =b_nmf*f_nmf(:, frame_range(1):frame_range(2));
                else
                    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
                    YmAC = double(Ypatch) - A_all{mpatch}*C_all{mpatch};
                    b_svd = b_{mpatch};
                    f_svd = f_{mpatch};
                    b0_svd = mean(YmAC, 2);
                    Ybg_patch = bsxfun(@plus, b_svd*f_svd(:, frame_range(1):frame_range(2)), b0_svd);
                end
                RSS{mpatch} = sum((YmAC(:)-Ybg_patch(:)).^2);
%                 temp = YmAC - Ybg_patch; 
%                 RSS{mpatch} = sum(temp(:).^2)- size(temp,2)*sum(mean(temp, 2).^2); 
                
            end
            temp = cell2mat(RSS);
            RSS_total = sum(temp(:));
            obj.P.RSS = RSS_total; 
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
        function [center, Cn, pnr] = pickNeurons(obj, Y, patch_par, seed_method, debug_on, K)
            if ~exist('patch_par', 'var')||isempty(patch_par)
                patch_par = [3,3];
            end
            if ~exist('seed_method', 'var')||isempty(seed_method)
                seed_method = 'auto';
            end
            if ~exist('debug_on', 'var')||isempty(debug_on)
                debug_on = false;
            end
            if ~exist('K', 'var') || isempty(K)
                K = [];
            end
            neuron = obj.copy();
            neuron.options.seed_method = seed_method;
            neuron.options.gSig = 1;
            neuron.options.center_psf = 0;
            [center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, false);
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
        function tags_ = tag_neurons_parallel(obj, min_pnr)
            if ~exist('min_pnr', 'var') || isempty(min_pnr)
                min_pnr = 3;
            end
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
                nz_spikes = full(sum(S_(:,2:end)>0, 2));
                tags_ = tags_ + uint16(nz_spikes<1)*2;
                
                tmp_std = std(obj.C_raw-obj.C, 0, 2);
                % check the noise level, it shouldn't be 0
                temp= tmp_std./GetSn(obj.C_raw);
                tags_ = tags_ + uint16(temp<0.1)*4;
                
                % check PNR of neural traces
                pnrs = max(obj.C, [], 2)./tmp_std;
                tags_ = tags_+uint16(pnrs<min_pnr)*8;
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
            evalin('base', sprintf('save(''%s'', ''neuron'', ''save_*'', ''show_*'', ''use_parallel'', ''with_*'', ''-v7.3''); ', file_path));
            try
                fp = fopen(obj.P.log_file, 'a');
                fprintf(fp, '\n--------%s--------\n[%s]\bSave the current workspace into file \n\t%s\n\n', get_date(), get_minute(), file_path);
                fprintf('The current workspace has been saved into file \n\t%s\n\n', file_path);
                fp.close();
            end
            
        end
        
        %% save results for all batches
        function file_paths = save_workspace_batch(obj, log_folder)
            if ~exist('log_folder', 'var')||isempty(log_folder)||(~exist(log_folder, 'dir'))
                log_folder = [cd(), filesep];
            end
            obj.P.log_folder = log_folder;
            nbatches = length(obj.batches);
            file_paths = cell(nbatches, 1);
            
            for mbatch=1:nbatches
                batch_k = obj.batches{mbatch};
                neuron = batch_k.neuron;
                
                fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);
                
                % update background
                neuron.compress_results();
                file_path = [neuron.P.log_folder,  strrep(get_date(), ' ', '_'), '.mat'];
                
                save(file_path, 'neuron');
                try
                    fp = fopen(neuron.P.log_file, 'a');
                    fprintf(fp, '\n--------%s--------\n[%s]\bSave the current workspace into file \n\t%s\n\n', get_date(), get_minute(), file_path);
                    fprintf('The current workspace has been saved into file \n\t%s\n\n', file_path);
                    fp.close();
                end
                file_paths{mbatch} = file_path;
                
            end
            obj.P.file_paths = file_paths;
            obj.save_workspace();
        end
        
        %% clean the log folders
        function clean_results(obj)
            if ~isempty(obj.batches)
                obj.clean_results_bach();
            end
            try
                rmdir(obj.P.log_folder, 's');
                fprintf('results have been cleaned\n');
            catch
                fprintf('file deletion failed\n');
            end
        end
        
        %% clean the log folders for all patch
        function clean_results_batch(obj)
            nbatches = length(obj.batches);
            for mbatch=1:nbatches
                batch_k = obj.batches{mbatch};
                neuron_k = batch_k.neuron;
                
                fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);
                
                % update background
                neuron_k.clean_results();
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
        
        %% compute correlation image and PNR image
        function [Cn, PNR] = correlation_pnr(obj, Y)
            [Cn, PNR] = correlation_image_endoscope(Y, obj.options);
        end
        
        %% compute correlation image and peak to noise ratio in parallel
        [Cn, PNR] = correlation_pnr_parallel(obj, frame_range)
        
        %% compute correlation image and PNR image for each batch
        function correlation_pnr_batch(obj)
            nbatches = length(obj.batches);
            for mbatch=1:nbatches
                fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);
                batch_k = obj.batches{mbatch};
                neuron_k = batch_k.neuron;
                [neuron_k.Cn, neuron_k.PNR] = neuron_k.correlation_pnr_parallel();
                batch_k.neuron = neuron_k;
                obj.batches{mbatch} = batch_k;
            end
            
        end
        
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
            else
                Coor = obj.Coor;
            end
            figure('papersize', [obj.options.d2, obj.options.d1]/40);
            init_fig;
            plot_contours(obj.A(:, ind), img, thr,with_label, [], obj.Coor(ind), 2);
            colormap gray;
            try
                file_path = [obj.P.log_folder,  'contours_neurons', strrep(get_date(), ' ', '_'), '.pdf'];
                saveas(gcf, file_path);
            end
        end
        
        %% get contours of the all neurons
        function Coor = get_contours(obj, thr, ind_show)
            A_ = obj.A;
            if exist('ind_show', 'var')
                A_ = A_(:, ind_show);
            else
                ind_show = 1:size(A_, 2);
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
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            %             tmp_kernel = strel('square', 3);
            for m=1:num_neuron
                % smooth the image with median filter
                A_temp = obj.reshape(full(A_(:, m)),2);
                % find the threshold for detecting nonzero pixels
                
                A_temp = A_temp(:);
                [temp,ind] = sort(A_temp(:).^2,'ascend');
                temp =  cumsum(temp);
                ff = find(temp > (1-thr)*temp(end),1,'first');
                thr_a = A_temp(ind(ff));
                A_temp = obj.reshape(A_temp,2);
                
                % crop a small region for computing contours
                [tmp1, tmp2, ~] = find(A_temp);
                if isempty(tmp1)
                    Coor{m} = zeros(2,1);
                    continue;
                end
                rmin = max(1, min(tmp1)-3);
                rmax = min(d1, max(tmp1)+3);
                cmin = max(1, min(tmp2)-3);
                cmax = min(d2, max(tmp2)+3);
                A_temp = A_temp(rmin:rmax, cmin:cmax);
                
                l = bwlabel(medfilt2(A_temp>thr_a));
                l_most = mode(l(l>0));
                ind = (l==l_most); 
                A_temp(ind) =  max(A_temp(ind), thr_a); 
                A_temp(~ind) = min(A_temp(~ind), thr_a*0.99); 
                
                pvpairs = { 'LevelList' , thr_a, 'ZData', A_temp};
                h = matlab.graphics.chart.primitive.Contour(pvpairs{:});
                temp = h.ContourMatrix;
                if isempty(temp)
                    temp = obj.get_contours((thr+1)/2, ind_show(m));
                    Coor{m} = temp{1};
                    continue;
                else
                    temp(:, 1) = temp(:, 2);
                    temp = medfilt1(temp')';
                    temp(:, 1) = temp(:, end);
                    Coor{m} = bsxfun(@plus, temp, [cmin-1; rmin-1]);
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
        %% determine the number of neurons within each patch 
        function nn = ids_per_patch(obj)
            % calculate how many neurons in each patch 
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
                % number of patches
                [nr_patch, nc_patch] = size(patch_pos);
            catch
                error('No data file selected');
            end
            
            nn = cell(size(obj.b)); 
            for mpatch=1:(nr_patch*nc_patch)
                tmp_patch = patch_pos{mpatch};
                r0 = tmp_patch(1);
                r1 = tmp_patch(2);
                c0 = tmp_patch(3);
                c1 = tmp_patch(4);
                ind = false(d1,d2); 
                ind(r0:r1, c0:c1) = true;
                nn{mpatch} = obj.ids(sum(obj.A(ind, :), 1)>0); 
            end
            obj.neurons_per_patch = nn; 
        end
        
        
     
        
    end
    
end
