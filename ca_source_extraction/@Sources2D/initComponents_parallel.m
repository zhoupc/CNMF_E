function [center, Cn, PNR] = initComponents_parallel(obj, K, frame_range, save_avi)
%% initializing spatial/temporal components for miceoendoscopic data
%% input:
%   K:  scalar, maximum number of neurons
%   frame_range: 1 X 2 vector indicating the starting and ending frames 
%   save_avi: save the video of initialization procedure
%% Output:
%   center: d*2 matrix, centers of all initialized neurons.
%   Cn:     correlation image
%   PNR:    peak to noise ratio
%% Author: Pengcheng Zhou, Columbia University, 2017 
%% email: zhoupc1988@gmail.com

%% process parameters
get_date = @() datestr(datetime(), 'mmm-dd_HH'); 
get_minute = @() datestr(datetime(), 'HH_MM_SS'); 

try 
    % map data 
    mat_data = obj.mat_data; 
    mat_file = mat_data.Properties.Source; 
    tmp_dir = fileparts(mat_file); 
    log_folder = [tmp_dir, filesep, 'logs_', get_date(), filesep]; 
    if ~exist('log_folder', 'dir')
        mkdir(log_folder); 
    end 
    log_file = [log_folder, 'logs.txt'];
    obj.P.log_file = log_file; 
    
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
if ~exist('frame_range', 'var') || isempty(frame_range) 
    frame_range = [1, T]; 
else
    frame_range(frame_range<1) = 1; 
    frame_range(frame_range>T) = T; 
end
T = diff(frame_range) + 1; 

% maximum neuron number in each patch
if (~exist('K', 'var')) || (isempty(K))
    % if K is not specified, use a very large number as default
    K = round((d1*d2));
end

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false; %don't save initialization procedure
end


% parameters for detrending the data first. 
if isfield(obj.options, 'nk') % number of knots for creating spline basis
    nk = obj.options.nk;
else
    nk = 1;
end

% parameter for avoiding using boundaries pixels as seed pixels 
options = obj.options;
if ~isfield(options, 'bd') || isempty(options.bd')
    options.bd = options.gSiz;   % boundary pixesl to be ignored during the process of detecting seed pixels
end
bd = options.bd; 

%% start initialization 
% save the log infomation 
tmp_file = [log_folder, get_minute(), '.mat']; 
save(tmp_file, 'options'); 
flog = fopen(log_file, 'w'); 
fprintf(flog, 'Data: %s\n\n', mat_file); 

fprintf(flog, '%s\b', get_minute()); 
fprintf(flog, 'Start running source extraction......\nThe options are saved into \n %s\n\n', tmp_file); 

fprintf(flog, '%s\b', get_minute()); 
fprintf(flog, 'Start initializing neurons...\n\n'); 

Ain = cell(nr_patch, nc_patch); % save spatial components of neurons in each patch
Cin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, denoised trace
Sin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, deconvolved trace
Cin_raw = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, raw trace
kernel_pars = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch
center = cell(nr_patch, nc_patch);     % save centers of all initialized neurons
Cn = zeros(d1, d2);
PNR = zeros(d1, d2);
default_kernel = obj.kernel; 

results = cell(nr_patch*nc_patch, 1); 
parfor mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch 
    [r, c] = ind2sub([nr_patch, nc_patch], mpatch); 
    tmp_patch = patch_pos{r, c};  %#ok<PFBNS>
    tmp_block = block_pos{r, c};  %#ok<PFBNS>
    
    % boundaries pixels to be avoided for detecting seed pixels 
    tmp_options = options;  
    tmp_options.visible_off = true; 
    tmp_options.bd = max([(tmp_patch-tmp_block).*[1, -1, 1, -1]; bd, bd, bd, bd], [], 1); 
    
    % patch dimension
    tmp_options.d1 = diff(tmp_block(1:2))+1;
    tmp_options.d2 = diff(tmp_block(3:4))+1;
    
    % parameter for calcium indicators. This one may not be used if the
    % selected deconvolution algorithm doesn't need it
    tmp_options.kernel = default_kernel;
    
    % file names for saving avi file
    if save_avi
        tmp_save_avi = sprintf('%sinitialization_%d_%d_%d_%d.avi', log_folder, tmp_block(1), tmp_block(2), tmp_block(3), tmp_block(4));
    else
        tmp_save_avi = save_avi;
    end
       
    % load the patch data 
    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
    Ypatch = double(reshape(Ypatch, [], T)); 
    if nk>1
        Ypatch_dt = detrend_data(Ypatch, nk); % detrend data
        [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch_dt, K, tmp_options, [], tmp_save_avi); 
    else
        [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch, K, tmp_options, [], tmp_save_avi);
    end
    
    % put everthing into one struct variable 
    tmp_results.center = tmp_center; 
    tmp_results.Cn = tmp_Cn; 
    tmp_results.PNR = tmp_PNR; 
    results{mpatch} = tmp_results; 
%     eval(sprintf('results_patch_%d=tmp_results;', mpatch));  %#ok<PFBFN>
    % display initialization progress
    fprintf('Patch (%d, %d) is done. %d X %d patches in total. \n\n', r, c, nr_patch, nc_patch);
end 

%% collect results 
for mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch
    [mr, mc] = ind2sub([nr_patch, nc_patch], mpatch);
    tmp_patch = patch_pos{mr, mc};
    tmp_block = block_pos{mr, mc};
    r0 = tmp_block(1); 
    r1 = tmp_block(2); 
    c0 = tmp_block(3); 
    c1 = tmp_block(4);  
    ind_patch = true(r1-r0+1, c1-c0+1); 
    ind_patch((tmp_patch(1):tmp_patch(2))-r0+1, (tmp_patch(3):tmp_patch(4))-c0+1) = false; 
    
    % unpack results
    tmp_results = results{mpatch}; 
%     eval(sprintf('tmp_results=results_patch_%d;', mpatch));
    tmp_Ain = tmp_results.Ain;
    tmp_Ain(ind_patch, :) = 0;
    tmp_Cin = tmp_results.Cin;
    tmp_Cin_raw = tmp_results.Cin_raw;
    tmp_center = tmp_results.center ;
    tmp_Cn = tmp_results.Cn;
    tmp_PNR = tmp_results.PNR;
    if options.deconv_flag
        tmp_Sin = tmp_results.Sin;
        tmp_kernel_pars = tmp_results.kernel_pars;
    end
    tmp_K = size(tmp_Ain, 2);   % number of neurons within the selected patch
    [tmp_d1, tmp_d2] = size(tmp_Cn); 
    
    temp = zeros(d1, d2, tmp_K);  % spatial components of all neurons
    temp(r0:r1, c0:c1, :) = reshape(full(tmp_Ain), tmp_d1, tmp_d2, []);
    Ain{mr, mc} = reshape(temp, d1*d2, tmp_K);
    Cin{mr, mc} = tmp_Cin;      % temporal components of all neurons
    Cin_raw{mr, mc} = tmp_Cin_raw;
    if options.deconv_flag
        Sin{mr, mc} = tmp_Sin;
        kernel_pars{mr,mc} = tmp_kernel_pars;
    end
    center{mr, mc} = bsxfun(@plus, tmp_center, [r0-1, c0-1]); % centers
    Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn);
    PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR);
end

%% export the results
Ain = cell2mat(reshape(Ain, 1, []));
Cin = cell2mat(reshape(Cin, [], 1));
Cin_raw = cell2mat(reshape(Cin_raw, [], 1));
center = cell2mat(reshape(center, [], 1));
obj.A = sparse(Ain);
obj.C = Cin;
obj.C_raw = Cin_raw;
if options.deconv_flag
    Sin = cell2mat(reshape(Sin, [], 1));
    kernel_pars = cell2mat(reshape(kernel_pars, [], 1));
    obj.S = sparse(Sin);
    obj.P.kernel_pars = kernel_pars;
else
    obj.S = zeros(size(obj.C));
end
obj.Cn = Cn;

%% save the results to log 
neuron = obj.obj2struct();  %#ok<NASGU>
tmp_file = [log_folder, get_minute(), '.mat']; 
save(tmp_file, 'neuron'); 
fprintf(flog, '%s\b', get_minute()); 
fprintf(flog, 'Finished the initialization procedure.\n'); 
fprintf(flog, 'In total, %d neurons were initialized. \n', size(Ain,2));
fprintf(flog, 'The initialization results were saved into \n %s\n\n', tmp_file); 
fclose(flog); 