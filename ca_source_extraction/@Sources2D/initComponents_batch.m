function initComponents_batch(obj, K, save_avi, use_parallel)
%% initializing spatial/temporal components for calcium imaging data
%% input:
%   K:  scalar, maximum number of neurons
%   save_avi: save the video of initialization procedure
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Output:
%   center: d*2 matrix, centers of all initialized neurons.
%   Cn:     correlation image
%   PNR:    peak to noise ratio
%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% process parameters
% maximum neuron number in each patch
if (~exist('K', 'var')) || (isempty(K))
    % if K is not specified, use a very large number as default
    K = [];
end

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false; %don't save initialization procedure
elseif save_avi
    use_parallel = false;
end

% use parallel or not
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true; %don't save initialization procedure
end

%% initializing neurons from all batches.
nbatches = length(obj.batches);
for mbatch=1:nbatches
    batch_k = obj.batches{mbatch};
    neuron_k = batch_k.neuron;
    neuron_k.options = obj.options;
    fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);
    
    if mbatch==1
        % initialization neurons from nothing
        neuron_k.initComponents_parallel(K, neuron_k.frame_range, save_avi, use_parallel);
        
        W_init = neuron_k.W;
        b_init = neuron_k.b;
        f_init = neuron_k.f;
        b0_init = neuron_k.b0;
    else
        % copy previous initialization of neuron shape
        neuron_k.A = obj.A;
        neuron_k.P.k_ids = obj.P.k_ids;
        neuron_k.W = obj.W;
        neuron_k.b = obj.b;
        neuron_k.ids = obj.ids;
        neuron_k.tags = obj.tags;
        
        neuron_k.W = W_init;
        neuron_k.b = b_init;
        neuron_k.f = f_init;
        neuron_k.b0 = b0_init;
        
        neuron_k.initTemporal(neuron_k.frame_range, use_parallel);
    end
    % update background
    neuron_k.update_background_parallel(use_parallel);
    
    % pick more neurons from the residual
    neuron_k.initComponents_residual_parallel([], save_avi, use_parallel);
    
    % collect results
    batch_k.neuron = neuron_k;
    obj.A = neuron_k.A;
    obj.W = neuron_k.W;
    obj.b = neuron_k.b;
    obj.f = neuron_k.f;
    obj.b0 = neuron_k.b0;
    obj.ids = neuron_k.ids;
    obj.tags = neuron_k.tags;
    obj.P.k_ids = neuron_k.P.k_ids; 
    obj.batches{mbatch} = batch_k;
end

%% spread the same A to all batches and update temporal components
K = size(obj.A, 2);
nbatches = length(obj.batches);
for mbatch=1:nbatches
    batch_k = obj.batches{mbatch};
    neuron_k = batch_k.neuron;
    neuron_k.options = obj.options;
    
    if size(neuron_k.A,2) == size(obj.A,2)
        break;
    end
    neuron_k.A = obj.A;
    neuron_k.P.k_ids = obj.P.k_ids;
    neuron_k.ids = obj.ids;
    neuron_k.tags = obj.tags;
    
    fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);
    [tmp_K, T] = size(neuron_k.C);
    if K>tmp_K
        neuron_k.C = [neuron_k.C; zeros(K-tmp_K, T)];
    end
    % update temporal components
    neuron_k.update_temporal_parallel(use_parallel);
    
    % collect results
    batch_k.neuron = neuron_k;
    obj.batches{mbatch} = batch_k;
end
