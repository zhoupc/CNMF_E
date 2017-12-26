function update_temporal_batch(obj, use_parallel)
%% update the the spatial components for all patches
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% update spatial components for all batches 
nbatches = length(obj.batches); 

if isempty(obj.A)
    fprintf('You have to initialize neurons first\n'); 
    return; 
end
K = size(obj.A, 2); 

for mbatch=1:nbatches
    batch_k = obj.batches{mbatch}; 
    neuron_k = batch_k.neuron; 
    
    fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches); 
    [tmp_K, T] = size(neuron_k.C); 
    if K>tmp_K
        neuron_k.C = [neuron_k.C; zeros(K-tmp_K, T)]; 
    end
    % update temporal components 
    neuron_k.update_temporal_parallel(use_parallel); 
    
    batch_k.neuron = neuron_k; 
    obj.batches{mbatch} = batch_k; 
end