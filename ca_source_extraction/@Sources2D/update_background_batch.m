function update_background_batch(obj, use_parallel)
%% update  background for all patches
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% update spatial components for all batches 
nbatches = length(obj.batches); 

for mbatch=1:nbatches
    batch_k = obj.batches{mbatch}; 
    neuron_k = batch_k.neuron; 
    
    fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches); 

    % update background
    neuron_k.update_background_parallel(use_parallel); 
    
    batch_k.neuron = neuron_k; 
    obj.batches{mbatch} = batch_k; 
end