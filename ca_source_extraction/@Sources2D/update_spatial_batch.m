function update_spatial_batch(obj, use_parallel)
%% update the the spatial components for all batches and then summarize the results 
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% update spatial components for all batches 
nbatches = length(obj.batches); 

[d, K] = size(obj.A); 
if K==0
    fprintf('You have to initialize neurons first\n'); 
    return; 
end
A = zeros(d, K); 
cc = zeros(nbatches, K); 

for mbatch=1:nbatches
    batch_k = obj.batches{mbatch}; 
    neuron_k = batch_k.neuron; 
    
    fprintf('\nprocessing batch %d/%d\n', mbatch, nbatches);

    % update temporal components 
    neuron_k.update_spatial_parallel(use_parallel); 
    cc(mbatch, :) = sum(neuron_k.C.^2, 2); 
    A = A + bsxfun(@times, neuron_k.A, cc(mbatch, :));    
end

obj.A = bsxfun(@times, A, 1./sum(cc, 1)); 
obj.post_process_spatial(); 

%% spread the same A to all batches and update temporal components
for mbatch=1:(length(obj.batches)-1)
    batch_k = obj.batches{mbatch}; 
    neuron_k = batch_k.neuron; 

    neuron_k.A = obj.A;

    batch_k.neuron = neuron_k;
    obj.batches{mbatch} = batch_k;
end