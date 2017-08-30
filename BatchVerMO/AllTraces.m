function [AllC, boundary]=AllTraces(neuron_batch,type)
if nargin<2||isempty(type);
    type='rawsignal';
end

AllC=[];
boundary=[];
if strcmp(type,'rawsignal')
    for i=1:length(neuron_batch)
        AllC=[AllC neuron_batch(i).rawsignal];
        boundary=[boundary size(neuron_batch(i).rawsignal,2)];
    end
elseif strcmp(type,'signal')
    for i=1:length(neuron_batch)
        AllC=[AllC neuron_batch(i).signal];
        boundary=[boundary size(neuron_batch(i).signal,2)];
    end
end
boundary=cumsum(boundary);
end
% %%
% framesize=[];
% framesizeC=[];
% for i=1:length(neuron_batch)
%     framesize=[framesize size(neuron_batch(i).signal,2)];
%     framesizeC=[framesizeC size(neuron_batch(i).neuron.C,2)];   
% end