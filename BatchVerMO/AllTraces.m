function [AllC, boundary]=AllTraces(neuron_batch)
AllC=[];
boundary=[];
for i=1:length(neuron_batch)
    AllC=[AllC neuron_batch(i).signal];
    boundary=[boundary size(neuron_batch(i).signal,2)];
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