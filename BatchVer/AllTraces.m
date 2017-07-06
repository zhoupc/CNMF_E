function AllC=AllTraces(neuron_batch)
AllC=[];
for i=1:length(neuron_batch)
    AllC=[AllC neuron_batch(i).signal];
end