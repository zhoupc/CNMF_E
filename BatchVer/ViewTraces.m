K=size(neuron_batch(1).neuron.A,2);
AllC=AllTraces(neuron_batch);
for i=1:K
    subplot(K,1,i)
    hold on
    plot(AllC(i,:))
end