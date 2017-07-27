% K=size(neuron_batch(1).neuron.A,2);
% AllC=AllTraces(neuron_batch);
% for i=1:K
%     subplot(K,1,i)
%     hold on
%     plot(AllC(i,:))
% end

% newIDs=sampleresult.newIDs;
% A0s=sampleresult.A0s;
% [~,size2]=cellfun(@size,A0s);
% sizecumsum=cumsum(size2);

%C=cat(2,ACS.Cin);
neuronums=[12,60];
K=length(neuronums);
for i=1:K
    I=neuronums(i);
    subplot(K,1,i)
    hold on
    plot(C(I,:))
    title(['original neuron ' num2str(I)])
end
%%
%A=cat(2,A0s{:});
neuronums=[20,25,61,63];
K=length(neuronums);
for i=1:K
    I=neuronums(i);
    subplot(K,1,i)
    hold on
    Areshaped=reshape(A(:,I),300,400);
    imagesc(Areshaped)
    title(['original neuron ' num2str(I)])
end