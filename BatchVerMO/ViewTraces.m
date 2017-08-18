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
neuronums=[1,60];
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
[AllC, ~]=AllTracesfromACS(ACS);
AllC=AllC(:,:);

[AllC_raw, ~]=AllTracesfromACS(ACS,true);
AllC_raw=AllC_raw(:,:);

%newIDs=Day4_newIDs;
newIDs_nums=[1];
neuronums=[];
for i=1:length(newIDs_nums)
    neuronums=[neuronums; newIDs{newIDs_nums(i)}];
end
%%
neuronums=[1 57 19 60];
K=length(neuronums);
figure
for i=1:K
    I=neuronums(i);
    subplot(K,1,i)
    plot(AllC_raw(I,:),'r')
    hold on
    plot(AllC(I,:),'b')
    title(['original neuron ' num2str(I)])

end


%%
K=length(neuronums);
for i=1:K
    I=neuronums(i);
    subplot(K,1,i)
    hold on
    Areshaped=reshape(A(:,I),300,400);
    imagesc(Areshaped)
    title(['original neuron ' num2str(I)])
end
%%
newIDs_nums=[20,42];
newIDs=DayAll_newIDs;
cum_sum=cumsum(cellfun(@(x) size(x,2),M1{1}));
file_num=[];
actual_num=[];
neuronums=[];
for i=1:length(newIDs_nums)
    neuronums=[neuronums; newIDs{newIDs_nums(i)}];
    for j=1:length(newIDs{newIDs_nums(i)})
        x_=newIDs{newIDs_nums(i)};
        x=x_(j);
        k=find((cum_sum>=x),1);
        file_num=[file_num k];
        if k>1        
            x_actual=x-cum_sum(k-1);
        else
            x_actual=x;
        end
        actual_num=[actual_num x_actual];
    end
end
%%

K=length(neuronums);
for i=1:K
    I=neuronums(i);
    subplot(K,1,i)
    hold on
    plot(AllC(I,:))
    title(['original neuron ' num2str(I)])
end

%% Just plotting from ACS
neurons=[20,21];
K=length(neurons);
for i =1:K
    subplot(K,1,i)
    plot(ACS(1).Cin(neurons(i),:))
end


