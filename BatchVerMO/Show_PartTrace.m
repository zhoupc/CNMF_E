neuronnum=3;
figure
ax(1)=subplot(2,1,1);
plot(PartC_raw{neuronnum})
hold on
plot(boundary_raw{neuronnum}'*[1 1], [0 max(PartC_raw{neuronnum})], ':', 'color', .7*[1 1 1]);
% plot(AllC(neuronnum,:))
% hold on
% plot(boundary_raw{neuronnum}'*[1 1], [0 max(AllC(neuronnum,:))], ':', 'color', .7*[1 1 1]);

ax(2)=subplot(2,1,2);
plot(PartC{neuronnum})
hold on
plot(boundary_raw{neuronnum}'*[1 1], [0 max(PartC{neuronnum})], ':', 'color', .7*[1 1 1]);
linkaxes(ax)