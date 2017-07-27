[AllC,boundary]=AllTraces(neuron_batch);
%Boundary=boundary;

figure; 
plot(AllC(1,:))
hold on
plot(boundary'*ones(1,2), [0 max(AllC(1,:),[],2)], 'k:')

    PartC{ni}=C;
    boundary{ni}=cumsum(bound);
    

figure; 
plot(PartC{14})
hold on
plot(boundary{14}'*ones(1,2), [0 max(PartC{14},[],2)], 'k:')