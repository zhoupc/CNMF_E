function weightedA=ReducingA(Aunique,STDunique)
    A=cat(2,Aunique{:});    
        STDsum=sum(cat(3,STDunique{:}),3);
        aveSTD=cellfun(@rdivide,STDunique,STDsum,'UniformOutput', false);
        aveSTD=cat(2,aveSTD{:});
    catSTD=diag(aveSTD);
    
    weightedA=A*catSTD;
        d=size(A,1);
        K=size(Aunique{1},2);
        
    weightedA=reshape(A,d,K,[]);
    weightedA=sum(weightedA,3);
end
