function weightedA=ReducingA(File)
    A=cat(2,File.Aunique);
    
        STD=cat(2,File.STDunique);
        STDsum=sum(cat(3,File.STDunique));
        aveSTD=bsxfun(@rdivide,STD,STDsum);
    catSTD=diag(aveSTD);
    
    weightedA=A*catSTD;
        d=size(A,1);
        K=size(File(1).Aunique,2);
        
    weightedA=reshape(A,d,K,[]);
    weightedA=sum(weightedA,3);
end
