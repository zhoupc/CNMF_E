function weightedA=ReducingA(Aunique,STDunique)
%Description: Reduce A for those that cannot be merged but has different values in each
%             sampled file. General method is use std as weight to merge
%             these.
%Input: two cells each contains Aunique,STDunique of each file.
%Output:One big weightedA.
%author: Shijie Gu, techel@live.cn, (ShanghaiTech University; Fee Lab at MIT-BCS)
    A=cat(2,Aunique{:});    
        STDsum=sum(cat(3,STDunique{:}),3);
        STDcat=cat(2,STDunique{:});
        STDsumcat=repmat(STDsum,1,length(STDunique));
        aveSTD=STDcat./STDsumcat;
        %aveSTD=cellfun(@rdivide,STDunique,STDsum,'UniformOutput', false);
    catSTD=diag(aveSTD);
    
    weightedA=A*catSTD;
        d=size(A,1);
        K=size(Aunique{1},2);
        
    weightedA=reshape(weightedA,d,K,[]);
    weightedA=sum(weightedA,3);
end
