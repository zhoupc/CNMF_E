function Atemp=centralA(Atemp)
for i=1:size(Atemp,2)
    ai=Atemp(:,i);
    temp = full(ai>quantile(ai, 0.5, 1));
    ai(~temp(:)) = 0;
    Atemp(:,i)=ai;
end
end