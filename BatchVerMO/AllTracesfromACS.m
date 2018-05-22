function [AllC, boundary]=AllTracesfromACS(ACS,rawOrNot)
if nargin<2
    rawOrNot=false;
end
AllC=[];
boundary=[];
for i=1:length(ACS)
    if rawOrNot
        AllC=[AllC ACS(i).Cin_raw];
        boundary=[boundary size(ACS(i).Cin_raw,2)];
    else
        AllC=[AllC ACS(i).Cin];
        boundary=[boundary size(ACS(i).Cin,2)];
    end
    boundary=cumsum(boundary);
end

    
end