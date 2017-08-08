function [AllC, boundary]=AllTracesfromACS(ACS)
AllC=[];
boundary=[];
for i=1:length(ACS)
    AllC=[AllC ACS(i).Cin];
    boundary=[boundary size(ACS(i).Cin,2)];
end
boundary=cumsum(boundary);
end