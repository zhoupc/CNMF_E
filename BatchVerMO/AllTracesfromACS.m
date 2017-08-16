function [AllC, boundary]=AllTracesfromACS(ACS)
AllC=[];
boundary=[];
for i=1:length(ACS)
    AllC=[AllC ACS(i).Cin_raw];
    boundary=[boundary size(ACS(i).Cin_raw,2)];
end
boundary=cumsum(boundary);
end