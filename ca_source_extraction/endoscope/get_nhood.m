function [r_shift, c_shift] = get_nhood(radius)
% find the pixel locations 
%% determine neibours of each pixel
rsub = (-radius):(radius);      % row subscript
csub = rsub;      % column subscript
[cind, rind] = meshgrid(csub, rsub);
R = sqrt(cind.^2+rind.^2);
neigh_kernel = (R>=radius) .* (R<radius+1);  % kernel representing the selected neighbors 

[r_shift, c_shift] = find(neigh_kernel);
r_shift = reshape(r_shift - radius -1, 1, []);
c_shift = reshape(c_shift - radius - 1, 1, []);
