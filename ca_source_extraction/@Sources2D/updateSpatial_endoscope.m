function updateSpatial_endoscope(obj, Y, numIter)
% number of iterations
if ~exist('numIter', 'var')||isempty(numIter)
    numIter = 10;
end
% determine the search locations
method = obj.options.search_method; 
params = obj.options; 
IND = logical(determine_search_location(obj.A, method, params));
% idx = find(sum(IND,2)>0); 
% CC = obj.C*obj.C'; 
% Cy = obj.C*Y(idx,:)'; 
% for m=1:length(idx)
%     ind = IND(idx(m), :); 
%     obj.A(idx(m),ind) = nnls(CC(ind,ind), Cy(ind'), obj.A(idx(m), ind)', [], 3); 
% end 
% run HALS
obj.A = HALS_spatial(Y, obj.A, obj.C, IND, numIter);

% thresholding the minimum number of neurons
obj.delete(sum(obj.A, 1)<=obj.options.min_pixel);
end
