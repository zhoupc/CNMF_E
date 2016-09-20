function updateSpatial_endoscope(obj, Y, num, method)
%% udpate spatial components

%% inputs:
%   Y: d X T matrix, data
%   num: scalar. If method=='hals', then num is the number of iterations to
%       update A; If method=='nnls', then num is the maximum number of neurons
%       overlapping at one pixel
%   method: method for updating the spatial components {'hals', 'nnls'}.
%       default: 'nnls'

%% Author: Pengcheng Zhou, Carnegie Mellon University.

%% input parameters number of iterations
if ~exist('method', 'var')||isempty(method)
    method = 'nnls';
end
if ~exist('num', 'var')||isempty(num)
    if strcmpi(method, 'nnls')
        num=5;
    else
        num = 10;
    end
end

%% determine the search locations
search_method = obj.options.search_method;
params = obj.options;
if strcmpi(search_method, 'dilate')
    obj.options.se = [];
end
IND = logical(determine_search_location(obj.A, search_method, params));

%% update spatial components
if strcmpi(method, 'hals')
    obj.A = HALS_spatial(Y, obj.A, obj.C, IND, num);
else
    obj.A = nnls_spatial(Y, obj.A, obj.C, IND, num);
end

%% thresholding the minimum number of neurons
obj.delete(sum(obj.A, 1)<=obj.options.min_pixel);

%% post-process
% obj.post_process_spatial();
% if strcmpi(method, 'nnls')
%     IND = bsxfun(@gt, obj.A, max(obj.A, [], 1)/100);
%     obj.A = nnls_spatial(Y, obj.A, obj.C, IND, num);
% end

end
