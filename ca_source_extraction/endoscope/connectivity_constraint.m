function img = connectivity_constraint(img, thr, sz)
% remove small nonzero pixels
if nargin<2;    thr = 0.01; end;
if nargin<3;    sz = 5; end;
se = strel('square', sz);
ai_open = imopen(img, se);

temp = full(ai_open>max(img(:))*thr);
l = bwlabel(temp, 4);   % remove disconnected components
[~, ind_max] = max(ai_open(:));

img(l(:)~=l(ind_max)) = 0;