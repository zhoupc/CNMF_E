function Y = get_patch_data(mat_data, patch_pos, frame_range, with_overlap)
%% pull data for a mat file where the video data is saved into multiple blocks.
%% inputs:
%   patch_pos: 1x4 vector, [r0, r1, c0, c1], the patch location [r0:r1,
%   c0:c1]


%% author: Pengcheng Zhou

%% parameters
if ~exist('with_overlap', 'var') || isempty(with_overlap)
    with_overlap = false;
end

% check whether the code is running on a worker (for parallel processing)
try
    isOnWorker = ~isempty(getCurrentTask());
catch
    isOnWorker = false;
end

%% check whether the mat_data has been loaded into the base space 
mat_file = mat_data.Properties.Source();
data_nam = sprintf('mat_data_%d', string2hash(mat_file));
if isOnWorker || isempty(evalin('base', sprintf('whos(''%s'')', data_nam)))
    ws = 'caller';   % use the current workspace
    data_nam = 'mat_data'; % load data from mat file
    isOnWorker = true; 
else
    ws = 'base';
end

dims = mat_data.dims;
d1 = dims(1);
d2 = dims(2);
T = dims(3); 
if ~exist('patch_pos', 'var') || isempty(patch_pos)
    patch_pos = [1, d1, 1, d2]; 
end 
if ~exist('frame_range', 'var') || isempty(frame_range)
    frame_range = [1, T]; 
end 
w_overlap = mat_data.w_overlap;

idx_r = mat_data.block_idx_r;
idx_c = mat_data.block_idx_c;
dtype = mat_data.dtype;

% determine the smallest block that includes the batch 
r0_patch = patch_pos(1); 
r1_patch = patch_pos(2); 
c0_patch = patch_pos(3); 
c1_patch = patch_pos(4); 
ind_r0 = find(idx_r<=r0_patch, 1, 'last'); 
ind_r1 = find(idx_r>=r1_patch, 1, 'first'); 
ind_c0 = find(idx_c<=c0_patch, 1, 'last'); 
ind_c1 = find(idx_c>=c1_patch, 1, 'first'); 

% indices for block coordinates
block_idx_r = idx_r(ind_r0:ind_r1);
block_idx_c = idx_c(ind_c0:ind_c1); 

%% pre-allocate a matrix
nr = block_idx_r(end)-block_idx_r(1) + 1;
nc = block_idx_c(end)-block_idx_c(1) + 1;
nframes = diff(frame_range)+1;
nr_block = length(block_idx_r) - 1;
nc_block = length(block_idx_c) - 1;

Y = zeros(nr, nc, nframes, 'like', cast(0,dtype));
block_rstart = block_idx_r(1)-1;
block_cstart = block_idx_c(1)-1;
for m=1:nr_block
    r0 = block_idx_r(m);
    r1 = block_idx_r(m+1);
    nr = r1-r0+1; %#ok<*NASGU>
    for n=1:nc_block
        c0 = block_idx_c(n);
        c1 = block_idx_c(n+1);
        nc = c1-c0+1;
        if isOnWorker
            Y((r0:r1)-block_rstart, (c0:c1)-block_cstart, :) =  eval(sprintf('mat_data.Y_%d_%d_%d_%d(1:nr, 1:nc, %d:%d);', r0, r1, c0, c1, frame_range(1), frame_range(2)));
        else
            tmp_str = sprintf('%s.Y_%d_%d_%d_%d(1:%d, 1:%d, %d:%d);', data_nam, r0, r1, c0, c1, nr, nc, frame_range(1), frame_range(2));
            Y((r0:r1)-block_rstart, (c0:c1)-block_cstart, :) = evalin(ws, tmp_str);
        end
    end
end

%% remove the overlapping area
if ~with_overlap
    Y = Y((r0_patch:r1_patch)-block_rstart, (c0_patch:c1_patch)-block_cstart, :);
end


%% convert string to hash values 
function hash=string2hash(str,type)
% This function generates a hash value from a text string
%
% hash=string2hash(str,type);
%
% inputs,
%   str : The text string, or array with text strings.
% outputs,
%   hash : The hash value, integer value between 0 and 2^32-1
%   type : Type of has 'djb2' (default) or 'sdbm'
%
% From c-code on : http://www.cse.yorku.ca/~oz/hash.html 
%
% djb2
%  this algorithm was first reported by dan bernstein many years ago 
%  in comp.lang.c
%
% sdbm
%  this algorithm was created for sdbm (a public-domain reimplementation of
%  ndbm) database library. it was found to do well in scrambling bits, 
%  causing better distribution of the keys and fewer splits. it also happens
%  to be a good general hashing function with good distribution.
%
% example,
%
%  hash=string2hash('hello world');
%  disp(hash);
%
% Function is written by D.Kroon University of Twente (June 2010)


% From string to double array
str=double(str);
if(nargin<2), type='djb2'; end
switch(type)
    case 'djb2'
        hash = 5381*ones(size(str,1),1); 
        for i=1:size(str,2), 
            hash = mod(hash * 33 + str(:,i), 2^32-1); 
        end
    case 'sdbm'
        hash = zeros(size(str,1),1);
        for i=1:size(str,2), 
            hash = mod(hash * 65599 + str(:,i), 2^32-1);
        end
    otherwise
        error('string_hash:inputs','unknown type');
end
