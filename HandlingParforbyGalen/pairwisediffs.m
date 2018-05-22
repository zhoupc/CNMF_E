function diffMatOrCell = pairwisediffs(ptA, ptB)
%PAIRWISEDIFFS finds all pairwise differences between two matrices
%	For two matrices with the same number of rows, PAIRWISEDIFFS finds
%	the distance between all pairs of points from matrix A and matrix B.
%
%   [diffMat] = PAIRWISEDIFFS(ptA) returns the difference between each
%   pair of points a row of the matrix ptA, excluding the difference
%   between each point and itself.
%
%   [diffMat} = PAIRWISEDIFFS(ptA, ptB) returns the differences
%   between every point in a row of ptA and all the other points in the
%   corresponding row of ptB. ptA and ptB must have the same number of
%   rows, or one must only have one row.
%
%   Written by Galen Lynch 10/04/2014
%   email: glynch@mit.edu
nRepA = size(ptA, 1);
nPtA = size(ptA, 2);
if ~exist('ptB', 'var') %Auto
    assert(~iscell(ptA), 'expecting a matrix');
    if nPtA < 2
        diffMatOrCell = nan(nRepA, 0);
        return
    end
    nDiffs = nchoosek(nPtA, 2);
    diffMatOrCell = nan(nRepA, nDiffs);
    pos = 1;
    for ptNo = 1:(nPtA - 1)
        thisNDiff = nPtA - ptNo;
        diffMatOrCell(:, pos:(pos+ thisNDiff - 1)) = ...
            bsxfun(@minus, ptA(:, (ptNo + 1):end), ptA(:, ptNo));
        pos = pos+ thisNDiff;
    end
    diffMatOrCell = abs(diffMatOrCell);
else
    nRepB = size(ptB, 1);
    hasSingleton = nRepB == 1 || nRepA == 1;
    assert(nRepA == nRepB || hasSingleton,...
        'matrices must have same size in the first dimension');
    diffMatOrCell = reshape(bsxfun(@minus, ptA,...
        permute(ptB, [1, 3, 2])), nRepA, []);
end
end