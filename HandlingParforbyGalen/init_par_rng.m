function init_par_rng(seedNo, varargin)
persistent p;
if isempty(p)
    p = inputParser();
    addParameter(p, 'method', 'simdTwister');
end
parse(p, varargin{:});
Opt = p.Results;

rng(seedNo, Opt.method);
%% Create workers if needed
poolObj = gcp();
poolSize = poolObj.NumWorkers;
badSeed = true;
while badSeed
    randSeeds = randi(intmax, 1, poolSize);
    seedDiffs = pairwisediffs(randSeeds);
    badSeed = any(seedDiffs == 0);
end
randOffsets = distributed(randSeeds);
spmd
    rng(getLocalPart(randOffsets), Opt.method);
end
end