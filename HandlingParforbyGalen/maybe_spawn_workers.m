function [poolObj, ownPool, poolSize] = maybe_spawn_workers(nWorkersStrOrNumeric)
% written by Galen Lynch, fixes issue with parfor on openmind.
%% Create workers if needed
poolObj = gcp('nocreate');
if isempty(poolObj)
    %% Check to see if numWorkers needs to be converted from string to double
    if ischar(nWorkersStrOrNumeric)
        nWorkers = str2double(nWorkersStrOrNumeric);
    else
        nWorkers = nWorkersStrOrNumeric;
    end
    if nWorkers > 1
        ownPool = true;
        c = parcluster();
        t = tempname();
        mkdir(t);
        c.JobStorageLocation = t;
        poolObj = parpool(c, nWorkers);
        poolSize = poolObj.NumWorkers;
    else
        ownPool = false;
        poolSize = 0;
    end
else
    ownPool = false;
    poolSize = poolObj.NumWorkers;
end
end