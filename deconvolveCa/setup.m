oasis_folder = fileparts(mfilename('fullpath')); 
addpath(sprintf('%s', oasis_folder)); 
addpath(sprintf('%s%sfunctions', oasis_folder, filesep)); 
addpath(sprintf('%s%soasis', oasis_folder, filesep)); 
addpath(sprintf('%s%soasis_kernel', oasis_folder, filesep)); 
addpath(sprintf('%s%sMCMC', oasis_folder, filesep)); 
addpath(sprintf('%s%sMCMC%sutilities', oasis_folder, filesep, filesep)); 

%% install convex optimization solvers
optimization_folder = sprintf('%s%soptimization', oasis_folder, filesep); 
if ~exist(optimization_folder, 'dir'); 
    mkdir(optimization_folder);
end

% install cvx 
if isempty(which('cvx_begin.m'))
    if ~exist('cvx', 'dir')
        %install cvx
        if ismac
            cvx_url = 'http://web.cvxr.com/cvx/cvx-maci64.zip';
        elseif isunix
            cvx_url = 'http://web.cvxr.com/cvx/cvx-a64.zip';
        elseif ispc
            cvx_url = 'http://web.cvxr.com/cvx/cvx-w64.zip';
        else
            fprints('Your platform is not supported by CVX\n');
            return;
        end
        fprintf('Downloading CVX...\n');
        unzip(cvx_url, optimization_folder);
    end
    run(sprintf('%s%scvx%scvx_setup', optimization_folder, filesep, filesep));
end

%% save the current path 
%savepath(); 
