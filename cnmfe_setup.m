cnmfe_folder = fileparts(mfilename('fullpath'));

if ~exist('cnmfe_loaded', 'var') || ~cnmfe_loaded
    addpath(fullfile(cnmfe_folder, 'ca_source_extraction'));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'utilities')));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'endoscope')));
    addpath(fullfile(cnmfe_folder, 'GUI'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'gui_callbacks'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'modules'));
    addpath(fullfile(cnmfe_folder, 'cnmfe_scripts'));
    cnmfe_loaded = true;
end

%% install deconvolution package 
run(fullfile(cnmfe_folder, 'OASIS_matlab', 'oasis_setup.m'));

%% install convex optimization solvers
% by default, we don't install cvx any more. if you want to install cvx,
% then set install_cvx = true and then run oasis-setup.m 
if exist('install_cvx', 'var') && install_cvx && isempty(which('cvx_begin.m'))
    if ~exist(fullfile(oasis_folder, 'packages', 'cvx'), 'dir')
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
        unzip(cvx_url, fullfile(oasis_folder, 'packages'));
    end
    run(fullfile(oasis_folder, 'packages', 'cvx', 'cvx_setup.m'));
end
%% save path
%savepath();

%% deconvolution 
