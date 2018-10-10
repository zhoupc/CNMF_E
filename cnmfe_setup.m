cnmfe_folder = fileparts(mfilename('fullpath'));

if ~exist('cnmfe_loaded', 'var') || ~cnmfe_loaded
    addpath(fullfile(cnmfe_folder, 'ca_source_extraction'));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'utilities')));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'endoscope')));
    addpath(fullfile(cnmfe_folder, 'GUI'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'gui_callbacks'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'modules'));
    addpath(fullfile(cnmfe_folder, 'scripts'));
    cnmfe_loaded = true;
end

%% install deconvolution package
oasis_folder = fullfile(cnmfe_folder, 'OASIS_matlab'); 
if exist(fullfile(oasis_folder, 'oasis_setup.m'), 'file')
    run(fullfile(cnmfe_folder, 'OASIS_matlab', 'oasis_setup.m'));
else
    oasis_url = 'https://github.com/zhoupc/OASIS_matlab/archive/master.zip'; 
    fprintf('downloading OASIS_matlab....\n'); 
    unzip(oasis_url, cnmfe_folder); 
    if exist(oasis_folder, 'dir')
        rmdir(oasis_folder, 's');
    end
    movefile(fullfile(cnmfe_folder, 'OASIS_matlab-master'), ...
        oasis_folder);
    fprintf('done!\n'); 
end
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
