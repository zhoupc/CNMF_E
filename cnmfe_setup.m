%% add path
CNMF_dir = fileparts(which('cnmfe_setup.m'));
addpath(fullfile(CNMF_dir, 'ca_source_extraction'));
addpath(genpath(fullfile(CNMF_dir, 'ca_source_extraction', 'utilities')));
addpath(genpath(fullfile(CNMF_dir, 'ca_source_extraction', 'endoscope')));
addpath(fullfile(CNMF_dir, 'GUI'));
addpath(fullfile(CNMF_dir, 'GUI', 'gui_callbacks'));
addpath(fullfile(CNMF_dir, 'GUI', 'modules'));
addpath(fullfile(CNMF_dir, 'cnmfe_scripts'));

%% setup cvx
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
            fprintf('Your platform is not supported by CVX\n');
            return;
        end
        fprintf('Downloading CVX...\n');
        unzip(cvx_url, CNMF_dir);
    end
    run(sprintf('%s%scvx%scvx_setup', CNMF_dir, filesep, filesep));
end
%% save path
%savepath();

%% deconvolution 
run(fullfile(CNMF_dir, 'deconvolveCa', 'setup.m'));
