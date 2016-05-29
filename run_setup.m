%% add path
CNMF_dir = fileparts(which('run_setup.m'));
addpath(sprintf('%s%sca_source_extraction%s', CNMF_dir, filesep, filesep));
addpath(sprintf('%s%sca_source_extraction%sutilities%s', CNMF_dir, filesep, filesep, filesep));
addpath(sprintf('%s%sca_source_extraction%sendoscope%s', CNMF_dir, filesep, filesep, filesep));

%% setup cvx
try 
    %test whether cvx has been installed 
    cvx_begin
    cvx_end
catch
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
    unzip(cvx_url, CNMF_dir);
    run(sprintf('%s%scvx%scvx_setup', CNMF_dir, filesep, filesep));
end
%% save path
savepath();