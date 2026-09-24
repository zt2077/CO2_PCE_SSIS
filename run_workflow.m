function output_dir = run_workflow(task, result_dir)
% RUN_WORKFLOW Run a packaged workflow without writing to the input directories.
%   run_workflow('figures') uses the bundled reference results.
%   run_workflow('analysis') fits surrogates and runs MCS and SSIS.
%   run_workflow('figures', RESULT_DIR) plots a newly generated analysis.
%   run_workflow('benchmark') runs the illustrative two-dimensional example.
%   run_workflow('forward') evaluates the bundled training inputs in COMSOL.
if nargin < 1, task = 'figures'; end
root = fileparts(mfilename('fullpath'));
if nargin < 2, result_dir = fullfile(root, 'data', 'reference_results'); end
result_dir = char(result_dir);
if ~isfolder(result_dir), error('Result directory not found: %s', result_dir); end
[ok, info] = fileattrib(result_dir);
assert(ok, 'Cannot resolve the result directory.');
result_dir = info.Name;
old_dir = pwd;
old_path = path;
cleanup = onCleanup(@() restore_environment(old_dir, old_path));
addpath(fullfile(root,'analysis'), fullfile(root,'visualization'), ...
    fullfile(root,'examples'), fullfile(root,'forward_model'), ...
    fullfile(root,'data','training'));
assert(license('test','Statistics_Toolbox'), ...
    'Statistics and Machine Learning Toolbox is required.');
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
output_dir = fullfile(root, 'outputs', [char(task) '_' stamp]);
assert(~isfolder(output_dir), 'Output directory already exists; retry.');
mkdir(output_dir);
cd(output_dir);
switch lower(char(task))
    case 'analysis'
        assert(exist('uqlab','file') == 2, ...
            'Install UQLab and add its core directory to the MATLAB path.');
        run_analysis;
    case 'figures'
        required = {'IS_results_r20.00.mat','MCS_results_r20.00.mat', ...
                    'comparison_data_r20.00.mat'};
        for k = 1:numel(required)
            assert(isfile(fullfile(result_dir,required{k})), ...
                'Missing result file: %s', required{k});
        end
        addpath(result_dir, '-begin');
        plot_comparison;
        plot_spatial;
        plot_gmm_slices;
    case 'benchmark'
        toy_example;
    case 'forward'
        generate_training_data(root, output_dir);
    otherwise
        error('Unknown task. Choose figures, analysis, benchmark, or forward.');
end
fprintf('\nResults written to:\n%s\n', output_dir);
end

function restore_environment(old_dir, old_path)
cd(old_dir);
path(old_path);
end
