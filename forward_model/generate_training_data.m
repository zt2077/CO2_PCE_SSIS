function generate_training_data(root, output_dir)
% Evaluate the archived input design with the supplied COMSOL model.
% Start MATLAB with COMSOL LiveLink and connect to a COMSOL server first.
assert(exist('mphopen','file') == 2, 'COMSOL LiveLink for MATLAB is required.');
global model
model = mphopen(fullfile(root,'models','co2_storage_modify1.mph'));
model.hist.disable;
inputs = load(fullfile(root,'data','training','X_data_sobol_188all1.mat'),'X');
coords = load(fullfile(root,'data','training','coordinate.mat'),'coordinate');
X = inputs.X;
save(fullfile(output_dir,'X_generated.mat'),'X');
support_unit(X, size(coords.coordinate,2), 20, output_dir);
fprintf('Forward responses saved separately from the archived training data.\n');
end
