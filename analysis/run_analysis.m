function run_analysis()
%% Analysis Script for SS+IS+GMM vs MCS Comparison
% This script runs the analysis and saves results for later visualization
% Run this script only when you need to perform new analysis


rng(100,'twister');
uqlab

%% 1. Load coordinate data
coords = load('coordinate.mat');
coords = coords.coordinate;  % 3 × Nnodes
xn = coords(1,:);  yn = coords(2,:);  zn = coords(3,:);
Nnodes = size(coords,2);
%Nnodes = 100;
use_coords = false; % Include node coordinates as PCE inputs.

%% 2. Analysis parameters
injection_ratios = [20];
nRatios = length(injection_ratios);
Pf_values = zeros(nRatios,Nnodes);

%% 3. Main analysis loop
for i = 1:nRatios
    r = injection_ratios(i);
    fprintf('Processing injection ratio r = %.2f\n', r);
    %Load data
    load('X_data_sobol_188all1.mat');
    load('Yall_data_sobol_188all1.mat');
    load('S1_sobol_188all1.mat');
    load('S3_sobol_188all1.mat'); 
    % load('X_data_combined.mat');
    % load('Yall_data_combined.mat'); 
    % load('S1_combined.mat');        
    % load('S3_combined.mat');        
    % X = X_comb;
    % Yall = Yall_comb;
    % dataS1_all = S1_comb;
    % dataS3_all = S3_comb;
    %NN = 64;
      NN = 100;
    X = X(1:NN,:);
    Yall = Yall(1:NN,:);
    dataS1_all = dataS1_all(1:NN,:);
    dataS3_all = dataS3_all(1:NN,:);
    YY = [Yall(:,1:Nnodes), dataS1_all(:,1:Nnodes), dataS3_all(:,1:Nnodes)];
    Y_only = Yall;
    
    %% 3.1 Build PCE surrogate models
    loo_max_per_node = nan(Nnodes,1);     % Maximum LOO error at each node.
    loo_all_per_node = nan(Nnodes,3);     % Per-output LOO errors at each node (three outputs).
    mse_max_per_node = nan(Nnodes,1);     % Maximum MSE at each node.
    mse_all_per_node = nan(Nnodes,3);     % Per-output MSE at each node.
    cv5_max_per_node = nan(Nnodes,1);     % Maximum five-fold CV error at each node.
    cv5_all_per_node = nan(Nnodes,3);     % Per-output five-fold CV errors at each node.
    hybridModels = cell(Nnodes, 1);      % Cell array of hybrid model structures.
    pce_degrees = nan(Nnodes,1);          % Selected PCE degree at each node.
    pce_terms = nan(Nnodes,1);            % Number of PCE terms at each node.
    
    % Store S3 normalization parameters.
    s3_normalization = struct();
    s3_normalization.mean_vals = nan(Nnodes,1);
    s3_normalization.std_vals = nan(Nnodes,1);
    
    fprintf('\nBuilding hybrid surrogate models for %d nodes...\n', Nnodes);
    fprintf('Modeling strategy: PCE for Y and S1; adaptive constant/linear/PCE for S3.\n');
    fprintf('S3 models are selected according to relative signal strength.\n');
    
    for j = 1:Nnodes
        if mod(j,10) == 1
            fprintf('Progress: %d/%d nodes\n', j, Nnodes);
        end
        
        if use_coords == true
            xj = coords(1, j); yj = coords(2, j); zj = coords(3, j);
            xyz = repmat([xj, yj, zj], size(X,1), 1);
            Xj = [X, xyz];
        else
            Xj = X;
        end
        
        % Extract raw data.
        Y_raw  = YY(:, j);                    % Displacement.
        S1_raw = YY(:, j+Nnodes);            % First principal stress.
        S3_raw = YY(:, j+2*Nnodes);          % Third principal stress.
        
        % % Standardize S3.
        % S3_mean = mean(S3_raw);
        % S3_std = std(S3_raw);
        
        % % Store normalization parameters.
        % s3_normalization.mean_vals(j) = S3_mean;
        % s3_normalization.std_vals(j) = S3_std;
        
        % % Standardize S3 while avoiding division by zero.
        % if S3_std > eps
        %     S3_normalized = (S3_raw - S3_mean) / S3_std;
        % else
        %     S3_normalized = S3_raw - S3_mean;  % Center only if the standard deviation is zero.
        % end
        
        % Assemble PCE training responses (S3 normalization is disabled above).
        Yj = [Y_raw, S1_raw, S3_raw];
        
        % Define input marginals
        PriorOpts = struct();
        PriorOpts.Marginals(1).Name = 'E';
        PriorOpts.Marginals(1).Type = 'Gaussian';
        PriorOpts.Marginals(1).Moments = [8.9e9 2e7];
        PriorOpts.Marginals(2).Name = 'nu';
        PriorOpts.Marginals(2).Type = 'Gaussian';
        PriorOpts.Marginals(2).Moments = [0.25 0.02];
        PriorOpts.Marginals(3).Name = 'cohesion';
        PriorOpts.Marginals(3).Type = 'Gaussian';
        PriorOpts.Marginals(3).Moments = [3e6 5e5];
        PriorOpts.Marginals(4).Name = 'internalphi';
        PriorOpts.Marginals(4).Type = 'Gaussian';
        PriorOpts.Marginals(4).Moments = [(3.14/180)*28.4 (3.14/180)*1];
        myPriorDist = uq_createInput(PriorOpts);
        
        % Build PCE models for Y and S1 and select an adaptive model for S3.
        % PCE options for Y and S1
        MetaOpts_PCE = struct();
        MetaOpts_PCE.Type = 'Metamodel';
        MetaOpts_PCE.MetaType = 'PCE';
        MetaOpts_PCE.Method = 'LARS';
        MetaOpts_PCE.TruncOptions.qNorm = 0.6;
        MetaOpts_PCE.TruncOptions.MaxInteraction = 2;
        %MetaOpts_PCE.Degree = 2:10;
        MetaOpts_PCE.Degree = 6:8;
        MetaOpts_PCE.ExpDesign.X = Xj;
        MetaOpts_PCE.ExpDesign.Y = [Y_raw, S1_raw];  % Use PCE for Y and S1 only.
        MetaOpts_PCE.Name = sprintf('PCE_YS1_Node_%d_%d', i,j);
        MetaOpts_PCE.Display = 'quiet';
        
        % Select the modeling method according to signal strength.
        S3_range = max(S3_raw) - min(S3_raw);
        S3_mean_abs = abs(mean(S3_raw));
        S3_cv = std(S3_raw) / S3_mean_abs;  % Coefficient of variation.
        S3_relative_range = S3_range / S3_mean_abs;  % Relative range.
        % fprintf('Node %d: S3 relative range: %.4f\n', j, S3_relative_range);
        
        % Assess S3 signal strength.
        is_weak_signal = S3_relative_range < 0.05;  % Relative variation below 5%.
        is_very_weak_signal = S3_relative_range < 0.02;  % Relative variation below 2%.
        
        % Select the S3 modeling method.
        if is_very_weak_signal
            % Very weak signal: use a stable mean model.
            S3_model_method = 'constant';
            S3_pred_value = mean(S3_raw);
            pceModel_S3 = [];
            % fprintf('Node %d: very weak signal; using the mean model.\n', j);
            
        elseif is_weak_signal
            % Weak signal: use simple, stable linear regression.
            S3_model_method = 'linear';
            % Linear regression Y = a*X + b.
            X_design = [ones(size(Xj,1),1), Xj];  % Add an intercept.
            S3_coeffs = X_design \ S3_raw;  % Least-squares solution.
            pceModel_S3 = [];
            % fprintf('Node %d: weak signal; using a linear model.\n', j);
        else
            % Normal signal: use conservative PCE settings.
            S3_model_method = 'PCE';
            MetaOpts_S3 = struct();
            MetaOpts_S3.Type = 'Metamodel';
            MetaOpts_S3.MetaType = 'PCE';
            MetaOpts_S3.Method = 'LARS';
            MetaOpts_S3.TruncOptions.qNorm = 0.5;  % More conservative truncation.
            MetaOpts_S3.TruncOptions.MaxInteraction = 2;
            MetaOpts_S3.Degree = 2:7;  % Lower polynomial degrees.
            MetaOpts_S3.ExpDesign.X = Xj;
            MetaOpts_S3.ExpDesign.Y = S3_raw;
            MetaOpts_S3.Name = sprintf('PCE_S3_Node_%d_%d', i,j);
            MetaOpts_S3.Display = 'quiet';
            pceModel_S3 = uq_createModel(MetaOpts_S3);
            fprintf('Node %d: normal signal; using PCE.\n', j);
            fprintf('Node %d: S3 relative range: %.4f\n', j, S3_relative_range);
        end
        
        % Create the PCE model for Y and S1.
        pceModel = uq_createModel(MetaOpts_PCE);
        
        % Report the selected PCE degree.
        try
            % Access PCE properties directly to avoid deeply nested checks.
            pce_data = pceModel.PCE;
            basis_data = pce_data.Basis;

            % Extract polynomial degree information.
            if isfield(basis_data, 'Degree')
                actual_degree = basis_data.Degree;
            elseif isfield(basis_data, 'MaxCompDeg')
                % Take the maximum of the per-variable degrees in MaxCompDeg.
                actual_degree = max(basis_data.MaxCompDeg);
            else
                % Calculate the maximum total degree from Indices.
                if isfield(basis_data, 'Indices')
                    indices = basis_data.Indices;
                    if isnumeric(indices)
                        actual_degree = max(sum(indices, 2));
                    else
                        actual_degree = NaN;
                    end
                else
                    actual_degree = NaN;
            end
        end
        
            % Extract the number of terms.
            if isfield(basis_data, 'Indices')
                indices = basis_data.Indices;
                if isnumeric(indices)
                    n_terms = size(indices, 1);
                else
                    n_terms = length(indices);
                end
            else
                n_terms = NaN;
            end

            % Store results.
            if ~isnan(actual_degree) && ~isnan(n_terms)
                pce_degrees(j) = actual_degree;
                pce_terms(j) = n_terms;
                if j <= 5 || mod(j,20) == 0
                    fprintf('Node %d: selected PCE degree=%d, terms=%d\n', j, actual_degree, n_terms);
                end
            else
                fprintf('Node %d: degree=%.1f, terms=%.1f (including NaN)\n', j, actual_degree, n_terms);
            end

        catch ME
            fprintf('Node %d: could not extract PCE degree information - %s\n', j, ME.message);
            pce_degrees(j) = NaN;
            pce_terms(j) = NaN;
        end
        
        % Create the adaptive model structure.
        hybridModels{j} = struct();
        hybridModels{j}.PCE = pceModel;          % PCE model for Y and S1.
        hybridModels{j}.PCE_S3 = pceModel_S3;    % S3 model (may be empty).
        hybridModels{j}.usePCK_S3 = false;
        hybridModels{j}.PCK_S3 = [];
        hybridModels{j}.model_type = 'Adaptive';  % Adaptive model type.
        hybridModels{j}.Y_S1_method = 'PCE';     % PCE for Y and S1.
        hybridModels{j}.S3_method = S3_model_method;  % Adaptive method for S3.
        hybridModels{j}.S3_relative_range = S3_relative_range;  % Record signal strength.
        
        % Store method-specific model parameters.
        if strcmp(S3_model_method, 'constant')
            hybridModels{j}.S3_constant = S3_pred_value;
        elseif strcmp(S3_model_method, 'linear')
            hybridModels{j}.S3_coeffs = S3_coeffs;
            hybridModels{j}.S3_X_design = X_design;
        end
        % % Retain S3 normalization fields for compatibility.
        % hybridModels{j}.S3_mean = S3_mean;
        % hybridModels{j}.S3_std = S3_std;
        
        % Retain the original variable reference for compatibility.
        mySurrogateModel{j} = pceModel;
        
        % Compute several error metrics.
        % 1. LOO error with model-specific evaluation.
        e_loo_pce = extract_LOO_vector(pceModel);    % LOO errors for Y and S1.
        
        % Evaluate the S3 error according to the model type.
        if strcmp(S3_model_method, 'constant')
            % Constant-model error: relative standard deviation.
            e_loo_s3 = std(S3_raw) / abs(mean(S3_raw));
            
        elseif strcmp(S3_model_method, 'linear')
            % Compute the linear-regression LOO error explicitly.
            n_samples = length(S3_raw);
            loo_errors = zeros(n_samples, 1);
            
            for k = 1:n_samples
                % Leave out sample k.
                train_idx = setdiff(1:n_samples, k);
                X_train = X_design(train_idx, :);
                Y_train = S3_raw(train_idx);
                X_test = X_design(k, :);
                Y_test = S3_raw(k);
                
                % Train a linear model.
                coeffs_loo = X_train \ Y_train;
                Y_pred_loo = X_test * coeffs_loo;
                
                % Calculate relative error.
                loo_errors(k) = abs(Y_test - Y_pred_loo) / abs(Y_test);
            end
            e_loo_s3 = mean(loo_errors);
            
        else
            % PCE LOO error.
            e_loo_s3_raw = extract_LOO_vector(pceModel_S3);
            if ~isempty(e_loo_s3_raw)
                e_loo_s3 = e_loo_s3_raw(1);
            else
                e_loo_s3 = NaN;
                end
            end
            
        % Combine the LOO errors.
        if ~isempty(e_loo_pce) && ~isnan(e_loo_s3)
            e_loo = [e_loo_pce, e_loo_s3];
        elseif ~isempty(e_loo_pce)
            e_loo = [e_loo_pce, NaN];
        elseif ~isnan(e_loo_s3)
            e_loo = [NaN, NaN, e_loo_s3];
        else
            e_loo = [];
        end
        
        % 2. Relative mean squared error.
        try
            % Predict Y and S1.
            Y_pred_pce = uq_evalModel(pceModel, Xj);      % Y and S1.
            
            % Predict S3 using its selected model.
            if strcmp(S3_model_method, 'constant')
                S3_pred = repmat(S3_pred_value, size(Xj,1), 1);
            elseif strcmp(S3_model_method, 'linear')
                S3_pred = X_design * S3_coeffs;
            else
                S3_pred = uq_evalModel(pceModel_S3, Xj);
            end
            
            % Assemble the complete prediction matrix.
            Y_pred_full = [Y_pred_pce, S3_pred];
            
            % Assemble the reference response matrix.
            Y_true = [Y_raw, S1_raw, S3_raw];
            
            if size(Y_pred_full, 2) == size(Y_true, 2)
                % Calculate relative MSE: MSE / var(Y_true).
                mse_abs = mean((Y_true - Y_pred_full).^2, 1);
                y_var = var(Y_true, 0, 1);  % Calculate the variance of each column.
                e_mse = mse_abs ./ max(y_var, eps);  % Use eps as a lower bound to avoid division by zero.
            else
                e_mse = [];
            end
        catch
            e_mse = [];
        end
        
        % 3. Five-fold CV is skipped for the hybrid model.
        e_cv5 = [];  % Skip CV5 because its hybrid-model implementation is complex.
        
        % Store LOO errors.
        if ~isempty(e_loo)
            nout = numel(e_loo);
            if size(loo_all_per_node,2) < nout
                loo_all_per_node(:, end+1:nout) = NaN;
            end
            loo_all_per_node(j,1:nout) = e_loo;
            loo_max_per_node(j) = max(e_loo, [], 'omitnan');
        else
            loo_max_per_node(j) = NaN;
        end
        
        % Store MSE errors.
        if ~isempty(e_mse)
            nout = numel(e_mse);
            if size(mse_all_per_node,2) < nout
                mse_all_per_node(:, end+1:nout) = NaN;
            end
            mse_all_per_node(j,1:nout) = e_mse;
            mse_max_per_node(j) = max(e_mse, [], 'omitnan');
        else
            mse_max_per_node(j) = NaN;
        end
        
        % Store five-fold CV errors.
        if ~isempty(e_cv5)
            nout = numel(e_cv5);
            if size(cv5_all_per_node,2) < nout
                cv5_all_per_node(:, end+1:nout) = NaN;
            end
            cv5_all_per_node(j,1:nout) = e_cv5;
            cv5_max_per_node(j) = max(e_cv5, [], 'omitnan');
        else
            cv5_max_per_node(j) = NaN;
        end
        
        % Store the corrected error metrics in the hybrid model.
        if ~isempty(e_loo) && length(e_loo) >= 3
            hybridModels{j}.corrected_loo_errors = e_loo;  % Store corrected LOO metrics (S3 uses the model-specific metric).
        end
        if ~isempty(e_mse) && length(e_mse) >= 3
            hybridModels{j}.relative_mse_errors = e_mse;   % Store relative MSE.
        end
        
        % Report adaptive modeling results for initial nodes and periodically.
        if j <= 5 || mod(j,1000) == 0
            fprintf('Node %d: adaptive modeling - Y&S1(PCE), S3(%s)\n', j, S3_model_method);
            fprintf('  S3 relative variation: %.4f%%, selected model: %s\n', S3_relative_range*100, S3_model_method);
            if ~isempty(e_loo) && length(e_loo) >= 3
                fprintf('  LOO error: Y=%.2e, S1=%.2e, S3=%.2e\n', e_loo(1), e_loo(2), e_loo(3));
            end
            if ~isempty(e_mse) && length(e_mse) >= 3
                fprintf('  Relative MSE: Y=%.4f, S1=%.4f, S3=%.4f\n', e_mse(1), e_mse(2), e_mse(3));
                if e_loo(3) < 0.05  % Judge quality using LOO rather than MSE.
                    fprintf('  -> S3 %s fit is excellent.\n', S3_model_method);
                elseif e_loo(3) < 0.1
                    fprintf('  -> S3 %s fit is good.\n', S3_model_method);
                elseif e_loo(3) > 0.2
                    fprintf('  -> S3 %s fit needs improvement.\n', S3_model_method);
                end
            end
        end
        
    
        %     if j <= 5 || mod(j,20) == 0
        %         current_error = mySurrogateModel{j}.Error.LOO;
        %         fprintf('Node %d model quality: LOO error = %.2e\n', j, max(current_error));
        %     end
    end
         % Compare summary statistics across error metrics.
         fprintf('\n=== Adaptive model error comparison (sample size=%d) ===\n', NN);
         fprintf('Model configuration: joint PCE for Y&S1 (degrees 6-8); adaptive constant/linear/PCE for S3.\n');
        
        % LOO error statistics.
        valid_loo = loo_max_per_node(isfinite(loo_max_per_node));
        if ~isempty(valid_loo)
            fprintf('LOO error statistics:\n');
            fprintf('  Mean: %.2e\n', mean(valid_loo));
            fprintf('  Median: %.2e\n', median(valid_loo));
            fprintf('  Standard deviation: %.2e\n', std(valid_loo));
            fprintf('  Maximum: %.2e\n', max(valid_loo));
            fprintf('  Minimum: %.2e\n', min(valid_loo));
            fprintf('  99th percentile: %.2e\n', prctile(valid_loo, 99));
        end
        
        % Relative MSE statistics.
        valid_mse = mse_max_per_node(isfinite(mse_max_per_node));
        if ~isempty(valid_mse)
            fprintf('Relative MSE statistics (MSE/Var):\n');
            fprintf('  Mean: %.3f\n', mean(valid_mse));
            fprintf('  Median: %.3f\n', median(valid_mse));
            fprintf('  Standard deviation: %.3f\n', std(valid_mse));
            fprintf('  Maximum: %.3f\n', max(valid_mse));
            fprintf('  Minimum: %.3f\n', min(valid_mse));
            fprintf('  99th percentile: %.3f\n', prctile(valid_mse, 99));
        end
        
        % Relative five-fold CV statistics.
        valid_cv5 = cv5_max_per_node(isfinite(cv5_max_per_node));
        if ~isempty(valid_cv5)
            fprintf('Relative five-fold CV statistics (MSE/Var):\n');
            fprintf('  Mean: %.3f\n', mean(valid_cv5));
            fprintf('  Median: %.3f\n', median(valid_cv5));
            fprintf('  Standard deviation: %.3f\n', std(valid_cv5));
            fprintf('  Maximum: %.3f\n', max(valid_cv5));
            fprintf('  Minimum: %.3f\n', min(valid_cv5));
            fprintf('  99th percentile: %.3f\n', prctile(valid_cv5, 99));
        end
        
        % Compare error metrics.
        fprintf('\nError metric comparison:\n');
        if ~isempty(valid_mse) && ~isempty(valid_cv5)
            mse_cv5_ratio = mean(valid_mse) / mean(valid_cv5);
            fprintf('  Relative MSE / relative CV5 ratio: %.3f\n', mse_cv5_ratio);
            if abs(mse_cv5_ratio - 1) < 0.1
                fprintf('  -> Relative MSE and CV5 are similar; model performance is consistent.\n');
            elseif mse_cv5_ratio > 1.2
                fprintf('  -> Relative MSE exceeds CV5; possible overfitting.\n');
            elseif mse_cv5_ratio < 0.8
                fprintf('  -> Relative MSE is below CV5; training performance is better.\n');
                end
            end
            
        if ~isempty(valid_mse)
            fprintf('  Interpretation of relative MSE:\n');
            fprintf('    - Near 0: accurate predictions.\n');
            fprintf('    - Near 1: prediction error equals the data variance.\n');
            fprintf('    - Above 1: worse than the constant-mean baseline.\n');
        end
        
        fprintf('===============================\n');
        
         % Detailed analysis of the adaptive S3 model.
         fprintf('\n=== Adaptive S3 model analysis ===\n');
         valid_s3_loo = loo_all_per_node(isfinite(loo_all_per_node(:,3)), 3);
         valid_s3_mse = mse_all_per_node(isfinite(mse_all_per_node(:,3)), 3);
         
         % Count the use of each modeling method.
         constant_count = 0; linear_count = 0; pce_count = 0;
         for k = 1:length(hybridModels)
             if isfield(hybridModels{k}, 'S3_method')
                 switch hybridModels{k}.S3_method
                     case 'constant', constant_count = constant_count + 1;
                     case 'linear', linear_count = linear_count + 1;
                     case 'PCE', pce_count = pce_count + 1;
                 end
                end
            end
            
         fprintf('S3 model distribution:\n');
         fprintf('  Constant model: %d nodes (%.1f%%) - very weak signal\n', constant_count, 100*constant_count/Nnodes);
         fprintf('  Linear model: %d nodes (%.1f%%) - weak signal\n', linear_count, 100*linear_count/Nnodes);
         fprintf('  PCE model: %d nodes (%.1f%%) - normal signal\n', pce_count, 100*pce_count/Nnodes);
         
         if ~isempty(valid_s3_loo)
             excellent_loo_count = sum(valid_s3_loo < 0.05);
             good_loo_count = sum(valid_s3_loo < 0.1);
             acceptable_loo_count = sum(valid_s3_loo < 0.2);
             fprintf('\nS3 performance based on the reported LOO metric:\n');
             fprintf('  Excellent-fit nodes (LOO<5%%): %d/%d (%.1f%%)\n', ...
                 excellent_loo_count, length(valid_s3_loo), 100*excellent_loo_count/length(valid_s3_loo));
             fprintf('  Good-fit nodes (LOO<10%%): %d/%d (%.1f%%)\n', ...
                 good_loo_count, length(valid_s3_loo), 100*good_loo_count/length(valid_s3_loo));
             fprintf('  Acceptable-fit nodes (LOO<20%%): %d/%d (%.1f%%)\n', ...
                 acceptable_loo_count, length(valid_s3_loo), 100*acceptable_loo_count/length(valid_s3_loo));
             fprintf('  LOO statistics: mean=%.4f, median=%.4f, maximum=%.4f\n', ...
                 mean(valid_s3_loo), median(valid_s3_loo), max(valid_s3_loo));
         end
         
         fprintf('\nS3 modeling selects a method according to signal strength; inspect the reported error metrics.\n');
         fprintf('============================\n\n');
    %% Hybrid-model LOO statistics and visualization.
    fprintf('\n=== Hybrid-model LOO error visualization ===\n');
    % 1) Collect nodewise errors: PCE LOO for Y/S1 and the selected S3 metric.
    [loo_hybrid_per_node, loo_hybrid_max] = collect_hybrid_LOO(hybridModels, Nnodes);
    
    % 2) Plot in the original style with HYBRID titles and filenames.
    labels    = {'Displacement (Y)', 'Principal stress 1 (S1)', 'Principal stress 3 (S3)'};
    var_names = {'Y','S1','S3'};
    plot_loo_suite(xn, yn, zn, loo_hybrid_max, loo_hybrid_per_node, r, labels, var_names, 'Hybrid');
    
    % 3) Export hybrid-model errors to CSV.
    T_h = table(xn(:), yn(:), zn(:), loo_hybrid_max(:), ...
                loo_hybrid_per_node(:,1), loo_hybrid_per_node(:,2), loo_hybrid_per_node(:,3), ...
                'VariableNames', {'X','Y','Z','LOO_max','LOO_Y','LOO_S1','LOO_S3'});
    csv_h = sprintf('HYBRID_LOO_errors_r%.2f.csv', r);
    writetable(T_h, csv_h);
    fprintf('Hybrid-model LOO error data exported: %s\n', csv_h);
    
    % 4) Summarize hybrid-model errors.
    valid_nodes_h = isfinite(loo_hybrid_max);
    n_valid_h = sum(valid_nodes_h);
    if n_valid_h > 0
        fprintf('\n=== Hybrid-model error summary ===\n');
        fprintf('Valid nodes: %d/%d (%.1f%%)\n', n_valid_h, Nnodes, 100*n_valid_h/Nnodes);
        fprintf('Maximum LOO error statistics:\n');
        fprintf('  Mean: %.2e\n', mean(loo_hybrid_max(valid_nodes_h)));
        fprintf('  Median: %.2e\n', median(loo_hybrid_max(valid_nodes_h)));
        fprintf('  Standard deviation: %.2e\n', std(loo_hybrid_max(valid_nodes_h)));
        fprintf('  Minimum: %.2e\n', min(loo_hybrid_max(valid_nodes_h)));
        fprintf('  Maximum: %.2e\n', max(loo_hybrid_max(valid_nodes_h)));
    
        dim_labels = {'Displacement (Y)', 'Principal stress 1 (S1)', 'Principal stress 3 (S3)'};
        for d = 1:min(3, size(loo_hybrid_per_node,2))
            valid_d = isfinite(loo_hybrid_per_node(:,d));
            if any(valid_d)
                fprintf('\n%s hybrid-model error statistics:\n', dim_labels{d});
                fprintf('  Mean: %.2e\n', mean(loo_hybrid_per_node(valid_d,d)));
                fprintf('  Median: %.2e\n', median(loo_hybrid_per_node(valid_d,d)));
                fprintf('  95th percentile: %.2e\n', prctile(loo_hybrid_per_node(valid_d,d), 95));
                fprintf('  High-error nodes (>1e-1): %d\n', sum(loo_hybrid_per_node(valid_d,d) > 1e-1));
            end
        end
    
        radial_dist_h = sqrt(xn(valid_nodes_h).^2 + yn(valid_nodes_h).^2);
        depth_vals_h  = zn(valid_nodes_h);
        [r_corr_h, r_p_h] = corr(radial_dist_h', loo_hybrid_max(valid_nodes_h), 'Type', 'Spearman');
        [z_corr_h, z_p_h] = corr(depth_vals_h',  loo_hybrid_max(valid_nodes_h), 'Type', 'Spearman');
        fprintf('\nSpatial distribution of hybrid-model errors:\n');
        fprintf('  Correlation with radial distance: %.3f (p=%.3f)\n', r_corr_h, r_p_h);
        fprintf('  Correlation with depth: %.3f (p=%.3f)\n', z_corr_h, z_p_h);
        high_error_ratio_h = sum(loo_hybrid_max(valid_nodes_h) > 1e-2) / n_valid_h;
        fprintf('\nSuggestions for improving the hybrid surrogate:\n');
        if high_error_ratio_h > 0.1
            fprintf('  High-error node fraction is large (%.1f%%); consider:\n', 100*high_error_ratio_h);
            fprintf('    - More training samples\n    - Higher PCE degrees or adjusted truncation\n    - Local models near hotspots\n');
            fprintf('    - PCK or local Kriging for S1/Y to reduce error\n');
        else
            fprintf('  Most nodes have acceptable error (%.1f%% low-error nodes).\n', 100*(1-high_error_ratio_h));
        end
    end
    %raise Exception("Stop the calculation")
    %% 3.2 MCS Analysis with field-level convergence tracking
    fprintf('\n=== Running MCS Analysis ===\n');
    LimitStateOpts.mHandle = @(X) limit_state_MohrCoulomb_hybrid(X, hybridModels);
    LimitStateOpts.isVectorized = true;
    LimitStateOpts.Name = sprintf('MohrCoulombLSF_%d', i);
    myLSF = uq_createModel(LimitStateOpts);

    MCSopt = struct();
    MCSopt.Type = 'Reliability';
    MCSopt.Method = 'MCS';
    MCSopt.Model = myLSF;
    MCSopt.Input = myPriorDist;
    MCSopt.Simulation.BatchSize = 1e4;
    MCSopt.Simulation.MaxSampleSize = 1e5;
    MCSopt.LimitState.Threshold = -1e-6;
    MCSopt.LimitState.CompOp = '>=';
    %MCSopt.LimitState.CompOp = '<=';
    tMCS = tic;   % << start MCS timer
    myMCS = uq_createAnalysis(MCSopt);
    resMCS = myMCS.Results;
    if isfield(resMCS, 'Pf') && ~isempty(resMCS.Pf)
        Pf_values(i,:) = resMCS.Pf;
        fprintf('MCS Pf results obtained successfully\n');
    else
        error('MCS analysis failed to produce Pf results.');
    end
    %actual_mcs_samples = resMCS.Simulation.BatchSize;
    actual_mcs_samples = max(resMCS.ModelEvaluations(:));
    %fprintf('MCS completed with %d samples\n', actual_mcs_samples);
    mcs_time_sec = toc(tMCS);   % << stop MCS timer
    fprintf('MCS wall time: %.2f s\n', mcs_time_sec);

    %% 3.3 SS+IS+GMM Analysis
    fprintf('\n=== Running SS+IS+GMM Analysis ===\n');
    tSSIS_total = tic;  % << start total SS+IS+GMM timer
    PriorOpts_T = PriorOpts; 
    myPriorDist_Transform = uq_createInput(PriorOpts_T);
    % Safe transform wrappers supporting different UQLab signatures.
    % Identify nonconstant dimensions and the effective dimension d_eff.
    d = numel(myPriorDist_Transform.Marginals);
    isConst = arrayfun(@(m) strcmpi(m.Type,'Constant'), myPriorDist_Transform.Marginals);
    d_eff   = sum(~isConst);
    StdInput = createStdNormalInput(d);              % Standard-normal input with Parameters fields in its marginals.
    T_X2U = @(X) X2U_safe(X, myPriorDist_Transform, StdInput);
    T_U2X = @(U) U2X_safe(U, myPriorDist_Transform, StdInput);
    d_full      = d;                 % Full input dimension.
    T_U2X_full  = T_U2X;             % Full-dimensional U-to-X transform.

    % Standard-normal log density.
    logphi = @(U) -0.5*sum(U.^2,2) - 0.5*size(U,2)*log(2*pi);

    % Scalar system limit state: k-of-n aggregation (kth largest).
    G_matrix = @(X) limit_state_MohrCoulomb_hybrid(X, hybridModels);    % N×Nnodes
    %k_for_SS = max(3, round(0.01*Nnodes));       % Approximately 1% of nodes, at least 3.
    % k_for_SS = 2;
    % h_sys    = @(X) kth_largest_rows(G_matrix(X), k_for_SS);         % N×1
    h_sys = @(X) quantile(G_matrix(X), 0.9, 2);
    
    % 3.3.1 Pilot sampling to estimate thr0 with exceedance probability p0_target.
    %b_final   = 1e6;      % Failure threshold aligned with MCS/IS.
    b_final   = -1e-6;
    %p0_target = 0.10; % Target conditional probability per level (often 0.1-0.3).
    p0_target = 0.20;
    Npilot    = 3000;                          % Pilot sample size.
    Xpilot    = uq_getSample(myPriorDist_Transform, Npilot, 'LHS');
    hvals     = h_sys(Xpilot);
    thr0      = quantile(hvals, 1 - p0_target);
    %thr0      = quantile(hvals, p0_target);
    fprintf('→ Testing scalar system LSF for SS...\n');
    fprintf('  ✓ h_sys range: [%.2e, %.2e], thr0=%.2e (p0≈%.2f)\n', ...
            min(hvals), max(hvals), thr0, p0_target);

    % 3.3.2 UQLab subset simulation with the MCS/IS final threshold.
    SSopt = struct();
    SSopt.Type = 'Reliability';
    SSopt.Method = 'Subset';
    SSopt.Model = uq_createModel(struct('mHandle',h_sys,'isVectorized',true, ...
                                        'Name',sprintf('MohrCoulomb_SYSTEM_h_%d', i)));
    SSopt.Input = myPriorDist_Transform;
    % SSopt.LimitState.Threshold = thr0;      % Pilot quantile threshold.
    SSopt.LimitState.Threshold = b_final;      % Final failure threshold; intermediate levels are selected by subset simulation.
    SSopt.LimitState.CompOp    = '>=';
    SSopt.Subset.p0            = p0_target;
    %SSopt.Simulation.BatchSize = 20000;      % Sufficient samples per level.
    %SSopt.Simulation.BatchSize = 10000;      % Sufficient samples per level.
    SSopt.Simulation.BatchSize = 5000;      % Samples per subset level.
    SSopt.Simulation.MaxSampleSize = 2e6;
    SSopt.Subset.MaxSubsets    = 20;
    SSopt.Subset.Proposal.Type = 'Gaussian';
    SSopt.Subset.Proposal.Parameters = 0.8; % Moderate proposal step to support mixing.
    SSopt.Subset.Proposal.Adaptive = true;
    SSopt.Subset.Proposal.AdaptivePeriod = 100;
    SSopt.Simulation.Alpha     = 0.05;

    fprintf('Running Subset Simulation...\n');
    %tSS = tic;  % << start Subset timer
    mySS = uq_createAnalysis(SSopt);
    S    = mySS.Results;
    % ss_time_sec = toc(tSS);  % << stop Subset timer
    % fprintf('Subset Simulation wall time: %.2f s\n', ss_time_sec);
    %% Construct U_hot strictly from subset-simulation history, without fallback.
    % 1) Read U from the subset-simulation history.
    assert(isfield(S,'History') && isfield(S.History,'U') && ~isempty(S.History.U), ...
    'SS.Results.History.U is missing or empty; check that history output is enabled.');
    U_cells  = S.History.U;
    m_levels = numel(U_cells);
    assert(m_levels >= 1, 'SS History.U is empty; no subsets were generated.');
    % 2) Use the final level and, if available, the penultimate level.
    pick_levels = m_levels;
    if m_levels >= 2
        pick_levels = [m_levels-1, m_levels];
    end
    U_hot_list = {};
    for kk = pick_levels
        U_k = U_cells{kk};
        % 3) Check full/effective dimensions and select the appropriate U-to-X transform.
        if size(U_k,2) == d_full
            toX   = @(Ufull) T_U2X_full(Ufull);
            toEff = @(Ufull) Ufull(:, ~isConst);
        elseif size(U_k,2) == d_eff
            toX   = @(Ueff)  Ueff_to_X_helper(Ueff, isConst, T_U2X_full);
            toEff = @(Ueff)  Ueff;
        else
            error('Level %d: History.U has %d columns, incompatible with d_full=%d/d_eff=%d.', ...
                  kk, size(U_k,2), d_full, d_eff);
        end
    
        % 4) Evaluate h_sys and select a band above b_final.
        X_k = toX(U_k);
        h_k = h_sys(X_k);
    
        % Check polarity: a substantial fraction of the final level should exceed b_final.
        if kk == m_levels
            prop = mean(h_k >= b_final);
            if prop < 0.05
                warning(['Only %.1f%% of final-level samples exceed the threshold (b_final=%.2e).' ...
                         'This usually indicates inconsistent threshold, polarity, or CompOp settings.'], 100*prop, b_final);
            end
        end
        % Distance above the threshold for CompOp = >=.
        dist = max(h_k - b_final, 0);
        % Use the selected positive-distance quantile as the near-boundary bandwidth.
        if ~any(dist > 0)
            error('Level %d: no samples exceed the threshold; cannot form the near-boundary set.', kk);
        end
        %Delta = quantile(dist(dist>0), 0.50);
        Delta = quantile(dist(dist>0), 0.60);
        mask = (h_k >= b_final) & (h_k <= b_final + Delta);
        % Fail explicitly if the current level has no near-boundary samples.
        n_sel = sum(mask);
        if n_sel == 0
            error('Level %d: empty near-boundary set (Delta=%.3e). Increase level size/p0 or adjust the proposal step.', kk, Delta);
        end
    
        U_hot_list{end+1} = toEff(U_k(mask, :)); 
        fprintf('SS L%d → near-boundary: %d/%d (Δ≈%.2e)\n', kk, n_sel, size(U_k,1), Delta);
    end
    
    % 5) Pool and subsample to reduce dependence.
    U_hot = vertcat(U_hot_list{:});
    max_hot = 20000;                       % Adjust as needed.
    if size(U_hot,1) > max_hot
        U_hot = U_hot(randperm(size(U_hot,1), max_hot), :);
    end
    
    % 6) Jitter nearly constant columns to avoid singular GMM covariance matrices.
    col_std = std(U_hot, 0, 1);
    if any(col_std < 1e-10)
        U_hot(:, col_std < 1e-10) = U_hot(:, col_std < 1e-10) + ...
                                    1e-6 * randn(size(U_hot,1), sum(col_std < 1e-10));
    end

    % Fit GMMs and scan BIC/DeltaBIC over candidate component counts.
    Kmax = 6;
    opts = statset('MaxIter',2000,'TolFun',1e-7);
    
    
    BICs  = nan(Kmax,1);
    GMMS  = cell(Kmax,1);
    ok    = false(Kmax,1);
    
    for Kc = 1:Kmax
        try
            gmm_c = fitgmdist(U_hot, Kc, ...
                'RegularizationValue', 1e-3, ...
                'Replicates', 10, ...
                'Options', opts, ...
                'Start', 'plus');   % k-means++ initialization.
            BICs(Kc) = gmm_c.BIC;
            GMMS{Kc} = gmm_c;
            ok(Kc)   = true;
        catch ME
            % Record failed or nonconverged fits as NaN and continue.
            BICs(Kc) = NaN;
            GMMS{Kc} = [];
            ok(Kc)   = false;
            fprintf('[WARN] Fitting failed for K=%d: %s\n', Kc, ME.message);
        end
    end
    
    % Select the best K while ignoring NaN entries.
    if ~any(ok)
        error('All GMM fits failed; check U_hot and regularization/initialization settings.');
    end
    [bestBIC, idxBest] = min(BICs(ok));
    K_candidates = find(ok);
    K = K_candidates(idxBest);
    bestGMM = GMMS{K};
    
    % Compute DeltaBIC and approximate BIC-based weights.
    detBIC = BICs - bestBIC;             % DeltaBIC.
    w = exp(-0.5 * detBIC);              % Conventional exponential information-criterion weights.
    w(~isfinite(w)) = 0;
    sw = sum(w);
    if sw>0, w = w/sw; end
    
    % Print the results table.
    fprintf('\n=== GMM BIC scan (same data and settings) ===\n');
    fprintf('  K      BIC         detBIC     weight    status\n');
    for Kc = 1:Kmax
        if ok(Kc)
            fprintf('%3d  %10.1f   %7.2f    %7.3f    OK\n', ...
                Kc, BICs(Kc), detBIC(Kc), w(Kc));
        else
            fprintf('%3d  %10s   %7s    %7s    FAIL\n', ...
                Kc, 'NaN', 'NaN', 'NaN');
        end
    end
    fprintf('=> Selected K = %d (minimum BIC).\n', K);
    % Set the final GMM and component count.
    gmm = bestGMM;
    % Save CSV.
    try
        T = table((1:Kmax)', BICs, detBIC, w, ok, ...
            'VariableNames', {'K','BIC','detBIC','Weight','OK'});
        csvname = sprintf('GMM_BIC_scan_r%.2f.csv', r);
        writetable(T, csvname);
        fprintf('DeltaBIC table saved to: %s\n', csvname);
    catch
        % Ignore optional export failures on an unwritable filesystem.
    end
    
    % Optional diagnostic plot.
    try
        hfig_bic = figure('Name','BIC / detBIC vs K','Position',[120 120 900 350]);
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    
        nexttile; 
        plot(1:Kmax, BICs, 'o-','LineWidth',1.5); grid on;
        xlabel('K'); ylabel('BIC'); title('BIC vs K');
    
        nexttile;
        plot(1:Kmax, detBIC, 'o-','LineWidth',1.5); hold on; grid on;
        yline(0,'k:'); % Baseline for the best model.
        xlabel('K'); ylabel('\DeltaBIC (detBIC)'); title('\DeltaBIC vs K');
    
        saveas(hfig_bic, sprintf('GMM_BIC_scan_r%.2f.png', r), 'png');
    catch
    end
    % ==================================================

    % Prepare the covariance-inflated auxiliary proposal.
    %tau = 1.8;         % Standard-deviation inflation factor (often 1.5-2.0).
    %eta = 0.15;         % Mixture weight allocated to the wide GMM.
    tau = 2.0;               
    eta = 0.20;         
    Kc  = gmm.NComponents;
    Sigma_wide = gmm.Sigma;
    for kk = 1:Kc
        Sigma_wide(:,:,kk) = tau^2 * Sigma_wide(:,:,kk);
    end

    % Importance Sampling
    %Nis = 50000;  % Increase the sample size.
    %Nis = 20000;
    Nis = 10000;
    eps_def = 0.20;
    assert(eps_def + eta < 1, 'eps_def + eta must be < 1');
    d_eff   = size(U_hot,2);

    % Proposal source labels: 1=GMM, 2=wide GMM, 3=standard normal.
    pi_mix = [max(1 - eps_def - eta, 0), eta, eps_def];  % [GMM, GMM_wide, STD]
    src = randsample(3, Nis, true, pi_mix);
    U_is = zeros(Nis, d_eff);
    idx1 = (src==1); n1 = sum(idx1); if n1>0, U_is(idx1,:) = random(gmm, n1);                      end
    idx2 = (src==2); n2 = sum(idx2); if n2>0, U_is(idx2,:) = mix_sample_wide(gmm, Sigma_wide, n2); end
    idx3 = (src==3); n3 = sum(idx3); if n3>0, U_is(idx3,:) = randn(n3, d_eff);                     end
    fprintf('IS split: GMM=%d, GMM_wide=%d, STD=%d, total=%d\n', n1, n2, n3, Nis);
    % ---- U→X ----
    X_is = Ueff_to_X_helper(U_is, isConst, T_U2X_full);
    
    % Evaluate the limit-state function at every node.
    G_all = G_matrix(X_is);                    % [Nis × Nnodes]
    I_all = (G_all >= b_final);                % Failure indicator matrix.
    
    % Evaluate the mixture proposal density using a stable log-space calculation.
    % log q_gmm, log phi
    log_q_gmm = log( max(pdf(gmm, U_is), realmin) );   % [Nis × 1]
    log_phi    = -0.5*sum(U_is.^2,2) - 0.5*d_eff*log(2*pi);                         % [Nis x 1], standard-normal density in d_eff dimensions.
    log_q_wide = log( max(pdf_gmm_wide(U_is, gmm, Sigma_wide), realmin) );
    a = log(max(1 - eps_def - eta, realmin)) + log_q_gmm;   % (1-ε-η)·GMM
    b = log(eta)                              + log_q_wide; % η·GMM_wide
    c = log(eps_def)                          + log_phi;    % ε·φ
    m = max([a,b,c], [], 2);
    log_q_mix = m + log( exp(a - m) + exp(b - m) + exp(c - m) );
    
    % Self-normalized importance weights: rho = pi/q_mix; w = rho/sum(rho).
    rho = exp(log_phi - log_q_mix);            % Unnormalized weights.
    rho = max(rho, realmin);
    w = rho / sum(rho);            % SNIS weights summing to one.
    
    % Nodewise SNIS failure probability: Pf_j = sum_i w_i * I_ij.
    Pf_IS = (w' * I_all);                     % [Nnodes × 1]
    diff = I_all - Pf_IS;                     % [Nis × Nnodes]
    se2  = (w.^2)' * (diff.^2);                % [1 × Nnodes]
    se_pf = sqrt(se2);                         % [Nnodes × 1]
    z = 1.96;
    Pf_low = max(Pf_IS - z*se_pf, 0);
    Pf_up  = Pf_IS + z*se_pf;
    
    % Diagnostics: effective sample size and weight dispersion.
    ESS = 1 / sum(w.^2);                       % Standard ESS definition for normalized weights.
    w_desc = sort(w, 'ascend');
    k1 = ceil(0.01*Nis);
    k5 = ceil(0.05*Nis);
    fprintf('IS (3-mix: std=%.0f%%, wide=%.0f%%): Nis=%d, ESS≈%.0f (%.1f%%)\n', ...
            100*eps_def, 100*eta, Nis, ESS, 100*ESS/Nis);
    % Robust weight-concentration diagnostics.
    w = w(:);                        % Force a column vector to avoid dimension errors.
    sumw = sum(w);
    w_desc = sort(w, 'descend');     % [N × 1]
    k1 = max(1, ceil(0.01 * numel(w_desc)));
    k5 = max(1, ceil(0.05 * numel(w_desc)));
    top1_share = sum(w_desc(1:k1));
    top5_share = sum(w_desc(1:k5));
    fprintf('sum(w)=%.6f, Top1%%=%.4f, Top5%%=%.4f, min(w)=%.3e, max(w)=%.3e\n', ...
            sumw, top1_share, top5_share, min(w), max(w));

    fprintf('Top1%% weight share=%.3f, Top5%%=%.3f\n', sum(w_desc(1:k1)), sum(w_desc(1:k5)));
    %% 3.4 Save Analysis Results
    % fprintf('\n=== All Analysis Completed ===\n');
    % fprintf('You can now run the visualization script without re-running analysis\n');
    fprintf('\n=== Saving Analysis Results ===\n');
    % MCS
    ssis_total_time_sec = toc(tSSIS_total);  % << stop total SS+IS+GMM timer
    fprintf('SS+IS+GMM total wall time: %.2f s\n', ssis_total_time_sec);

    MCS_results = struct();
    MCS_results.Pf = Pf_values(i,:);
    MCS_results.actual_samples = actual_mcs_samples;
    MCS_results.r = r;
    MCS_results.coordinates = coords;
    save(sprintf('MCS_results_r%.2f.mat', r), 'MCS_results');
    % 
    % % Save SS+IS+GMM results without the large simulation history.
    IS_results = struct();
    IS_results.Pf = Pf_IS;
    IS_results.CI = [Pf_low; Pf_up];
    IS_results.ESS = ESS;
    IS_results.Nis = Nis;
    IS_results.K = K;
    IS_results.r = r;
    IS_results.coordinates = coords;
    IS_results.U_hot = U_hot;
    % IS_results.gmm_mu      = gmm.mu;
    % IS_results.gmm_Sigma   = gmm.Sigma;
    % IS_results.gmm_weights = gmm.ComponentProportion;
    IS_results.gmm = gmm;
    save(sprintf('IS_results_r%.2f.mat', r), 'IS_results');

    % Comparison
    comparison_data = struct();
    comparison_data.Pf_MCS = Pf_values(i,:);
    comparison_data.Pf_IS  = Pf_IS;
    comparison_data.Diff   = Pf_IS - Pf_values(i,:);
    comparison_data.r      = r;
    comparison_data.coordinates = coords;
    comparison_data.MCS_samples = actual_mcs_samples;
    comparison_data.IS_samples  = Nis;
    comparison_data.ESS         = ESS;
    comparison_data.GMM_K       = K;
    comparison_data.s3_normalization = s3_normalization; % Save S3 normalization fields.
    save(sprintf('comparison_data_r%.2f.mat', r), 'comparison_data');

    % CSV (IS results + CI)
    csv_data = zeros(Nnodes, 7);
    csv_data(:, 1) = xn(1,1:Nnodes)';
    csv_data(:, 2) = yn(1,1:Nnodes)';
    csv_data(:, 3) = zn(1,1:Nnodes)';
    csv_data(:, 4) = Pf_IS';
    csv_data(:, 5) = Pf_low';
    csv_data(:, 6) = Pf_up';
    csv_data(:, 7) = (1:Nnodes)';

    csv_headers = {'X_Coordinate','Y_Coordinate','Z_Coordinate', ...
                   'Failure_Probability','CI_Lower','CI_Upper','Node_ID'};
    csv_filename = sprintf('IS_results_r%.2f.csv', r);

    fid = fopen(csv_filename, 'w');
    assert(fid ~= -1, 'Cannot open CSV: %s', csv_filename);
    fprintf(fid, '%s', csv_headers{1});
    for col = 2:numel(csv_headers)
        fprintf(fid, ',%s', csv_headers{col});
    end
    fprintf(fid, '\n');
    for row = 1:size(csv_data, 1)
        fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%d\n', ...
            csv_data(row,1), csv_data(row,2), csv_data(row,3), ...
            csv_data(row,4), csv_data(row,5), csv_data(row,6), ...
            csv_data(row,7));
    end
    fclose(fid);
    fprintf('CSV file saved: %s\n', csv_filename);

    fprintf('Analysis completed for r = %.2f\n', r);
    fprintf('Files saved:\n  - MCS_results_r%.2f.mat\n  - IS_results_r%.2f.mat\n  - comparison_data_r%.2f.mat\n  - IS_results_r%.2f.csv\n', r, r, r, r);
end

fprintf('\n=== All Analysis Completed ===\n');
fprintf('You can now run the visualization script without re-running analysis\n');

end

function y = kth_largest_rows(G, k)
% Return the kth largest entry of each row in G (k-of-n aggregation).
    if isempty(G)
        y = [];
        return;
    end
    k = max(1, min(k, size(G,2)));
    Gs = sort(G, 2, 'descend');
    y  = Gs(:, k);
end

function StdInput = createStdNormalInput(d)
% Create d standard-normal marginals with Parameters fields for compatibility.
    StdOpts = struct();
    for kk = 1:d
        StdOpts.Marginals(kk).Type      = 'Gaussian';
        StdOpts.Marginals(kk).Moments   = [0 1];   % uq_createInput also populates Parameters.
    end
    StdInput = uq_createInput(StdOpts);
end

function U = X2U_safe(X, Input, StdInput)
% Safe X-to-U transform: try supported signatures and report a clear error.
    X = double(X);
    d = numel(Input.Marginals);
    assert(size(X,2) == d, 'X columns (%d) ≠ Input dim (%d).', size(X,2), d);
    if any(~isfinite(X), 'all'), error('X contains NaN/Inf.'); end

    % Signature using marginal and copula arguments.
    try
        U = uq_GeneralIsopTransform(X, Input.Marginals, Input.Copula, ...
                                       StdInput.Marginals, StdInput.Copula);
        return;
    catch
        % Alternate signature used by some UQLab versions.
        try
            U = uq_GeneralIsopTransform(X, Input, 'original', 'standardnormal');
            return;
        catch ME
            error('X→U transform failed (both signatures). Detail: %s', ME.message);
        end
    end
end

function X = U2X_safe(U, Input, StdInput)
% Safe U-to-X transform with alternate-signature fallback.
    U = double(U);
    d = numel(Input.Marginals);
    assert(size(U,2) == d, 'U columns (%d) ≠ Input dim (%d).', size(U,2), d);
    if any(~isfinite(U), 'all'), error('U contains NaN/Inf.'); end
    try
        X = uq_GeneralIsopTransform(U, StdInput.Marginals, StdInput.Copula, ...
                                       Input.Marginals, Input.Copula);
        return;
    catch
        try
            X = uq_GeneralIsopTransform(U, Input, 'standardnormal', 'original');
            return;
        catch ME
            error('U→X transform failed (both signatures). Detail: %s', ME.message);
        end
    end
end


% Helper functions.
function rho1 = lag1_autocorr(x)
    x = x(:) - mean(x);
    if numel(x) < 3 || std(x)==0, rho1 = 0; return; end
    rho1 = sum(x(1:end-1).*x(2:end)) / sum(x.^2);
    rho1 = max(min(rho1, 0.999), -0.999);
end

function ess = eff_size_lag1(N, rho1)
    ess = N * (1 - rho1) / (1 + rho1);
    ess = max(1, min(ess, N));
end

function X = Ueff_to_X_helper(Ueff, isConst, T_U2X_full)
    N = size(Ueff,1); d_full = numel(isConst);
    U_full = zeros(N, d_full);
    U_full(:, ~isConst) = Ueff;        % Fill nonconstant dimensions; leave constant dimensions at zero in U space.
    X = T_U2X_full(U_full);
end

function U = mix_sample_wide(gmm, Sigma_wide, n)
    comps = randsample(gmm.NComponents, n, true, gmm.ComponentProportion);
    U = zeros(n, size(gmm.mu,2));
    for i = 1:n
        k = comps(i);
        U(i,:) = mvnrnd(gmm.mu(k,:), Sigma_wide(:,:,k), 1);
    end
end

function p = pdf_gmm_wide(U, gmm, Sigma_wide)
    N = size(U,1); K = gmm.NComponents;
    log_comp = zeros(N, K);
    for k = 1:K
        log_comp(:,k) = log_mvnpdf(U, gmm.mu(k,:), Sigma_wide(:,:,k)) + log(gmm.ComponentProportion(k));
    end
    m = max(log_comp, [], 2);
    p = exp( m + log(sum(exp(log_comp - m), 2)) );
end

function l = log_mvnpdf(X, mu, Sigma)
    d = size(X,2);
    XC = bsxfun(@minus, X, mu);
    [R,p] = chol(Sigma); if p>0, error('Sigma not PD'); end
    xR = XC / R;  quad = sum(xR.^2, 2);
    c = d*log(2*pi) + 2*sum(log(diag(R)));
    l = -0.5*(c + quad);
end

function G = limit_state_MohrCoulomb_hybrid(X, hybridModelsCell)
% Return an N-by-Nnodes limit-state matrix, one node per column.
    Nnodes = numel(hybridModelsCell);
    N = size(X,1);
    G = zeros(N, Nnodes);

    for j = 1:Nnodes
        M = hybridModelsCell{j};

        % Evaluate Y/S1 with PCE and S3 with its selected adaptive model.
        YS1_PCE = uq_evalModel(M.PCE, X);    % [N x 2] predictions for Y and S1.
        Ycol  = YS1_PCE(:,1);
        S1col = YS1_PCE(:,2);

        % Predict S3 according to the adaptive model type.
        if isfield(M, 'S3_method')
            switch M.S3_method
                case 'constant'
                    S3col = repmat(M.S3_constant, N, 1);
                case 'linear'
                    X_design_eval = [ones(N,1), X];  % Add an intercept.
                    S3col = X_design_eval * M.S3_coeffs;
                case 'PCE'
                    S3col = uq_evalModel(M.PCE_S3, X);
                otherwise
                    error('Node %d: unknown S3 model type: %s', j, M.S3_method);
            end
        elseif isfield(M, 'PCE_S3') && ~isempty(M.PCE_S3)
            S3col = uq_evalModel(M.PCE_S3, X); % Fall back to PCE.
        elseif isfield(M, 'Kriging_S3') && ~isempty(M.Kriging_S3)
            S3col = uq_evalModel(M.Kriging_S3, X); % Fall back to Kriging.
        else
            error('Node %d: no S3 model is available.', j);
        end

        % === Mohr-Coulomb failure criterion ===
        % Extract input parameters from X
        cohesion = X(:,3);   % N×1
        phi = X(:,4);        % N×1
        
        % Calculate tan(phi) for Mohr-Coulomb formula
        tanphi = tan(phi/2 + 45*pi/180);  % Vectorized N-by-1 calculation.
        
        % Mohr-Coulomb failure criterion
        % g = -sigma3 - (-sigma1 * tan^2(phi) + 2*c*tan(phi))
        g_j = -S3col - (-S1col.*tanphi.^2 + 2*cohesion.*tanphi); % Original sign convention.
        %g_j = (-S1col.*tanphi.^2 + 2*cohesion.*tanphi)-(-S3col) ;
        % ==========================================================================
        G(:,j) = g_j;
    end
    
    % Debug output for first call
    if N <= 10  % Only for small samples to avoid spam
        gmax = max(G(:));
        gmin = min(G(:));
        fprintf('[DEBUG] limit_state_MohrCoulomb_hybrid: N=%d, gvals range = [%.4e, %.4e]\n', N, gmin, gmax);
        failure_count = sum(G >= -1e-6, 'all');
        fprintf('[DEBUG] Samples with g >= -1e-6 (failure): %d/%d (%.2f%%)\n', failure_count, numel(G), 100*failure_count/numel(G));
    end
end


function e = extract_LOO_vector(uqModel)
    e = [];
    if isempty(uqModel), return; end
    hasErr = (isobject(uqModel) && isprop(uqModel, 'Error')) || ...
             (isstruct(uqModel) && isfield(uqModel, 'Error'));
    if ~hasErr, return; end
    Err = uqModel.Error;
    if isstruct(Err) && numel(Err) >= 1 && all(arrayfun(@(s) isfield(s,'LOO'), Err))
        e = double([Err.LOO]); e = e(:).';
    elseif isstruct(Err) && isfield(Err,'LOO') && isnumeric(Err.LOO)
        e = double(Err.LOO(:)).';
    elseif isstruct(Err) && isfield(Err,'ModifiedLOO') && isnumeric(Err.ModifiedLOO)
        e = double(Err.ModifiedLOO(:)).';
    end
end

function [loo_per_node, loo_max] = collect_hybrid_LOO(hybridModels, Nnodes)
% Collect nodewise errors: PCE LOO for Y/S1 and the available S3 error metric.
    loo_per_node = nan(Nnodes, 3);
    loo_max      = nan(Nnodes, 1);

    for j = 1:Nnodes
        M = hybridModels{j};

        % Obtain Y/S1 errors from PCE.
        if isfield(M, 'PCE') && ~isempty(M.PCE)
            eP = extract_LOO_vector(M.PCE);
        if ~isempty(eP)
            if numel(eP) >= 1, loo_per_node(j,1) = eP(1); end
            if numel(eP) >= 2, loo_per_node(j,2) = eP(2); end
            end
        end

        % Prefer stored S3 metrics, followed by separate PCE_S3 and other models.
        if isfield(M, 'corrected_loo_errors') && length(M.corrected_loo_errors) >= 3
            % Use corrected errors calculated in the main loop.
            loo_per_node(j,3) = M.corrected_loo_errors(3);
        elseif isfield(M, 'PCE_S3') && ~isempty(M.PCE_S3)
            % Use LOO from the separate PCE_S3 model.
            eS3 = extract_LOO_vector(M.PCE_S3);
            if ~isempty(eS3)
                loo_per_node(j,3) = eS3(1);
            end
        elseif isfield(M, 'Kriging_S3') && ~isempty(M.Kriging_S3)
            % Fall back to Kriging LOO.
            try
                eS3 = extract_LOO_vector(M.Kriging_S3);
                if ~isempty(eS3)
                    loo_per_node(j,3) = eS3(1);
                end
            catch
                loo_per_node(j,3) = NaN;
            end
        elseif isfield(M,'usePCK_S3') && M.usePCK_S3 && isfield(M,'PCK_S3') && ~isempty(M.PCK_S3)
            eS3 = extract_LOO_vector(M.PCK_S3);
            if ~isempty(eS3)
                loo_per_node(j,3) = eS3(1);
            end
        else
            % Finally try the joint PCE if it has a third output.
            if isfield(M, 'PCE') && ~isempty(M.PCE)
                eP = extract_LOO_vector(M.PCE);
            if ~isempty(eP) && numel(eP) >= 3
                loo_per_node(j,3) = eP(3);
                 end
            end
        end

        % Maximum LOO error at each node.
        loo_max(j) = max(loo_per_node(j,:), [], 'omitnan');
    end
end

function plot_loo_suite(xn, yn, zn, loo_max, loo_all, r, labels, var_names, tag)
% Shared spatial/statistical LOO plots; append tag to titles and filenames.
    valid = isfinite(loo_max);
    if any(valid)
        % Combined spatial and statistical distribution plots.
        hfig1 = figure('Position', [100, 100, 1400, 900]);
        % 3D scatter plot.
        subplot(2,3,1);
        scatter3(xn(valid), yn(valid), zn(valid), 36, loo_max(valid), 'filled');
        try set(gca,'ColorScale','log'); end
        cb = colorbar; cb.Label.String = sprintf('Max %s LOO Error', tag);
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title(sprintf('%s maximum error distribution, r=%.2f', tag, r));
        grid on; axis equal; view(3);

        % Histogram.
        subplot(2,3,2);
        histogram(log10(loo_max(valid)), 30, 'FaceColor', [0.7 0.7 0.9]);
        xlabel(sprintf('log_{10}(Max %s LOO Error)', tag));
        ylabel('Number of nodes'); title('Histogram of maximum errors'); grid on;

        % Radial distance.
        subplot(2,3,3);
        radial_dist = sqrt(xn(valid).^2 + yn(valid).^2);
        scatter(radial_dist, loo_max(valid), 20, 'filled');
        set(gca, 'YScale', 'log');
        xlabel('Radial distance (m)'); ylabel(sprintf('Max %s LOO Error', tag));
        title('Maximum error versus radial distance'); grid on;

        % Depth.
        subplot(2,3,4);
        scatter(zn(valid), loo_max(valid), 20, 'filled');
        set(gca, 'YScale', 'log');
        xlabel('Depth Z (m)'); ylabel(sprintf('Max %s LOO Error', tag));
        title('Maximum error versus depth'); grid on;

        % Statistics.
        subplot(2,3,[5,6]);
        err_stats = [mean(loo_max(valid)), median(loo_max(valid)), ...
                     std(loo_max(valid)), min(loo_max(valid)), max(loo_max(valid))];
        bar(err_stats, 'FaceColor', [0.8 0.6 0.9]);
        set(gca,'XTickLabel',{'Mean','Median','Std','Min','Max'}, 'YScale','log');
        ylabel('Error'); title(sprintf('%s maximum error statistics', tag)); grid on; xtickangle(45);

        fname = sprintf('%s_LOO_error_field_comprehensive_r%.2f', upper(tag), r);
        saveas(hfig1, [fname '.png']); savefig(hfig1, [fname '.fig']);
        fprintf('%s maximum-error summary plot saved: %s\n', tag, fname);
    end

    % Per-output plots for Y, S1, and S3.
    for d = 1:min(3, size(loo_all,2))
        valid_d = isfinite(loo_all(:,d));
        if ~any(valid_d), continue; end

        hfig2 = figure('Position', [100+d*100, 100+d*100, 1400, 900]);
        % 3D
        subplot(2,4,[1,2,5,6]);
        scatter3(xn(valid_d), yn(valid_d), zn(valid_d), 45, loo_all(valid_d,d), 'filled', 'MarkerEdgeColor','k','LineWidth',0.1);
        try set(gca,'ColorScale','log'); end
        cb = colorbar; cb.Label.String = sprintf('%s %s LOO Error', labels{d}, tag); cb.Label.FontSize = 11;
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title(sprintf('%s %s spatial error distribution (r=%.2f)', labels{d}, tag, r), 'FontSize', 13);
        grid on; axis equal; view(-37.5, 30);

        % Histogram.
        subplot(2,4,3);
        histogram(log10(loo_all(valid_d,d)), 25, 'FaceColor', [0.6 0.8 0.9], 'EdgeColor', 'k');
        xlabel(sprintf('log_{10}(%s Error)', tag)); ylabel('Number of nodes'); title(sprintf('%s error distribution', labels{d})); grid on;

        % Radial distance.
        subplot(2,4,4);
        radial_dist = sqrt(xn(valid_d).^2 + yn(valid_d).^2);
        scatter(radial_dist, loo_all(valid_d,d), 20, 'filled'); set(gca,'YScale','log');
        xlabel('Radial distance (m)'); ylabel(sprintf('%s %s Error', labels{d}, tag)); title('Radial dependence'); grid on;

        % Depth.
        subplot(2,4,7);
        scatter(zn(valid_d), loo_all(valid_d,d), 20, 'filled'); set(gca,'YScale','log');
        xlabel('Depth Z (m)'); ylabel(sprintf('%s %s Error', labels{d}, tag)); title('Depth dependence'); grid on;

        % Statistics.
        subplot(2,4,8);
        err_vals  = loo_all(valid_d,d);
        err_stats = [mean(err_vals), median(err_vals), std(err_vals), min(err_vals), max(err_vals)];
        bar(err_stats, 'FaceColor', [0.8 0.6 0.9]);
        set(gca,'XTickLabel',{'Mean','Median','Std','Min','Max'}, 'YScale','log');
        ylabel(sprintf('%s error', tag)); title(sprintf('%s error statistics', labels{d})); grid on; xtickangle(45);

        fname = sprintf('%s_LOO_error_field_%s_comprehensive_r%.2f', upper(tag), var_names{d}, r);
        saveas(hfig2, [fname '.png']); savefig(hfig2, [fname '.fig']);
        fprintf('Saved %s-%s error analysis plot: %s\n', labels{d}, tag, fname);
    end

    % Combined comparison in the original style.
    hfig3 = figure('Position', [50, 50, 1600, 1000]);

    for d = 1:min(3, size(loo_all,2))
        valid_d = isfinite(loo_all(:,d));
        if ~any(valid_d), continue; end
        subplot(2,4,d);
        scatter3(xn(valid_d), yn(valid_d), zn(valid_d), 25, loo_all(valid_d,d), 'filled');
        try set(gca,'ColorScale','log'); end
        cb = colorbar; cb.Label.String = sprintf('%s LOO Error (%s)', labels{d}, tag); cb.Label.FontSize = 10;
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title(sprintf('%s error distribution', labels{d}), 'FontSize', 12); grid on; view(-37.5, 30);
    end

    subplot(2,4,4);
    valid_max = isfinite(loo_max);
    if any(valid_max)
        scatter3(xn(valid_max), yn(valid_max), zn(valid_max), 25, loo_max(valid_max), 'filled');
        try set(gca,'ColorScale','log'); end
        cb = colorbar; cb.Label.String = sprintf('Max LOO Error (%s)', tag); cb.Label.FontSize = 10;
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title('Maximum error distribution', 'FontSize', 12); grid on; view(-37.5, 30);
    end

    subplot(2,4,[5,6]);
    all_errors = []; error_labels = {};
    for d = 1:min(3, size(loo_all,2))
        valid_d = isfinite(loo_all(:,d));
        if any(valid_d)
            all_errors = [all_errors; loo_all(valid_d,d)];
            error_labels = [error_labels; repmat(labels(d), sum(valid_d), 1)];
        end
    end
    if ~isempty(all_errors)
        boxplot(log10(all_errors), error_labels);
        ylabel(sprintf('log_{10}(%s LOO Error)', tag));
        title('Y/S1/S3 error distribution comparison', 'FontSize', 12); grid on; xtickangle(45);
    end

    subplot(2,4,7);
    means = nan(1,3);
    for d = 1:min(3, size(loo_all,2))
        valid_d = isfinite(loo_all(:,d));
        if any(valid_d), means(d) = mean(loo_all(valid_d,d)); end
    end
    bar(means, 'FaceColor', [0.5 0.7 0.9]); set(gca,'XTickLabel',{'Y','S1','S3'},'YScale','log');
    ylabel(sprintf('Mean %s LOO error', tag)); title('Comparison of mean errors'); grid on;

    subplot(2,4,8);
    high_error_counts = zeros(1,3);
    thresh = 1e-2;
    for d = 1:min(3, size(loo_all,2))
        valid_d = isfinite(loo_all(:,d));
        if any(valid_d)
            high_error_counts(d) = sum(loo_all(valid_d,d) > thresh);
        end
    end
    bar(high_error_counts, 'FaceColor', [0.9 0.5 0.5]);
    set(gca,'XTickLabel',{'Y','S1','S3'});
    ylabel('Number of high-error nodes'); title(sprintf('High-error nodes (>%g)', thresh)); grid on;

    fname_comp = sprintf('%s_LOO_error_comprehensive_comparison_r%.2f', upper(tag), r);
    saveas(hfig3, [fname_comp '.png']); savefig(hfig3, [fname_comp '.fig']);
    fprintf('Combined comparison plot saved: %s\n', fname_comp);
end
