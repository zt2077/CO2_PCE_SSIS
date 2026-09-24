function plot_comparison()
% Compare the saved nodal SSIS and MCS failure probabilities.
load('IS_results_r20.00.mat', 'IS_results');
load('MCS_results_r20.00.mat', 'MCS_results');
load('comparison_data_r20.00.mat', 'comparison_data');
Pf_IS = IS_results.Pf;
Pf_MCS = MCS_results.Pf;
Diff = comparison_data.Diff;
r = comparison_data.r;
assert(numel(Pf_IS) == numel(Pf_MCS), 'Result dimensions do not match.');
create_accuracy_analysis(Pf_IS, Pf_MCS, Diff, r);

end

function rel_errors  = create_accuracy_analysis(Pf_IS, Pf_MCS, Diff, r)
    % Calculate accuracy metrics
    mae = mean(abs(Diff));
    rmse = sqrt(mean(Diff.^2));
    mse = mean(Diff.^2);
    % Relative errors are undefined at nodes with zero MCS probability.
    positive_ref = isfinite(Pf_MCS) & Pf_MCS > 0 & isfinite(Diff);
    mape = mean(abs(Diff(positive_ref)./Pf_MCS(positive_ref))) * 100;
    fprintf('Relative-error metrics use %d/%d nodes with positive MCS probability.\n', ...
        nnz(positive_ref), numel(Pf_MCS));
    
    correlation = corrcoef(Pf_MCS, Pf_IS);
    corr_coeff = correlation(1,2);
    
    ss_res = sum(Diff.^2);
    ss_tot = sum((Pf_MCS - mean(Pf_MCS)).^2);
    r_squared = 1 - (ss_res / ss_tot);
    
    bias = mean(Diff);
    precision = std(Diff);
    
    % Create accuracy analysis figure
    figure('Position', [100, 100, 1200, 800], 'Name', 'SS+IS versus MCS accuracy assessment');
    
    % Subplot 1: Scatter plot with regression
    subplot(2, 3, 1);
    scatter(Pf_MCS, Pf_IS, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    plot([min(Pf_MCS), max(Pf_MCS)], [min(Pf_MCS), max(Pf_MCS)], 'r--', 'LineWidth', 2);
    p = polyfit(Pf_MCS, Pf_IS, 1);
    y_fit = polyval(p, Pf_MCS);
    plot(Pf_MCS, y_fit, 'g-', 'LineWidth', 2);
    xlabel('MCS failure probability', 'FontSize', 11);
    ylabel('SS+IS failure probability', 'FontSize', 11);
    title('Scatter plot and regression line', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Data', '1:1 line', 'Regression line', 'Location', 'best', 'FontSize', 10);
    grid on;
    
    % Add R² and regression equation
    text(0.1*max(Pf_MCS), 0.9*max(Pf_IS), sprintf('R² = %.4f\ny = %.2fx + %.2e', r_squared, p(1), p(2)), ...
         'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Subplot 2: Residual analysis
    subplot(2, 3, 2);
    residuals = Diff;
    scatter(Pf_MCS, residuals, 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    plot([min(Pf_MCS), max(Pf_MCS)], [0, 0], 'k--', 'LineWidth', 1.5);
    xlabel('MCS failure probability', 'FontSize', 11);
    ylabel('Residual (SS+IS - MCS)', 'FontSize', 11);
    title('Residual analysis', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Subplot 3: Residual histogram
    subplot(2, 3, 3);
    histogram(residuals, 30, 'FaceColor', [0.7 0.9 1.0], 'EdgeColor', 'black');
    xlabel('Residual', 'FontSize', 11);
    ylabel('Count', 'FontSize', 11);
    title('Residual distribution', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Subplot 4: Q-Q plot
    subplot(2, 3, 4);
    qqplot(residuals);
    xlabel('Standard normal quantiles');
    ylabel('Residual quantiles');
    title('Residual Q-Q plot (normality assessment)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Subplot 5: Relative error distribution
    subplot(2, 3, 5);
    rel_errors = (residuals(positive_ref) ./ Pf_MCS(positive_ref)) * 100;
    histogram(rel_errors, 30, 'FaceColor', [0.7 1.0 0.7], 'EdgeColor', 'black');
    xlabel('Relative error (%)', 'FontSize', 11);
    ylabel('Count', 'FontSize', 11);
    title('Relative-error distribution', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Subplot 6: Accuracy metrics summary
    subplot(2, 3, 6);
    metrics = [mae, rmse, mape, r_squared*100, abs(bias), precision];
    metric_names = {'MAE', 'RMSE', 'MAPE(%)', 'R²(%)', '|Bias|', 'Precision'};
    bar(metrics, 'FaceColor', [1.0 0.6 0.0], 'EdgeColor', 'black');
    set(gca, 'XTickLabel', metric_names);
    ylabel('Value', 'FontSize', 11);
    title('Accuracy metrics', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    xtickangle(45);
    
    sgtitle(sprintf('SS+IS versus MCS accuracy assessment (r = %.2f)', r), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save accuracy analysis figure
    accuracy_filename = sprintf('SSIS_vs_MCS_accuracy_r%.2f.png', r);
    saveas(gcf, accuracy_filename, 'png');
    fprintf('Accuracy assessment saved to: %s\n', accuracy_filename);
    
    % Print accuracy summary
    fprintf('\n=== Accuracy metrics ===\n');
    fprintf('Mean absolute error (MAE): %.2e\n', mae);
    fprintf('Root mean squared error (RMSE): %.2e\n', rmse);
    fprintf('Mean absolute percentage error (MAPE): %.2f%%\n', mape);
    fprintf('Correlation coefficient (R): %.4f\n', corr_coeff);
    fprintf('Coefficient of determination (R squared): %.4f\n', r_squared);
    fprintf('Systematic bias: %.2e\n', bias);
    fprintf('Precision: %.2e\n', precision);
end

