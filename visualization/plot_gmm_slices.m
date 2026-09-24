function plot_gmm_slices()
rng(100);
%% 1. Load saved results
fprintf('Loading saved analysis results...\n');

% Load coordinate data
coords = load('coordinate.mat');
coords = coords.coordinate;
xn = coords(1,:);  yn = coords(2,:);  zn = coords(3,:);
Nnodes = size(coords,2);

% Load analysis results
r = 20.00;
fprintf('Loading results for r = %.2f\n', r);

% Load IS results
is_file = sprintf('IS_results_r%.2f.mat', r);
if exist(is_file, 'file')
    load(is_file);
    fprintf('Loaded IS results: %s\n', is_file);
else
    error('IS results file not found: %s', is_file);
end
IS_data = IS_results;
% Use coordinates from the selected result, including newly generated runs.
coords = IS_data.coordinates;
% Use fixed RGB colors for components.
compCols = [ ...
    0.98 0.75 0.55;  % Soft light orange.
    0.98 0.78 0.85;  % Soft pink.
    0.65 0.90 0.70]; % Bright light green.
if exist('IS_data','var') && isfield(IS_data,'gmm') && ~isempty(IS_data.gmm)
    if IS_data.gmm.NumComponents > size(compCols,1)
        compCols = lines(IS_data.gmm.NumComponents);
    end
    % One-dimensional marginals: histogram, component curves, and mixture curve.
    d = size(IS_data.U_hot, 2);
    % Draw isosurfaces only.
    fprintf('\nGenerating isosurface plots for dimensions U1-U4...\n');
    
    for dim = 1:4  % U1, U2, U3, U4
        fprintf('Processing dimension U%d...\n', dim);
        
        % Generate isosurfaces at several quantiles of each fixed dimension.
        for fv = ["q25", "median", "q75"]
            fprintf('  - Generating isosurface for U%d=%s\n', dim, string(fv));
            % Generate only the isosurface plot.
            plot_gmm_3d_isosurface_vs_ellipsoids(IS_data, dim, fv, 0.90, compCols);
        end
    end

    fprintf('Isosurface plots for all dimensions are complete.\n');
else
    fprintf('Warning: no GMM data; skipping conditional 3D visualization.\n');
end

end

function create_gmm_conditional_slice(IS_data, fixed_dim, fixed_value, compCols)
% Plot 95% conditional ellipsoids in the remaining dimensions with matching colors.

    if nargin < 4 || isempty(compCols)
        compCols = lines(IS_data.gmm.NumComponents); % Fallback.
    end

    gmm = IS_data.gmm;
    U   = IS_data.U_hot;
    [~,d] = size(U);
    assert(d==4, 'This function assumes d=4. Extend selected_dims for other dimensions.');
    K    = gmm.NumComponents;
    baseW= gmm.ComponentProportion(:)';

    % Resolve the slice value.
    if ischar(fixed_value) || isstring(fixed_value)
        switch lower(string(fixed_value))
            case "median", a = median(U(:,fixed_dim));
            case "q25",    a = quantile(U(:,fixed_dim), 0.25);
            case "q75",    a = quantile(U(:,fixed_dim), 0.75);
            otherwise, error('Unrecognized fixed_value: %s', fixed_value);
        end
    else
        a = fixed_value;
    end

    % Select the remaining three dimensions.
    sel = setdiff(1:d, fixed_dim);   % Three dimensions.
    t = sqrt(chi2inv(0.95, 3));      % 95%：chi^2_3

    % Extract covariance matrices.
    function S = Sigma_k(k)
        if strcmpi(gmm.CovarianceType,'diagonal')
            if gmm.SharedCovariance, S = diag(gmm.Sigma(:));
            else,                    S = diag(gmm.Sigma(:,k));
            end
        else
            if gmm.SharedCovariance, S = gmm.Sigma(:,:,1);
            else,                    S = gmm.Sigma(:,:,k);
            end
        end
    end

    % Conditional means, covariances, and component weights.
    condMu  = zeros(3,K);
    condSig = zeros(3,3,K);
    like    = zeros(1,K);
    for k = 1:K
        mu = gmm.mu(k,:)';  S  = Sigma_k(k);
        mu_s = mu(sel);   mu_f = mu(fixed_dim);
        S_ss = S(sel, sel); S_sf = S(sel, fixed_dim);
        S_fs = S(fixed_dim, sel); S_ff = max(S(fixed_dim, fixed_dim), eps);

        condMu(:,k)    = mu_s + S_sf * (1/S_ff) * (a - mu_f);
        condSig(:,:,k) = S_ss - S_sf * (1/S_ff) * S_fs;

        like(k) = baseW(k) * normpdf(a, mu_f, sqrt(S_ff));
    end
    w = like / (sum(like) + eps);

    % Plot using the isosurface style.
    figure('Name', sprintf('Conditional GMM slice (dimension %d = %.3f)', fixed_dim, a), ...
           'Position', [50 50 1000 750]);
    ax = axes('Position',[0.08 0.08 0.88 0.88]); hold(ax,'on'); grid(ax,'on'); view(135,25);
    xlabel(sprintf('U_%d', sel(1)),'FontName','Times New Roman','FontSize',16); 
    ylabel(sprintf('U_%d', sel(2)),'FontName','Times New Roman','FontSize',16); 
    zlabel(sprintf('U_%d', sel(3)),'FontName','Times New Roman','FontSize',16);
    title(sprintf('Conditional slice: U_%d=%.3f (95%% ellipsoids in the remaining dimensions)', fixed_dim, a),'FontName','Times New Roman','FontSize',18);
    
    % Set the axis tick-label font size.
    set(ax,'FontName','Times New Roman','FontSize',14);

    [Xe,Ye,Ze] = sphere(36); E = [Xe(:)'; Ye(:)'; Ze(:)'];  % Moderate sphere resolution.
    fa = 0.03;  ea = 0.60;  % Transparent surfaces with visible edges.
    for k = 1:K
        Sk = condSig(:,:,k);
        [~,p] = chol(Sk);
        if p>0, Sk = Sk + 1e-9*eye(3); end
        [V,D] = eig(Sk);
        A = t * V * sqrt(max(D,0));
        E3 = A * E;
        Xp = reshape(E3(1,:) + condMu(1,k), size(Xe));
        Yp = reshape(E3(2,:) + condMu(2,k), size(Ye));
        Zp = reshape(E3(3,:) + condMu(3,k), size(Ze));

        % Improve 3D appearance through lighting and shading.
        surf(Xp, Yp, Zp, 'FaceColor', compCols(k,:), 'FaceAlpha', fa, ...
             'EdgeColor', compCols(k,:)*0.5, ...  % Lighter edges.
             'EdgeAlpha', ea, ...
             'FaceLighting', 'phong', ...        % Realistic lighting.
             'SpecularStrength', 0.4, ...         % Increase specular reflection.
             'DiffuseStrength', 0.8, ...          % Diffuse lighting strength.
             'AmbientStrength', 0.3);             % Ambient lighting strength.
        
        plot3(condMu(1,k), condMu(2,k), condMu(3,k), 'p', ...
              'Color', compCols(k,:)*0.8, ...
              'MarkerFaceColor', compCols(k,:), ...
              'MarkerSize', 12, ...
              'MarkerEdgeColor', compCols(k,:)*0.6, ...
              'LineWidth', 1.5);
    end

    % Use dummy legend handles to preserve colors and show conditional weights.
    hl = gobjects(1,K);
    for k=1:K
        hl(k) = plot3(nan,nan,nan,'p','MarkerFaceColor',compCols(k,:), ...
                      'MarkerEdgeColor',compCols(k,:)*0.6,'MarkerSize',12);
    end
    legtxt = arrayfun(@(k)sprintf('Component %d  (w|slice=%.1f%%)', k, 100*w(k)), 1:K, 'uni',0);
    legend(hl, legtxt, 'Location','northeastoutside','FontName','Times New Roman','FontSize',16);

    % Overlay samples near the slice for visual checking.
    delta = 0.15 * std(U(:,fixed_dim));
    idx = abs(U(:,fixed_dim) - a) <= delta;
    if any(idx)
        scatter3(U(idx, sel(1)), U(idx, sel(2)), U(idx, sel(3)), ...
                 8, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    end

    axis vis3d; material dull; 
    % Enhance lighting.
    camlight headlight; 
    camlight right; 
    camlight left;
    lighting phong;

    % Save a high-resolution figure.
    outname = sprintf('GMM_slice_dim%d_%s.png', fixed_dim, string(fixed_value));
    print(gcf, outname, '-dpng', '-r300');
    fprintf('Conditional ellipsoid plot saved: %s (300 DPI)\n', outname);
end

%% Resolve q25, median, or a numeric value to a slice coordinate.
function a = resolve_fixed_value(U, fixed_dim, fixed_value)
    if ischar(fixed_value) || isstring(fixed_value)
        switch lower(string(fixed_value))
            case "median", a = median(U(:,fixed_dim));
            case "q10",    a = quantile(U(:,fixed_dim), 0.10);
            case "q25",    a = quantile(U(:,fixed_dim), 0.25);
            case "q75",    a = quantile(U(:,fixed_dim), 0.75);
            case "q90",    a = quantile(U(:,fixed_dim), 0.90);
            otherwise, error('Unrecognized fixed_value: %s', fixed_value);
        end
    else
        a = fixed_value;
    end
end

%% Extract component covariance for full/diagonal and shared/unshared storage.
function S = Sigma_k_from_gmm(gmm, k)
    if strcmpi(gmm.CovarianceType,'diagonal')
        if gmm.SharedCovariance,  S = diag(gmm.Sigma(:));
        else,                     S = diag(gmm.Sigma(:,k));
        end
    else
        if gmm.SharedCovariance,  S = gmm.Sigma(:,:,1);
        else,                     S = gmm.Sigma(:,:,k);
        end
    end
end

%% Calculate the conditional GMM parameters for a fixed dimension and value.
function [muC, SigC, wC, sel, a] = compute_conditional_gmm(IS_data, fixed_dim, fixed_value)
    gmm = IS_data.gmm; U = IS_data.U_hot;
    a   = resolve_fixed_value(U, fixed_dim, fixed_value);
    d   = size(U,2);   sel = setdiff(1:d, fixed_dim);
    K = gmm.NumComponents; w0 = gmm.ComponentProportion(:)';

    muC  = zeros(3,K);
    SigC = zeros(3,3,K);
    like = zeros(1,K);

    for k=1:K
        mu = gmm.mu(k,:)'; S = Sigma_k_from_gmm(gmm,k);
        mu_s = mu(sel); mu_f = mu(fixed_dim);
        S_ss = S(sel, sel); S_sf = S(sel, fixed_dim);
        S_fs = S(fixed_dim, sel); S_ff = max(S(fixed_dim, fixed_dim), eps);
        muC(:,k)    = mu_s + S_sf*(1/S_ff)*(a - mu_f);
        SigC(:,:,k) = S_ss - S_sf*(1/S_ff)*S_fs;
        like(k)     = w0(k) * normpdf(a, mu_f, sqrt(S_ff));
    end
    wC = like / (sum(like) + eps);
end

%% Latent-space embedding (PCA/t-SNE) with 95% ellipse outlines.
function plot_latent_embedding(IS_data, compCols, method, pointSize)
    if nargin<3 || isempty(method), method='pca'; end
    if nargin<4 || isempty(pointSize), pointSize = 20; end

    U   = IS_data.U_hot;
    gmm = IS_data.gmm;
    K   = gmm.NumComponents;

    % Distinct fallback colors when compCols is missing or too short.
    if nargin<2 || isempty(compCols) || size(compCols,1)<K
        compCols = [ 0.20 0.47 0.75;   % Blue.
                     0.90 0.25 0.20;   % Red.
                     0.18 0.64 0.24;   % Green.
                     0.93 0.69 0.13;   % Yellow.
                     0.55 0.34 0.69;   % Purple.
                     0.40 0.40 0.40];  % Gray fallback.
        compCols = compCols(1:K,:);
    end

    % Reduce to two dimensions.
    switch lower(method)
        case 'tsne'
            try
                Y = tsne(U, 'NumDimensions',2, 'Perplexity',30, 'Standardize',true);
                A = []; muU = [];  % t-SNE has no linear projection matrix.
            catch
                warning('t-SNE is unavailable; falling back to PCA.'); method='pca';
            end
    end
    if ~exist('Y','var')   % Use PCA by default.
        muU = mean(U,1);
        [coeff, score] = pca(U, 'Centered',true);
        Y = score(:,1:2);           % Projected sample coordinates.
        A = coeff(:,1:2);           % Linear projection matrix (d-by-2).
    end

    % Color samples by the component with maximum posterior probability.
    P = posterior(gmm, U); [~, cid] = max(P,[],2);

    % Create the plot.
    f = figure('Name',sprintf('Latent embedding (%s)', upper(method)), ...
               'Position',[60 60 900 780], 'Color','w');
    ax = axes('Parent',f); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    % Larger filled points with transparency.
    for k=1:K
        scatter(ax, Y(cid==k,1), Y(cid==k,2), pointSize, ...
            'MarkerFaceColor', compCols(k,:), 'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha', 0.55);
    end

    % Project covariance ellipses only for the linear PCA embedding.
    if strcmpi(method,'pca') && ~isempty(A)
        t2 = chi2inv(0.95, 2);                  % 95% confidence ellipse.
        th = linspace(0, 2*pi, 200);
        unitCircle = [cos(th); sin(th)];        % 2×T

        for k=1:K
            % Mean and covariance in the original space.
            mu_k = gmm.mu(k,:);                 % 1×d
            Sk   = Sigma_k_from_gmm(gmm, k);    % d-by-d covariance matrix.

            % PCA projection: mu2d = (mu - mean(U))*A; S2d = A''*S*A.
            mu2d = (mu_k - muU) * A;            % 1×2
            S2d  = A' * Sk * A;                 % 2×2

            % Ellipse parameterization: mu + sqrt(t2)*V*sqrt(D)*circle.
            [V,D] = eig(S2d);  T = sqrt(t2) * V * sqrt(max(D,0));
            E2 = (T * unitCircle) + mu2d';      % 2×T

            plot(ax, E2(1,:), E2(2,:), '-', 'Color', compCols(k,:), 'LineWidth', 2.0);
        end
    end

    xlabel(ax, sprintf('%s-1', upper(method))); 
    ylabel(ax, sprintf('%s-2', upper(method)));
    title(ax, sprintf('Latent embedding (%s), colored by argmax posterior', upper(method)));

    % Preserve point and ellipse colors in the legend.
    hLeg = gobjects(1,K);
    for k=1:K
        hLeg(k) = plot(ax, nan, nan, '-', 'Color', compCols(k,:), 'LineWidth', 2.0);
    end
    legend(hLeg, arrayfun(@(k)sprintf('Comp %d',k),1:K,'uni',0), 'Location','northeastoutside');

    % Add some margin around the axes.
    axis(ax,'tight'); xlim(ax, expandlim(xlim(ax), 0.06)); ylim(ax, expandlim(ylim(ax), 0.06));

    % Save a high-resolution figure.
    outname = sprintf('GMM_embedding_%s_with_ellipses.png', lower(method));
    print(f, outname, '-dpng', '-r300');
    fprintf('Embedding plot with 95%% ellipses saved: %s (300 DPI)\n', outname);
end

% Helper for expanding the axis limits by a 6% margin.
function L = expandlim(L, frac)
    c = mean(L); r = range(L)/2; if r==0, r=1; end
    L = [c-(1+frac)*r, c+(1+frac)*r];
end

%% Conditional mean-shift arrows from from_value to to_value in the remaining 3D space.
function plot_conditional_shift_vectors(IS_data, fixed_dim, from_value, to_value, compCols)
    [muA, SigA, ~, sel, ~] = compute_conditional_gmm(IS_data, fixed_dim, from_value);
    [muB, ~, ~, ~,   ~] = compute_conditional_gmm(IS_data, fixed_dim, to_value);
    gmm = IS_data.gmm; K = gmm.NumComponents; if size(compCols,1)<K, compCols = lines(K); end

    t3 = sqrt(chi2inv(0.95,3)); [Xe,Ye,Ze]=sphere(48); E=[Xe(:)';Ye(:)';Ze(:)'];  % Increase resolution.
    f = figure('Name','Conditional shift vectors','Position',[50 50 1000 750],'Color','w');
    ax = axes('Parent',f,'Position',[0.08 0.08 0.88 0.88]); hold(ax,'on'); grid(ax,'on'); view(130,25); axis('vis3d');

    % Draw pale reference ellipsoids at the initial value.
    for k=1:K
        [V,D]=eig(SigA(:,:,k)); A=t3*V*sqrt(max(D,0)); E3=A*E;
        Xp=reshape(E3(1,:)+muA(1,k),size(Xe));
        Yp=reshape(E3(2,:)+muA(2,k),size(Ye));
        Zp=reshape(E3(3,:)+muA(3,k),size(Ze));
        surf(ax,Xp,Yp,Zp,'FaceColor',compCols(k,:),'FaceAlpha',0.12,...
            'EdgeColor',compCols(k,:)*0.4,'EdgeAlpha',0.40,...
            'FaceLighting','phong','SpecularStrength',0.3,...
            'DiffuseStrength',0.7,'AmbientStrength',0.4);
    end

    % Draw displacement arrows from muA to muB.
    for k=1:K
        u1 = muA(:,k); u2 = muB(:,k); du = u2 - u1;
        quiver3(ax, u1(1),u1(2),u1(3), du(1),du(2),du(3), 0, ...
            'Color', compCols(k,:), 'LineWidth', 2.0, 'MaxHeadSize', 0.5);
        plot3(ax, u1(1),u1(2),u1(3), 'o', 'Color',compCols(k,:), 'MarkerFaceColor',compCols(k,:), 'MarkerSize',6);
        plot3(ax, u2(1),u2(2),u2(3), 's', 'Color',compCols(k,:), 'MarkerFaceColor',compCols(k,:), 'MarkerSize',6);
    end

    xlabel(ax,sprintf('U_%d',sel(1)),'FontName','Times New Roman','FontSize',16); 
    ylabel(ax,sprintf('U_%d',sel(2)),'FontName','Times New Roman','FontSize',16); 
    zlabel(ax,sprintf('U_%d',sel(3)),'FontName','Times New Roman','FontSize',16);
    title(ax, sprintf('Shift of conditional means: U_%d: %s → %s', fixed_dim, string(from_value), string(to_value)),'FontName','Times New Roman','FontSize',18);
    legend(ax, arrayfun(@(k)sprintf('Comp %d',k),1:K,'uni',0), 'Location','northeastoutside','FontName','Times New Roman','FontSize',16);
    
    % Set the axis tick-label font size.
    set(ax,'FontName','Times New Roman','FontSize',14);
    out = sprintf('GMM_conditional_shift_U%d_%s_to_%s.png', fixed_dim, string(from_value), string(to_value));
    print(f, out, '-dpng', '-r300'); fprintf('Mean-shift arrow plot saved: %s (300 DPI)\n', out);
end

%% Isodensity shells supporting multiple probability masses, such as 90/95/97%.
% Same interface; target_mass can be a scalar or vector.

function plot_gmm_3d_isosurface_vs_ellipsoids(IS_data, fixed_dim, fixed_value, target_mass, compCols)
    if nargin<5 || isempty(compCols), compCols = lines(IS_data.gmm.NumComponents); end
    if isscalar(target_mass), target_mass = target_mass(:)'; end
    target_mass = sort(target_mass, 'ascend');  % Sort in ascending order.
    %isoFace = [0.75 0.55 0.90]; isoEdge = [0.45 0.20 0.65];
    %isoFace = [0.95 0.70 0.80]; isoEdge = [0.80 0.40 0.55];   % Light pink with dark pink edges.
    %isoFace = [0.55 0.75 0.95]; isoEdge = [0.25 0.45 0.75];   % Light blue with dark blue edges.
    isoFace = [0.85 0.92 0.98];   isoEdge = [0.50 0.70 0.92];   % Lighter blue.
    p.FaceAlpha = 0.10;
    p.EdgeAlpha = 0.30;
    gmm = IS_data.gmm; U = IS_data.U_hot;
    [muC,SigC,wC,sel,a] = compute_conditional_gmm(IS_data, fixed_dim, fixed_value);
    K = gmm.NumComponents;
    mu_mix = muC * wC(:);   % Column vector in the three selected coordinates.
    % Grid and mixture density.
    pad=1.0; m=mean(U(:,sel),1); s=std(U(:,sel),0,1);
    x1 = unique(double(linspace(m(1)-pad*3*s(1), m(1)+pad*3*s(1), 90)));
    x2 = unique(double(linspace(m(2)-pad*3*s(2), m(2)+pad*3*s(2), 90)));
    x3 = unique(double(linspace(m(3)-pad*3*s(3), m(3)+pad*3*s(3), 90)));
    [X,Y,Z] = meshgrid(x1,x2,x3); pts=[X(:) Y(:) Z(:)];
    PDF = zeros(size(pts,1),1);
    for k=1:K
        mu_k = muC(:,k)'; S_k = SigC(:,:,k);
        diff = pts - mu_k;
        [C,p]=chol(S_k,'lower'); if p>0, S_k=S_k+1e-9*eye(3); C=chol(S_k,'lower'); end
        y = diff / C';
        PDF = PDF + wC(k) * exp(-0.5*sum(y.^2,2) - sum(log(diag(C))) - 0.5*3*log(2*pi));
    end
    PDF3 = reshape(PDF, size(X));

    % Estimate density thresholds by sampling.
    nsamp = 6e4; num_k=max(1,round(nsamp*wC)); num_k(end)=nsamp-sum(num_k(1:end-1));
    pdf_samps=zeros(nsamp,1); cnt=0;
    for k=1:K
        nk=num_k(k); if nk<=0, continue; end
        Xk = mvnrnd(muC(:,k)', SigC(:,:,k), nk);
        pk=zeros(nk,1);
        for j=1:K, pk=pk+wC(j)*mvnpdf(Xk, muC(:,j)', SigC(:,:,j)); end
        pdf_samps(cnt+1:cnt+nk)=pk; cnt=cnt+nk;
    end
    pdf_samps = sort(pdf_samps(1:cnt));
    thr = arrayfun(@(q) quantile(pdf_samps, 1-q), target_mass);

    % Compact plot layout with reduced margins.
    f = figure('Name','Iso shells + conditional ellipsoids','Position',[50 50 1000 750],'Color','w');
    ax = axes('Parent',f,'Position',[0.16 0.18 0.70 0.72]); hold(ax,'on'); grid(ax,'on'); view(ax,130,25); axis(ax,'vis3d');

    % Draw the ellipsoids beneath the isosurfaces.
    t3 = sqrt(chi2inv(0.95,3)); [Xe,Ye,Ze]=sphere(48); E=[Xe(:)';Ye(:)';Ze(:)'];  % Increase resolution.
    fa=0.03; ea=0.60;  % Transparent surfaces with visible edges.
    for k=1:K
        [V,D]=eig(SigC(:,:,k)); A=t3*V*sqrt(max(D,0)); E3=A*E;
        Xp=reshape(E3(1,:)+muC(1,k),size(Xe));
        Yp=reshape(E3(2,:)+muC(2,k),size(Ye));
        Zp=reshape(E3(3,:)+muC(3,k),size(Ze));
        % Improve 3D appearance through lighting and shading.
        surf(ax,Xp,Yp,Zp,'FaceColor',compCols(k,:),'FaceAlpha',fa, ...
            'EdgeColor',compCols(k,:)*0.5,'EdgeAlpha',ea,...
            'FaceLighting','phong','SpecularStrength',0.4,...
            'DiffuseStrength',0.8,'AmbientStrength',0.3);
        % plot3(ax,muC(1,k),muC(2,k),muC(3,k),'p','Color',compCols(k,:)*0.8, ...
        %       'MarkerFaceColor',compCols(k,:), 'MarkerSize',12,...
        %       'MarkerEdgeColor',compCols(k,:)*0.6,'LineWidth',1.5);
    end

    % Draw the isosurfaces on top for emphasis.
    alphas = linspace(0.02, 0.08, numel(thr));  % Transparent surfaces emphasize the blue shell.
    ealph  = linspace(0.40, 0.65, numel(thr));   % Emphasize the edges.
    p90 = [];  % Store the handle of the 90% shell.
    for t = 1:numel(thr)
        p = patch(ax, isosurface(X,Y,Z,PDF3,thr(t)));
        isonormals(X,Y,Z,PDF3,p);
        set(p,'FaceColor',isoFace,'FaceAlpha',alphas(t), ...
              'EdgeColor',isoEdge,'EdgeAlpha',ealph(t),'LineStyle','-');
        uistack(p,'top');  % Move to the top layer.
        if abs(target_mass(t) - 0.90) < 1e-6   % Select the 90% shell.
            p90 = p;
        end
    end
    if ~isempty(p90)
        log_shell_centers_csv('shell_centers_90.csv', fixed_dim, a, target_mass, ...
                          sel, muC, wC, X, Y, Z, PDF3, thr, p90);
    end
    xlabel(ax,sprintf('U_%d',sel(1)),'FontName','Times New Roman','FontSize',26); 
    ylabel(ax,sprintf('U_%d',sel(2)),'FontName','Times New Roman','FontSize',26); 
    zlabel(ax,sprintf('U_%d',sel(3)),'FontName','Times New Roman','FontSize',26);
    
    % Set the axis tick-label font size.
    set(ax,'FontName','Times New Roman','FontSize',22);
    % Display only the analytical mixture center.
    h_center = plot3(ax, mu_mix(1), mu_mix(2), mu_mix(3), ...
    'o', 'MarkerSize', 9, 'LineWidth', 1.8, ...
    'MarkerFaceColor',[0.20 0.35 0.80], 'MarkerEdgeColor',[0.10 0.20 0.50]);
    % Include ellipsoids and isosurfaces in the legend.
    hLeg = gobjects(1,K+1);
    for k=1:K
        hLeg(k) = plot3(ax, nan,nan,nan,'p','MarkerFaceColor',compCols(k,:), ...
                      'MarkerEdgeColor',compCols(k,:)*0.6,'MarkerSize',12);
    end
    %hLeg(K+1) = plot3(ax, nan,nan,nan,'s','MarkerFaceColor',isoFace, 'MarkerEdgeColor',isoEdge);
    h_iso = plot3(ax, nan,nan,nan,'s','MarkerFaceColor',isoFace,'MarkerEdgeColor',isoEdge);
    % legend(ax, [h_iso, h_center], [arrayfun(@(k)sprintf('Component %d',k),1:K,'uni',0), {'Isosurface'}], ...
    %        'Location','northeastoutside','FontName','Times New Roman','FontSize',24);
    
    % title(ax, sprintf('Slice: U_%d = %.3f | shells: %s', fixed_dim, a, ...
    %       strjoin(compose('%.0f%%', 100*target_mass), ', ')),'FontName','Times New Roman','FontSize',26);
    outname = sprintf('GMM_isosurface_shells_dim%d_%s.png', fixed_dim, string(fixed_value));
    print(f, outname, '-dpng', '-r300'); fprintf('Multi-threshold isosurface plot saved: %s (300 DPI)\n', outname);
    % Calculate three centers of the 90% shell and export them to CSV.
    function log_shell_centers_csv(csvfile, fixed_dim, a, target_mass, sel, muC, wC, ...
                                   X, Y, Z, PDF3, thr, p_handle)
        % 1) Analytical center: conditional mixture mean.
        mu_mix = muC * wC(:);   % 3x1
    
        % 2) Geometric center: centroid of isosurface vertices.
        V = p_handle.Vertices;                  % N×3
        geo_cent = mean(V, 1).';                % 3x1
    
        % 3) Mass center: probability-weighted centroid of voxels above the threshold.
        mask = (PDF3 >= thr(1));
        PX = PDF3 .* mask;
        S  = sum(PX(:)) + eps;
        cx = sum(PX(:) .* X(:)) / S;
        cy = sum(PX(:) .* Y(:)) / S;
        cz = sum(PX(:) .* Z(:)) / S;
        mass_cent = [cx; cy; cz];               % 3x1
    
        % Append one row to the CSV.
        header = ['fixed_dim,fixed_value,target_mass,axis1,axis2,axis3,' ...
                  'mu_mix_1,mu_mix_2,mu_mix_3,geo_1,geo_2,geo_3,mass_1,mass_2,mass_3'];
        if ~isfile(csvfile)
            fid = fopen(csvfile, 'w'); fprintf(fid, '%s\n', header); fclose(fid);
        end
        fid = fopen(csvfile, 'a');
        fprintf(fid, '%d,%.10g,%.3f,U_%d,U_%d,U_%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            fixed_dim, a, target_mass(1), sel(1), sel(2), sel(3), ...
            mu_mix(1), mu_mix(2), mu_mix(3), ...
            geo_cent(1), geo_cent(2), geo_cent(3), ...
            mass_cent(1), mass_cent(2), mass_cent(3));
        fclose(fid);
    end

end

%% Conditional component weights as functions of the slice value.
% values: numeric coordinates or quantile labels such as q25, median, and q75.
function plot_slice_weight_curves(IS_data, fixed_dim, values, compCols)
    U = IS_data.U_hot; gmm = IS_data.gmm; K=gmm.NumComponents;
    if size(compCols,1) < K, compCols = lines(K); end

    % Convert values to numeric slice coordinates and create axis labels.
    if iscellstr(values) || (isstring(values) && ~isscalar(values))
        values = string(values);
    end
    if isnumeric(values)
        a_list = values(:)';
        xlab   = arrayfun(@(x)sprintf('%.3g',x), a_list, 'uni',0);
    else
        a_list = zeros(1, numel(values)); xlab = cell(1,numel(values));
        for i=1:numel(values)
            a_list(i) = resolve_fixed_value(U, fixed_dim, values(i));
            xlab{i}   = char(values(i));
        end
    end

    W = zeros(numel(a_list), K);
    for i=1:numel(a_list)
        [~, ~, wC] = compute_conditional_gmm(IS_data, fixed_dim, a_list(i));
        W(i,:) = wC;
    end

    f = figure('Name','Conditional weights vs slice','Position',[80 60 900 420],'Color','w');
    ax = axes('Parent',f); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for k=1:K
        plot(ax, 1:numel(a_list), W(:,k), '-o', 'Color', compCols(k,:), ...
            'MarkerFaceColor', compCols(k,:), 'LineWidth', 2);
    end
    xlim(ax, [0.8, numel(a_list)+0.2]); set(ax,'XTick',1:numel(a_list),'XTickLabel',xlab);
    ylabel(ax, 'w_k | slice'); xlabel(ax, sprintf('U_%d slice', fixed_dim));
    title(ax, sprintf('Conditional component weights across U_%d slices', fixed_dim));
    legend(ax, arrayfun(@(k)sprintf('Comp %d',k),1:K,'uni',0), 'Location','northeastoutside');
    out = sprintf('GMM_slice_weights_U%d.png', fixed_dim);
    print(f, out, '-dpng', '-r300'); fprintf('Conditional-weight plot saved: %s (300 DPI)\n', out);
end
