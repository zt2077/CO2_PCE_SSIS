function toy_example()
% Stage A shows successive SS levels approaching the limit state.
% A fitted GMM improves the visual comparison between Stages B and C.
%% Benchmark-1: SS → GMM → IS three-stage evolution (fixed)
rng(42);

%% ----- Problem definition -----
gfun = @(x1,x2) 5 - x2 - 0.5*(x1 - 0.1).^2;   % limit-state
logp  = @(X) -0.5*sum(X.^2,2) - log(2*pi);    % log pdf of N(0,I)
N0 = 2000; p0 = 0.25; rho = 0.9;               % SS params

%% ----- Subset Simulation -----
X = randn(N0,2); 
g = gfun(X(:,1),X(:,2));
b = quantile(g,p0); 
idx = g<=b; 
seeds = X(idx,:);                             % Initial-level seeds.
levels = struct('X',X,'g',g,'b',b);
fprintf('L0 threshold b=%.3f (%d/%d)\n',b,sum(idx),N0);

% Sample each chain separately instead of passing all seeds to a single chain.
for L = 2:4
    Nchains   = size(seeds,1);
    NperChain = floor(N0 / Nchains);
    Xc = rw_mh_all(seeds, b, rho, NperChain); % Multi-chain sampler.
    g  = gfun(Xc(:,1), Xc(:,2));
    bq = quantile(g, p0);
    % Stop adaptively when the threshold reaches zero.
    if bq <= 0
        b = 0;
        levels(L) = struct('X', Xc, 'g', g, 'b', b);
        seeds     = Xc(g <= 0, :);           % Failure seeds at the final level.
        fprintf('L%d threshold b=%.3f -> clipped to 0, stop.\n', L-1, bq);
        break
    else
        b = bq;
        levels(L) = struct('X', Xc, 'g', g, 'b', b);
        seeds     = Xc(g <= b, :);
        fprintf('L%d threshold b=%.3f (%d chains × %d)\n', L-1, b, Nchains, NperChain);
    end
    % b  = quantile(g, p0);
    % levels(L) = struct('X', Xc, 'g', g, 'b', b);
    % seeds     = Xc(g<=b, :);
    % fprintf('L%d threshold b=%.3f (%d chains × %d)\n', L-1, b, Nchains, NperChain);
end
X_last = levels(end).X; 
g_last = levels(end).g;

%% ----- Fit GMM (K=2) to last-level samples -----
K=2;
opts=statset('MaxIter',2000,'TolFun',1e-7);
gmm=fitgmdist(X_last,K,'RegularizationValue',1e-3,...
    'Options',opts,'Start','plus');
fprintf('GMM fitted: %d comps.\n',K);

%% ----- Importance Sampling -----
Nis=1e5;
Xq=random(gmm,Nis);
I_F = gfun(Xq(:,1),Xq(:,2))<=0;
logw = log(I_F + 1e-300) + logp(Xq) - log(pdf(gmm,Xq));
w    = exp(logw - max(logw));                  % Numerical stabilization.
% Restore the log-weight scale for the ordinary IS probability estimate.
% Keep w scaled for the unchanged visualization and resampling below.
Pf_IS = exp(max(logw)) * mean(w);
fprintf('Pf(IS)=%.3e\n',Pf_IS);

%% ----- Visualization (three-stage evolution) -----
[x1,x2]=meshgrid(linspace(-6,6,250));
p  = exp(-0.5*(x1.^2+x2.^2))/(2*pi);
g0 = gfun(x1,x2);

figure('Color','w','Position',[100 100 1200 430]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
% Use Times New Roman consistently.
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');
set(groot, 'DefaultLegendFontName', 'Times New Roman');
set(groot, 'DefaultColorbarFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 12);
set(groot, 'DefaultTextFontSize', 12);
% Use TeX for subscripts while retaining the selected font.
set(groot, 'DefaultTextInterpreter', 'tex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
set(groot, 'DefaultLegendInterpreter', 'tex');

% Stage A: SS last level
nexttile;
axis equal; xlim([-6 6]); ylim([-6 6]); hold on; %grid on;
xl = xlim; yl = ylim;
% True-color grayscale background independent of other colormaps.
% Assume x1, x2, and nonnegative p have already been evaluated.
cmap = grayGrad(256);                        % Custom grayscale gradient defined below.
pN = (p - min(p(:))) / (max(p(:)) - min(p(:)) + eps);  % [0,1]
idx = max(1, min(256, round(pN*255) + 1));   % Map values to indices 1 through 256.
bgRGB = ind2rgb(idx, cmap);                  % True-color RGB.
hbg = image(x1(1,:), x2(:,1), bgRGB);        % Alternatively use imagesc.
set(gca,'YDir','normal');
set(hbg, 'AlphaData', 0.95);                 % Semi-transparent.
uistack(hbg, 'bottom');                      % Send to the bottom layer.

nAvail = numel(levels);
Luse = min(max(nAvail-1,0), 2);              % Use at most level L2.
%if Luse >= 1
if Luse >= 0
    %b_list = arrayfun(@(k) levels(k).b, 2:(1+Luse));  % b^(1), b^(2), ...
    b_list = arrayfun(@(k) levels(k).b, 1:(1+Luse));
    % Light blue/orange fills with blue/orange boundaries.
    % fillCols = [0.82 0.90 1.00;   1.00 0.92 0.80];
    % edgeCols = [0.00 0.45 0.74;   0.85 0.33 0.10];
    % fillCols = [0.90 0.95 1.00;   0.82 0.90 1.00;   1.00 0.92 0.80];
    % edgeCols = [0.50 0.50 0.50;   0.00 0.45 0.74;   0.85 0.33 0.10];
    fillCols = [0.90 1.00 0.90;   0.82 0.90 1.00;   1.00 0.92 0.80];
    edgeCols = [0.20 0.60 0.20;   0.00 0.45 0.74;   0.85 0.33 0.10];
    %for ell = 1:Luse
    for ell = 1:(Luse+1)
        bL = b_list(ell);
        % Shade F^(ell) = {g<=b^(ell)}, the region above the parabola.
        b_plot = max(bL,0);
        %shade_subset_level(gca, xlim, ylim, bL, fillCols(ell,:), 0.22);
        shade_subset_level(gca, xl, yl, b_plot, fillCols(ell,:), 0.22);
        % Solid boundary g(x)=b^(ell).
        %xs = linspace(xlim(1), xlim(2), 1000);
        xs = linspace(xl(1), xl(2), 1000);
        %yb = 5 - bL - 0.5*(xs - 0.1).^2;
        yb = 5 - b_plot - 0.5*(xs - 0.1).^2;
        plot(xs, yb,'Color', edgeCols(ell,:), 'LineWidth', 2.0);
        % Label b^(ell) at the vertex.
        [~, iApex] = max(yb);
        text(xs(iApex), yb(iApex)+0.25, sprintf('{\\it b}\\rm^{(%d)}', ell-1), ...
            'Color', edgeCols(ell,:), 'FontSize', 13, ...
            'HorizontalAlignment','center','Interpreter','tex');
        % text(xs(iApex), yb(iApex)+0.25, sprintf('b^{(%d)}', ell), ...
        %     'Color', edgeCols(ell,:), 'FontSize', 13, ...
        %     'HorizontalAlignment','center','Interpreter','tex');
    end
end

%shade_failure_region(gca, xlim, ylim);    % Draw the gray failure region first.
contour(x1,x2,p,10,'k:'); hold on;
contour(x1,x2,g0,[0 0],'k--','LineWidth',1.5);
scatter(X_last(:,1),X_last(:,2),3,g_last,'filled');
ax = gca; ax.Box = 'on'; ax.Layer = 'top';
ax.LineWidth = 0.75; ax.TickDir = 'In'; ax.FontSize = 14;
%title('Stage A: SS last-level'); axis equal; xlim([-6 6]); ylim([-6 6]);
xlabel('\it x\rm_1', 'Interpreter','tex','FontSize',16);
ylabel('\it x\rm_2', 'Interpreter','tex','FontSize',16);

% Stage B: GMM proposal
nexttile;
axis equal; xlim([-6 6]); ylim([-6 6]); hold on; %grid on;
% True-color grayscale background independent of other colormaps.
% Assume x1, x2, and nonnegative p have already been evaluated.
cmap = grayGrad(256);                        % Custom grayscale gradient defined below.
pN = (p - min(p(:))) / (max(p(:)) - min(p(:)) + eps);  % [0,1]
idx = max(1, min(256, round(pN*255) + 1));   % Map values to indices 1 through 256.
bgRGB = ind2rgb(idx, cmap);                  % True-color RGB.
hbg = image(x1(1,:), x2(:,1), bgRGB);        % Alternatively use imagesc.
set(gca,'YDir','normal');
set(hbg, 'AlphaData', 0.95);                 % Semi-transparent.
uistack(hbg, 'bottom');                      % Send to the bottom layer.

shade_failure_region(gca, xlim, ylim);    % Gray failure region.
contour(x1,x2,p,10,'k:'); hold on;
contour(x1,x2,g0,[0 0],'k--','LineWidth',1.5);
Zgmm = reshape(pdf(gmm,[x1(:) x2(:)]), size(x1));
levB = linspace(min(Zgmm(:)), max(Zgmm(:)), 12);   % Define contour levels for Stage B.
levB = levB(2:end);     % Remove the last contour level.
contour(x1,x2,Zgmm,levB,'LineWidth',1.0);
ax = gca; ax.Box = 'on'; ax.Layer = 'top';
ax.LineWidth = 0.75; ax.TickDir = 'In'; ax.FontSize = 14;
%title('Stage B: GMM proposal (K=2)'); axis equal; xlim([-6 6]); ylim([-6 6]);
xlabel('\it x\rm_1', 'Interpreter','tex','FontSize',16);
ylabel('\it x\rm_2', 'Interpreter','tex','FontSize',16);
%draw_gmm_major_axes(gca, gmm, 2.8, 'Color',[0 0 0], 'LineStyle','-', 'LineWidth',1.6);

%% ===== Stage C：IS reweighted density (weighted-resampling GMM version) =====
nexttile;
axis equal; xlim([-6 6]); ylim([-6 6]); hold on;

% True-color grayscale background independent of other colormaps.
cmap = grayGrad(256);
pN = (p - min(p(:))) / (max(p(:)) - min(p(:)) + eps);
idx = max(1, min(256, round(pN*255) + 1));
bgRGB = ind2rgb(idx, cmap);
hbg = image(x1(1,:), x2(:,1), bgRGB);
set(gca,'YDir','normal');
set(hbg, 'AlphaData', 0.95);
uistack(hbg, 'bottom');

% Gray-blue shading of the failure region.
xl = xlim; yl = ylim;
shade_failure_region(gca, xl, yl);

% Extract failure samples and weights.
Xf = Xq(I_F, :);
wf = w(I_F);
if isempty(Xf)
    warning('No failure samples were obtained; increase Nis or inflate the GMM covariance and retry.');
    return
end
wf = wf / sum(wf);
% Generate approximate samples by weighted resampling.
wf = wf / sum(wf);
Nrep = round(wf * 10000);   % Approximately 1e4 samples; adjust as needed.
Nrep(Nrep==0) = 1;          % Retain at least one copy of each failure sample.
Xrep = repelem(Xf, Nrep, 1);
% Refit a GMM to the weighted resample.
opts = statset('MaxIter',2000,'TolFun',1e-7);
gmm_IS = fitgmdist(Xrep, 2, ...
    'RegularizationValue', 1e-3, ...
    'Options', opts, ...
    'Start', 'plus');
% Draw the Stage C contours.
Zc = reshape(pdf(gmm_IS, [x1(:) x2(:)]), size(x1));
% Reuse the Stage B levels if available; otherwise generate them.
if exist('levB','var')
    levC = levB;
else
    levC = linspace(min(Zc(:)), max(Zc(:)), 12);
    levC = levC(2:end);  % Remove the outermost contour.
end
contour(x1, x2, Zc, levC, 'LineWidth', 1.2);
% Overlay the background and limit-state boundary.
contour(x1, x2, p, 10, 'k:'); hold on;
contour(x1, x2, g0, [0 0], 'k--', 'LineWidth', 1.5);

% Axis formatting.
axis equal; xlim([-6 6]); ylim([-6 6]);
ax = gca;
ax.Box = 'on';
ax.Layer = 'top';
ax.LineWidth = 0.75;
ax.TickDir = 'In';
ax.FontSize = 14;
xlabel('\it x\rm_1', 'Interpreter','tex','FontSize',16);
ylabel('\it x\rm_2', 'Interpreter','tex','FontSize',16);
%title('Stage C: IS-weighted (resampled GMM)', 'FontSize',13);
%draw_gmm_major_axes(gca, gmm_IS, 2.8, 'Color',[0 0 0], 'LineStyle','-', 'LineWidth',1.6);
%print(gcf, '-dpng', '-r300', 'three_stage_evolution.png');
print(gcf, '-dtiff', '-r300', 'three_stage_evolution.tiff');

%% ===== Helper functions =====
end

function Xall = rw_mh_all(seeds, b, rho, Nper)
    % Concatenate samples from multiple chains.
    Nchains = size(seeds,1);
    d = size(seeds,2);
    Xall = zeros(Nchains*Nper, d);
    k = 1;
    for c = 1:Nchains
        Xc = ss_chain(seeds(c,:), b, rho, Nper);    % Pass a single 1-by-2 seed to each chain.
        Xall(k:k+Nper-1, :) = Xc;
        k = k + Nper;
    end
end

function Xc = ss_chain(x0, b, rho, Nper)
    % SS random-walk Metropolis-Hastings chain constrained to {g<=b}.
    step = sqrt(1 - rho^2);
    d = numel(x0); 
    Xc = zeros(Nper, d);
    x = x0(:)';                                    % Ensure a 1-by-2 row vector.
    for i = 1:Nper
        xprop = rho*x + step*randn(1,d);           % 1×2
        if 5 - xprop(2) - 0.5*(xprop(1)-0.1)^2 <= b   % g(xprop) <= b
            x = xprop;
        end
        Xc(i,:) = x;                               % Assign the 1-by-2 sample.
    end
end

function shade_failure_region(ax, xlimv, ylimv)
% Shade g(x)<=0 above the parabola on the specified axes.
    x = linspace(xlimv(1), xlimv(2), 800);
    y = 5 - 0.5*(x - 0.1).^2;          % g(x)=0: x2 = 5 - 0.5*(x1-0.1)^2
    % Clip to the y limits and construct a polygon above the parabola.
    y = max(y, ylimv(1));              % Respect the lower limit.
    y = min(y, ylimv(2));              % Respect the upper limit.
    xpoly = [x, fliplr(x)];
    ypoly = [y,  ylimv(2)*ones(size(y))];
    % h = patch('Parent', ax, ...
    %           'XData', xpoly, 'YData', ypoly, ...
    %           'FaceColor', [0.75 0.75 0.75], ...   % Light gray.
    %           'EdgeColor', 'none', ...
    %           'FaceAlpha', 0.25);                  % Adjustable opacity, typically 0.15-0.35.
     h = patch('Parent', ax, ...
              'XData', xpoly, 'YData', ypoly, ...
              'FaceColor', [0.70 0.85 1.00], ...  % Light blue.
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.25);                 % Adjustable opacity, typically 0.15-0.35.
    uistack(h, 'bottom');                          % Place beneath the contours.
end

function cmap = blueGrad(n)
% Smooth gradient from light blue outside to dark blue at the center.
    if nargin < 1, n = 256; end
    stops = [0;   0.25; 0.50; 0.75; 1.00];
    cols  = [0.90 0.95 1.00;   % Very light blue at the outside.
             0.80 0.90 1.00;
             0.65 0.80 1.00;
             0.40 0.65 0.98;
             0.20 0.45 0.95];  % Deep blue at the center.
    xi    = linspace(0,1,n)';
    cmap  = interp1(stops, cols, xi, 'pchip');
end

function cmap = grayGrad(n)
% Smooth gradient from light gray outside to dark gray at the center.
    if nargin < 1, n = 256; end
    stops = [0;   0.25; 0.50; 0.75; 1.00];
    cols  = [0.92 0.92 0.92;   % Very light gray at the outside.
             0.85 0.85 0.85;
             0.72 0.72 0.72;
             0.45 0.45 0.45;
             0.18 0.18 0.18];  % Deep gray at the center.
    xi   = linspace(0,1,n)';
    cmap = interp1(stops, cols, xi, 'pchip');
end

% Shade each intermediate subset {g(x) <= b}.
function shade_subset_level(ax, xlimv, ylimv, b, faceColor, alphaVal)
    % g(x)=b: x2 = 5 - b - 0.5*(x1-0.1)^2
    x = linspace(xlimv(1), xlimv(2), 1200);
    yb = 5 - b - 0.5*(x - 0.1).^2;
    % Clip to the current y limits.
    yb = max(yb, ylimv(1)); 
    yb = min(yb, ylimv(2));
    % Polygon above the parabola up to the top boundary.
    xpoly = [x, fliplr(x)];
    ypoly = [yb, ylimv(2)*ones(size(yb))];
    h = patch('Parent', ax, 'XData', xpoly, 'YData', ypoly, ...
              'FaceColor', faceColor, 'EdgeColor', 'none', ...
              'FaceAlpha', alphaVal);
    uistack(h, 'bottom');
end

function draw_gmm_major_axes(ax, gm, scale, varargin)
% Draw the major-axis direction of each GMM component.
% Support cell, 2D, and 3D covariance storage in older MATLAB versions.

    hold(ax, 'on');
    K = gm.NumComponents;
    MU = gm.mu;

    for k = 1:K
        % Handle the three covariance storage formats.
        if iscell(gm.Sigma)
            S = gm.Sigma{k};
        elseif ndims(gm.Sigma) == 3
            S = gm.Sigma(:,:,k);
        else
            S = gm.Sigma;  % Single-component case.
        end

        S = (S + S.') / 2;               % Enforce numerical symmetry.
        [V, D] = eig(S);
        [lamMax, idx] = max(diag(D));    % Direction associated with the largest eigenvalue.
        v = V(:, idx) / norm(V(:, idx)); % Normalize the direction vector.

        a = scale * sqrt(max(lamMax, 0));  % Half-length scale.
        c = MU(k, :)';                     % Component center.
        p1 = c - a * v; 
        p2 = c + a * v;

        plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], varargin{:});
    end
end
