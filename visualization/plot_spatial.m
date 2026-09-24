function plot_spatial()
% Draw comparable spatial maps using coordinates stored with each result.
a = load('MCS_results_r20.00.mat','MCS_results');
b = load('IS_results_r20.00.mat','IS_results');
results = {a.MCS_results, b.IS_results};
names = {'MCS','SSIS'};
all_pf = [a.MCS_results.Pf(:); b.IS_results.Pf(:)];
limits = [min(all_pf), max(all_pf)];
if limits(1) == limits(2), limits(2) = limits(1) + eps; end
for k = 1:2
    res = results{k};
    xyz = res.coordinates;
    assert(size(xyz,2) == numel(res.Pf), 'Coordinate and probability sizes differ.');
    f = figure('Color','w','Position',[100,100,1000,750]);
    scatter3(xyz(1,:),xyz(2,:),xyz(3,:),22,res.Pf(:),'filled');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    cb = colorbar; cb.Label.String = 'Failure probability';
    clim(limits); axis equal; view(3); grid on;
    title(sprintf('%s, injection rate = %.0f kg/s', names{k}, res.r));
    exportgraphics(f, sprintf('Spatial_Pf_%s.png',names{k}), 'Resolution',300);
end
end
