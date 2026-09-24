%% The model function
function [Yall,dataS1_all,dataS3_all] = support_unit(X,Nnodes,r_injection,output_dir)
    import com.comsol.model.*             % bring in the core COMSOL model classes  
    import com.comsol.model.util.*        % bring in ModelUtil, etc. 
    [Nsim, k] = size(X);
    %Yfail = zeros(Nsim,1);
    kk = Nnodes;
    %Y = zeros(Nsim,kk);
    %selected_nodes = [10, 100, 200];   % Customize the selected nodes.
    selected_nodes = [10];
    pp=length(selected_nodes);
    Yall = zeros(Nsim,kk);
    dataS1_all = zeros(Nsim,kk);
    dataS3_all = zeros(Nsim,kk);
    tau_max_all = zeros(Nsim,kk);
    tau_MC_all = zeros(Nsim,kk);
    Y = zeros(Nsim,pp);
    global model
    %Y = zeros(kk,Nsim);
    dir_str = [char(output_dir), filesep];
    if ~isfolder(output_dir), mkdir(output_dir); end
    for a = 1:Nsim
        a
        % 1) Extract KLE weights
        xi = X(a,:).';               % [k×1]
        %%  Assign to material and solve
        %model.param.set('r_injection', '15[kg/s]', 'Injection rate');
        model.param.set('r_injection', sprintf('%g[kg/s]', r_injection), 'Injection rate');
        model.component('comp1').physics('solid').feature('lemm1').set('E', xi(1));
        model.component('comp1').physics('solid').feature('lemm1').set('nu', xi(2));
        model.component('comp1').physics('solid').feature('lemm1').feature('soil1').set('cohesion', xi(3));
        model.component('comp1').physics('solid').feature('lemm1').feature('soil1').set('internalphi', xi(4));
        ModelUtil.showProgress(true);
        model.sol('sol5').runAll;
        model.sol('sol6').runAll;
        exportName = 'data1';

        % Remove an existing export object if present.
        try
            model.result.export.remove(exportName);
        catch
            % Nothing to remove if the object does not exist.
        end
        % Create a new export object.
        model.result.export.create(exportName, 'Data');
        model.result.export(exportName).set('data', 'dset6');
        model.result.export(exportName).set('expr', {'w', 'solid.sp1', 'solid.sp3'});
        model.result.export(exportName).set('unit', {'m', 'N/m^2', 'N/m^2'});
        % Set the output file path.
        % Ensure the data_all directory exists.
        data_dir = fullfile(dir_str, 'data_all');
        if ~exist(data_dir, 'dir')
            mkdir(data_dir);
        end
        filename = sprintf('data_all_%d.txt', a);
        % filename = 'data_all.txt';
        full_filename = fullfile(dir_str, 'data_all', filename);
        model.result.export(exportName).set('filename', full_filename);
        % Execute the export.
        model.result.export(exportName).run();
        %data_all = mphinterp(model, 'w', 'coord', coords_lagrange);
        %data_all = mpheval(model,'w','dataset','dset5');
        data_all = readmatrix(full_filename, 'CommentStyle', '%');
        assert(size(data_all,2) >= 12 && size(data_all,1) == Nnodes, ...
            'Unexpected COMSOL export layout; verify dataset dset6 and node ordering.');
        Yall(a,:) = data_all(:,10)';
        Y(a,:) = Yall(a, selected_nodes);
        % dataS1_all(a,:) = mpheval(model,'solid.sp1').dl;
        % dataS3_all(a,:) = mpheval(model,'solid.sp3').dl;
        % dataS1_all(a,:) = mpheval(model,'solid.sp1').d1;
        dataS1_all(a,:) =data_all(:,11)';
        size(dataS1_all)
        dataS3_all(a,:) =data_all(:,12)';
        %tau_max_all(a,:) = (dataS1_all(a,:) - dataS3_all(a,:))/2;
        %sigma_n_all(a,:) = (dataS1_all(a,:) + dataS3_all(a,:))/2;
        %tau_max(a,:) = tau_max_all(a, selected_nodes);
        %sigma_n(a,:) = sigma_n_all(a,selected_nodes);
        c = xi(3);
        phi = xi(4);
        %tanphi = tan(phi/2+45/180*pi);
        %sigma_n_eff = min(sigma_n_all(a,:), 0);  % compression only
        %tau_MC_all(a,:) = c - sigma_n_eff .* tan(phi);
        %tau_MC(a,:) = tau_MC_all(a, selected_nodes);
    end
    %Pf = sum(Yfail)/Nsim;
    %mphlaunch(model)
    save([dir_str,'Yall_generated.mat'], 'Yall');
    save([dir_str,'S1_generated.mat'], 'dataS1_all');
    save([dir_str,'S3_generated.mat'], 'dataS3_all');
    writematrix(Yall, [dir_str,'Yall_generated.csv']);
    writematrix(dataS1_all, [dir_str,'S1_generated.csv']);
    writematrix(dataS3_all, [dir_str,'S3_generated.csv']);
    filename = 'data_all_1.txt';
    data_all = readmatrix(fullfile(dir_str, 'data_all', filename), 'CommentStyle', '%');
    coordinate = zeros(3, size(data_all, 1));
    coordinate(1,:) = data_all(:,1)';
    coordinate(2,:) = data_all(:,2)';
    coordinate(3,:) = data_all(:,3)';
    % Save as a MAT file.
    save([dir_str,'coordinate.mat'], 'coordinate');
end
