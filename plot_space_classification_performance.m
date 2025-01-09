function plot_space_classification_performance(class_metric, axis,varargin)
% PLOT_SPACE_CLASSIFICATION_PERFORMANCE Plot accuracy or bias across varying levels of spatial sampling
% for whole brain and ROI analyses
%
% Use as plot_space_classification_performance('pial',1)
% where the first argument is the classification metric and the second the
% number of measurement axes
%
%   plot_space_classification_performance(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

stem_dir = '/data/pt_np-helbling/layer_opm_sim/';
save_dir = '/data/pt_np-helbling/layer_opm_sim/results_figures_spmdev/';

space = 25:10:55; % inter-sensor distances in mm
invfoi = [10 30]; % frequency band 
SNRs = [-40 -30 -20 -10 -5]; % trial-wise signal-to-noise-ratios in dB

methodnames = {'EBB','MSP'}; % empirical Bayesian beamformer and Multiple Sparse Priors source reconstruction approaches

% parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', '/data/pt_user-helbling_ticket017439/helbling/NormativeMEG/Data/Freesurfer6.0.0_Recons/');  % define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

% define critical t-value
dof = 199; % degrees of freedom for 20 trials
alpha = 1.0-0.05/2; % significance level of alpha = 5%, two-tailed test
t_thresh = tinv(alpha, dof);

% original and downsampled white matter surface
orig_white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','white.hc_PDw2.surf.gii');
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');

orig_pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','pial.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

simmeshes = {white_mesh,pial_mesh};
Nmesh = length(simmeshes);

% set colour scheme
addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
cm_all = brewermap(14, 'BuPu');
cm = colormap(cm_all(3:end-2,:));

for s = 1:length(SNRs) % loop across SNRs
    SNR = SNRs(s);
    for methind = 1:length(methodnames) % loop across spurce reconstruction approaches
        method = methodnames{methind};
        disp([method ' ' num2str(SNR)]);

        figure();
        hold on;

        perc_nmb_unthresholded = zeros(1,length(space)); % allocate space for mean classification performance in percent (of either classification correct or biased towards the pial surface)
        stderr_perc_nmb_unthresholded = zeros(1,length(space)); % allocate space for the standard deviation of classification performance in percent
        perc_nmb_significant = zeros(1,length(space)); % allocate space for percentage of simulations yielding a significant classifications
        
        disp('whole brain');
        nmb_unthresholded_per_space = zeros(1,Nmesh*params.nsims*length(space));
        space_label = ordinal([ones(1,Nmesh*params.nsims)*1 ones(1,Nmesh*params.nsims)*2 ones(1,Nmesh*params.nsims)*3 ones(1,Nmesh*params.nsims)*4]);

        for i = 1:length(space)
            % load whole-brain results
            data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d/',space(i),axis),...
                sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
                invfoi(1),invfoi(2),SNR,params.dipole_moment));
            load(data_file,'allcrossF');

            nmb_unthresholded = zeros(1,Nmesh*params.nsims);
            nmb_significant = zeros(1,Nmesh*params.nsims);
            for simmeshind = 1:Nmesh
                if strcmp(class_metric,'correct')
                    class_performance = squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind)-allcrossF(simmeshind,1:params.nsims,(2-simmeshind)+1,methind));
                elseif strcmp(class_metric,'pial')
                    class_performance = squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
                end
                nmb_unthresholded((simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = class_performance>0;
                nmb_significant((simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = abs(class_performance)>3;
            end
            nmb_unthresholded_per_space((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;
            nmb_unthresholded_per_space_f = nmb_unthresholded_per_space;
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);
            % classification performance significantly different than
            % chance level?
            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('Spatial sampling (mm) = %.2f, %s = %.4f, p = %.5f\n', space(i), class_metric, perc_nmb_unthresholded(i), pout);
        end
        plot_fading_line(space, perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '-');
        % logistic regression across inter-sensor distances
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_space,space_label,'model','ordinal');
            disp('logistic regression free energy');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
        end
        

        disp('ROI');
        perc_nmb_unthresholded = zeros(1,length(space));
        stderr_perc_nmb_unthresholded = zeros(1,length(space));
        perc_nmb_significant = zeros(1,length(space));
        nmb_unthresholded_per_space = zeros(1,Nmesh*params.nsims*length(space));

        for i = 1:length(space)
            SNR = SNRs(s);
            data_dir = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d/',space(i),axis),...
                sprintf('f%d_%d_SNR%d_dipolemoment%d',invfoi(1),invfoi(2),SNR,params.dipole_moment));
            if strcmp(class_metric,'correct')
                class_t = get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
                    white_mesh, orig_pial_mesh, orig_white_mesh, 'recompute_trials', false);
                nmb_unthresholded = [class_t(1:params.nsims)<0; class_t(params.nsims+1:Nmesh*params.nsims)>0];

            elseif strcmp(class_metric,'pial')
                class_t = get_wmpial_t(data_dir, method, params.nsims, pial_mesh, ...
                    white_mesh, orig_pial_mesh, orig_white_mesh, 'recompute_trials', false);
                nmb_unthresholded = class_t>0;
            end
            nmb_unthresholded_per_space((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;
            nmb_unthresholded_per_space_roi = nmb_unthresholded_per_space;

            nmb_significant = [abs(class_t(1:params.nsims))>t_thresh; abs(class_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);

            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('Spatial sampling (mm) = %.2f, perc = %.4f, p = %.5f\n', space(i), perc_nmb_unthresholded(i), pout);
        end
        % plot results
        plot_fading_line(space, perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '--');
        % logistic regression across inter-sensor distances
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_space,space_label,'model','ordinal');
            disp('logistic regression roi');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
        end

        hold off;
        title([method ' (SNR ' num2str(SNR) ' dB)'],'FontSize',16,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
        xlabel('Spatial sampling (mm)','FontSize',14)
        set(gca,'Xtick',space,'XLim',[20 60])
        if strcmp(class_metric,'correct')
            ylabel('% Correct','FontSize',14)
            set(gca,'Ytick',[0, 0.2, 0.4, 0.6, 0.8,1]*100,'YLim',[35 105])
        elseif strcmp(class_metric,'pial')
            ylabel('% Classified Pial','FontSize',14)
            set(gca,'Ytick',[0.3, 0.5, 0.7, 0.9]*100,'YLim',[25 105])
        end
        set(gcf,'Colormap',cm)
        c = colorbar;
        c.Title.String = '% Significant';
        c.Title.FontSize = 10;
        c.Ticks = [0, 0.2, 0.4, 0.6, 0.8,1];
        c.TickLabels = {'0','20','40','60','80','100'};
        set(gcf,'color','w')

        if strcmp(class_metric,'correct')
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_correct.fig', save_dir, SNR, axis,method))
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_correct.svg', save_dir, SNR, axis,method),'svg')
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_correct.jpg', save_dir, SNR, axis,method),'jpg')
        elseif strcmp(class_metric,'pial')
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_pial.fig', save_dir, SNR, axis,method))
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_pial.svg', save_dir, SNR, axis,method),'svg')
            saveas(gcf,sprintf('%sspace_SNR_%d_axis_%d_ds_%s_perc_pial.jpg', save_dir, SNR, axis,method),'jpg')

        end
        % Exact McNemar tests to evaluate differences between whole-brain and
        % ROI-based analyses
        for ii = 1:length(space)
            fprintf('Spatial sampling (mm) = %.2f', space(ii));
            both_wrong = sum(~nmb_unthresholded_per_space_f((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims) & ~nmb_unthresholded_per_space_roi((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims));
            f_wrong_not_roi = sum(~nmb_unthresholded_per_space_f((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims) & nmb_unthresholded_per_space_roi((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims));
            roi_wrong_not_f = sum(nmb_unthresholded_per_space_f((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims) & ~nmb_unthresholded_per_space_roi((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims));
            both_right = sum(nmb_unthresholded_per_space_f((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims) & ~nmb_unthresholded_per_space_roi((ii-1)*Nmesh*params.nsims+1:ii*Nmesh*params.nsims));
            McNemarextest([both_wrong,f_wrong_not_roi,roi_wrong_not_f,both_right],2,0.05);
        end
    end
end 
