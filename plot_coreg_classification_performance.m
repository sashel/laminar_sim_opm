function plot_coreg_classification_performance(class_metric,space,axis,varargin)
% PLOT_COREG_CLASSIFICATION_PERFORMANCE Plot accuracy or bias across varying fiducial errors
% for whole brain and ROI analyses
% Adapted from the functions plot_snr_perc_corrrect and plot_snr_perc_pial
% by James Bonaiuto (https://github.com/jbonaiuto/laminar_sim)
%
% Use as plot_coreg_classification_performance('pial',35,1)
% where the first argument is the classification metric and the second the 
% number of measurement axes
%
%   plot_coreg_classification_performance(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

stem_dir = '/data/pt_np-helbling/layer_opm_sim/';
save_dir = '/data/pt_np-helbling/layer_opm_sim/results_figures_spmdev/';

invfoi = [10 30];
coreg_errs = {'','_coreg_1','_coreg_2','_coreg_3','_coreg_4'};

SNRs = [-20 -10];
methodnames = {'EBB','MSP'};

% parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', '/freesurfer_recons/');  % define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

% define critical t-value
dof = 199;
alpha = 1.0-0.05/2;
t_thresh = tinv(alpha, dof);

% original and downsampled white matter and pial surfaces
orig_white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','white.hc_PDw2.surf.gii');
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');

orig_pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','pial.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

simmeshes = {white_mesh,pial_mesh};
Nmesh = length(simmeshes);

% set colour scheme
addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
cm_all = brewermap(20, 'BuPu');
cm = colormap(cm_all(8:end-3,:));
for s = 1:length(SNRs)
    SNR = SNRs(s);
    for methind = 1:length(methodnames)
        method = methodnames{methind};
        sprintf('%s %d',method, SNR);

        % whole-brain analysistfrnakfurt 
        perc_nmb_unthresholded = zeros(1,length(coreg_errs)); % allocate space for mean classification performance in percent
        stderr_perc_nmb_unthresholded = zeros(1,length(coreg_errs)); % allocate space for standard deviation of classification performance in percent
        perc_nmb_significant = zeros(1,length(coreg_errs)); % allocate space for percentage of simulations yielding a significant classification

        disp('whole brain');
        nmb_unthresholded_per_coreg = zeros(1,Nmesh*params.nsims*length(coreg_errs));
        coreg_label = ordinal([ones(1,Nmesh*params.nsims)*1 ones(1,Nmesh*params.nsims)*2 ones(1,Nmesh*params.nsims)*3 ones(1,Nmesh*params.nsims)*4 ones(1,Nmesh*params.nsims)*5]);


        for i = 1:length(coreg_errs)
            % load whole-brain results
            data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_%s/',space,axis,coreg_errs{i}),...
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

            nmb_unthresholded_per_coreg((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);
            % classification performance significantly different from
            % chance level?
            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('If coreg error (mm) = %s, perc = %.4f, p = %.5f\n', coreg_errs{i}, perc_nmb_unthresholded(i), pout);
        end
        % plot results
        figure()
        hold on
        plot_fading_line(1:length(coreg_errs), perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '-');
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_coreg,coreg_label,'model','ordinal');
            disp('logistic regression free energy');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
            disp('logistic regression free energy');
        end

        % ROI-based analysis
        disp('ROI');
        perc_nmb_unthresholded = zeros(1,length(coreg_errs));
        stderr_perc_nmb_unthresholded = zeros(1,length(coreg_errs));
        perc_nmb_significant = zeros(1,length(coreg_errs));
        nmb_unthresholded_per_coreg = zeros(1,Nmesh*params.nsims*length(coreg_errs));

        for i = 1:length(coreg_errs)
            SNR = SNRs(s);
            % load ROI-based results
            data_dir = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d%s/',space,axis,coreg_errs{i}),...
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
            nmb_unthresholded_per_coreg((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;

            nmb_significant = [abs(class_t(1:params.nsims))>t_thresh; abs(class_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);

            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('If coreg error (mm) = %s, perc = %.4f, p = %.5f\n', coreg_errs{i},perc_nmb_unthresholded(i), pout);

        end
        plot_fading_line(1:length(coreg_errs), perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '--');
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_coreg,coreg_label,'model','ordinal');
            disp('logistic regression roi');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
            disp('logistic regression roi');
        end

        hold off
            xlabel('Co-registration error (mm)','FontSize',14)
        set(gca,'Xtick',1:5,'XLim',[.5 5.5])
        set(gca,'XtickLabels',{'0','1','2','3','4'})
        if strcmp(class_metric,'correct')
            ylabel('% Correct','FontSize',14)
            set(gca,'Ytick',[0, 0.2, 0.4, 0.6, 0.8,1]*100,'YLim',[35 105])
        elseif strcmp(class_metric,'pial')
            ylabel('% Classified Pial','FontSize',14)
            set(gca,'Ytick',[0.3, 0.5, 0.7, 0.9]*100,'YLim',[25 105])
        end
        set(gcf,'Colormap',cm)
        c = colorbar;
        c.Title.String = '% signif';
        c.Title.FontSize = 18;
        c.Ticks = [0, 0.2, 0.4, 0.6, 0.8,1];
        c.TickLabels = {'0','20','40','60','80','100'};
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',21)
        set(gca,'linewidth',2.4)
        title([method ' (SNR ' num2str(SNR) ' dB)'],'FontSize',21,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
        set(gcf,'color','w')
        if strcmp(class_metric,'correct')
            saveas(gcf,sprintf('%scoreg_SNR_%d_space_%d_axis_%d_ds_%s_perc_correct.svg', save_dir, SNR, space, axis,method),'svg')
        elseif strcmp(class_metric,'pial')
            saveas(gcf,sprintf('%scoreg_SNR_%d_space_%d_axis_%d_ds_%s_perc_pial.svg', save_dir, SNR, space, axis,method),'svg')
        end
    end
end
