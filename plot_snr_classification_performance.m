function plot_snr_classification_performance(class_metric,space,axis,varargin)
% PLOT_SNR_CLASSIFICATION_PERFORMANCE  Plot accuracy over SNR levels for whole brain and
% ROI analysis
%
% Use as plot_snr_classification_performance ('pial',55,1)
% where the first argument is the classification metric, the second is the
% inter-sensor distance in mm, and the third the number of measurement axes
%
%   plot_SNR_classification_performance(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

stem_dir = '/data/pt_np-helbling/layer_opm_sim/';
save_dir = '/data/pt_np-helbling/layer_opm_sim/results_figures_spmdev/';

invfoi = [10 30];
SNRs = [ -50 -40 -30 -20 -10 -5];

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', '/data/pt_user-helbling_ticket017439/helbling/NormativeMEG/Data/Freesurfer6.0.0_Recons/');  % define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

% critical t value
dof = 199;
alpha = 1.0-0.05/2;
t_thresh = tinv(alpha, dof);

methodnames = {'EBB','MSP'};
orig_white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','white.hc_PDw2.surf.gii');
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');

orig_pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','pial.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

simmeshes = {white_mesh,pial_mesh};
Nmesh = length(simmeshes);

addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
cm_all = brewermap(14, 'BuPu');
cm = colormap(cm_all(3:end-2,:));

disp([num2str(space)]);
for methind = 1:length(methodnames)
    method = methodnames{methind};
    disp(method);

    figure();
    hold on;

    perc_nmb_unthresholded = zeros(1,length(SNRs));
    stderr_perc_nmb_unthresholded = zeros(1,length(SNRs));
    perc_nmb_significant = zeros(1,length(SNRs));

    disp('whole brain');
    nmb_unthresholded_per_snr = zeros(1,Nmesh*params.nsims*length(SNRs));
    snr_label = ordinal([ones(1,Nmesh*params.nsims)*1 ones(1,Nmesh*params.nsims)*2 ones(1,Nmesh*params.nsims)*3 ones(1,Nmesh*params.nsims)*4 ones(1,Nmesh*params.nsims)*5 ones(1,Nmesh*params.nsims)*6]);

    for s = 1:length(SNRs)
        SNR = SNRs(s);
        data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_ds_spmdev/',space,axis),...
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
        nmb_unthresholded_per_snr((s-1)*Nmesh*params.nsims+1:s*Nmesh*params.nsims) = nmb_unthresholded;
        perc_nmb_significant(s) = mean(nmb_significant);
        perc_nmb_unthresholded(s) = mean(nmb_unthresholded);
        stderr_perc_nmb_unthresholded(s) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));

        pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
        fprintf('SNR = %.2f dB, perc = %.4f, p = %.5f\n', SNR, perc_nmb_unthresholded(s), pout);

    end

    try
        [~,~,stats] = mnrfit(nmb_unthresholded_per_snr,snr_label,'model','ordinal');

        disp('logistic regression free energy');
        fprintf('beta = %.8f\n',stats.beta(end))
        fprintf('p = %.8f\n',stats.p(end))
    catch
    end
    plot_fading_line(SNRs, perc_nmb_unthresholded.*100, ...
        stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '-');

    disp('ROI');
    perc_nmb_unthresholded = zeros(1,length(SNRs));
    stderr_perc_nmb_unthresholded = zeros(1,length(SNRs));
    perc_nmb_significant = zeros(1,length(SNRs));
    nmb_unthresholded_per_snr = zeros(1,Nmesh*params.nsims*length(SNRs));

    for s = 1:length(SNRs)
        SNR = SNRs(s);
        data_dir = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_ds_spmdev/',space,axis),...
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
        nmb_unthresholded_per_snr((s-1)*Nmesh*params.nsims+1:s*Nmesh*params.nsims) = nmb_unthresholded;

        nmb_significant = [abs(class_t(1:params.nsims))>t_thresh; abs(class_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
        perc_nmb_unthresholded(s) = mean(nmb_unthresholded);
        stderr_perc_nmb_unthresholded(s) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));

        perc_nmb_significant(s) = mean(nmb_significant);

        pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
        fprintf('SNR = %.2f dB, perc = %.4f, p = %.5f\n', SNR, perc_nmb_unthresholded(s), pout);

    end
    plot_fading_line(SNRs, perc_nmb_unthresholded.*100, ...
        stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '--');

    try
        [~,~,stats] = mnrfit(nmb_unthresholded_per_snr,snr_label,'model','ordinal');
        disp('logistic regression roi');
        fprintf('beta = %.8f\n',stats.beta(end))
        fprintf('p = %.8f\n',stats.p(end))
    catch
        disp('logistic regression roi');
    end

    hold off;
    title(['Spatial sampling ' num2str(space) 'mm'],'FontSize',16,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')

    xlabel('SNR (dB)','FontSize',14)
    if strcmp(class_metric,'correct')
        ylabel('% Correct','FontSize',14)
    elseif strcmp(class_metric,'pial')
        ylabel('% Classified Pial','FontSize',14)
    end

    set(gca,'Xtick',SNRs,'XLim',[-55 0])
    set(gca,'Ytick',[0, 0.2, 0.4, 0.6, 0.8,1]*100,'YLim',[35 105])
    set(gcf,'Colormap',cm)
    c = colorbar;
    c.Title.String = '% Significant';
    c.Title.FontSize = 10;
    c.Ticks = [0, 0.2, 0.4, 0.6, 0.8,1];
    c.TickLabels = {'0','20','40','60','80','100'};

    set(gcf,'color','w')

    if strcmp(class_metric,'correct')
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_correct.fig', save_dir, space, axis, method))
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_correct.jpg', save_dir, space, axis, method),'jpg')
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_correct.svg', save_dir, space, axis, method),'svg')
    elseif strcmp(class_metric,'pial')
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_pial.fig', save_dir, space, axis, method))
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_pial.jpg', save_dir, space, axis, method),'jpg')
        saveas(gcf,sprintf('%ssnr_space_%d_axis_%d_ds_%s_perc_pial.svg', save_dir, space, axis, method),'svg')
    end
end
