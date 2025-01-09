function plot_offset_classification_performance(class_metric,space,axis,varargin)
% PLOT_OFFSET_CLASSIFICATION_PERFORMANCE Plot accuracy or bias across varying offset levels
% for whole brain and ROI analyses
% Adapted from the functions plot_snr_perc_corrrect and plot_snr_perc_pial
% by James Bonaiuto (https://github.com/jbonaiuto/laminar_sim)
%
% Use as plot_offset_classification_performance('pial',55,1)
% where the first argument is the classification metric, the second ,the inter-sensor distance in mm
% and the third the number of measurement axes
%
%   plot_offset_classification_performance(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

stem_dir = '/data/pt_np-helbling/layer_opm_sim/';
save_dir = '/data/pt_np-helbling/layer_opm_sim/results_figures_spmdev/';

invfoi = [10 30];
offsets = {'','_offset_20','_offset_30','_offset_40'};

SNRs = [-20 -10 -5];
methodnames = {'EBB','MSP'};

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', '/freesurfer_recons/');  % define default values
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

% original and downsampled white matter surface
orig_white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','white.hc_PDw2.surf.gii');
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');

orig_pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','pial.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

simmeshes = {white_mesh,pial_mesh};
Nmesh = length(simmeshes);

addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
cm_all = brewermap(14, 'BuPu');
cm = colormap(cm_all(3:end-2,:));

for s = 1:length(SNRs)
    SNR = SNRs(s);
    for methind = 1:length(methodnames)
        method = methodnames{methind};
        sprintf('%s %d',method, SNR);

        perc_nmb_unthresholded = zeros(1,length(offsets));
        stderr_perc_nmb_unthresholded = zeros(1,length(offsets));
        perc_nmb_significant = zeros(1,length(offsets));

        disp('whole brain');
        nmb_unthresholded_per_offset = zeros(1,Nmesh*params.nsims*length(offsets));
        offset_label = ordinal([ones(1,Nmesh*params.nsims)*1 ones(1,Nmesh*params.nsims)*2 ones(1,Nmesh*params.nsims)*3 ones(1,Nmesh*params.nsims)*4]);

        for i = 1:length(offsets)
            data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_%s/',space,axis,offsets{i}),...
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

            nmb_unthresholded_per_offset((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);

            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('If offset (mm) = %s, perc = %.4f, p = %.5f\n', offsets{i}, perc_nmb_unthresholded(i), pout);

        end
        figure()
        hold on
        plot_fading_line(1:length(offsets), perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '-');
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_offset,offset_label,'model','ordinal');
            disp('logistic regression free energy');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
            disp('logistic regression free energy');
        end

        disp('ROI');
        perc_nmb_unthresholded = zeros(1,length(offsets));
        stderr_perc_nmb_unthresholded = zeros(1,length(offsets));
        perc_nmb_significant = zeros(1,length(offsets));
        nmb_unthresholded_per_offset = zeros(1,Nmesh*params.nsims*length(offsets));

        for i = 1:length(offsets)
            SNR = SNRs(s);
            data_dir = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_%s/',space,axis,offsets{i}),...
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
            nmb_unthresholded_per_offset((i-1)*Nmesh*params.nsims+1:i*Nmesh*params.nsims) = nmb_unthresholded;

            nmb_significant = [abs(class_t(1:params.nsims))>t_thresh; abs(class_t(params.nsims+1:Nmesh*params.nsims))>t_thresh];
            perc_nmb_unthresholded(i) = mean(nmb_unthresholded);
            stderr_perc_nmb_unthresholded(i) = std(nmb_unthresholded)/sqrt(length(nmb_unthresholded));
            perc_nmb_significant(i) = mean(nmb_significant);

            pout = myBinomTest(sum(nmb_unthresholded),length(nmb_unthresholded),0.5,'two');
            fprintf('If offset (mm) = %s, perc = %.4f, p = %.5f\n', offsets{i}, perc_nmb_unthresholded(i), pout);

        end
        plot_fading_line(1:length(offsets), perc_nmb_unthresholded.*100, ...
            stderr_perc_nmb_unthresholded.*100, perc_nmb_significant, cm, '--');
        try
            [~,~,stats] = mnrfit(nmb_unthresholded_per_offset,offset_label,'model','ordinal');
            disp('logistic regression roi');
            fprintf('beta = %.3f\n',stats.beta(end))
            fprintf('p = %.3f\n',stats.p(end))
        catch
            disp('logistic regression roi');
        end

        hold off
        title([method ' (SNR ' num2str(SNR) ' dB)'],'FontSize',16,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
        xlabel('Scalp-sensor offset (mm)','FontSize',14)
        set(gca,'Xtick',1:4,'XLim',[.5 4.5])
        set(gca,'XtickLabels',{'6.5','20','30','40'})
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
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_correct.fig', save_dir, SNR, space, axis,method))
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_correct.svg', save_dir, SNR, space, axis,method),'svg')
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_correct.jpg', save_dir, SNR, space, axis,method),'jpg')
        elseif strcmp(class_metric,'pial')
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_pial.fig', save_dir, SNR, space, axis,method))
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_pial.svg', save_dir, SNR, space, axis,method),'svg')
            saveas(gcf,sprintf('%soffsets_SNR_%d_space_%d_axis_%d_ds_%s_perc_pial.jpg', save_dir, SNR, space, axis,method),'jpg')
        end
    end
end
