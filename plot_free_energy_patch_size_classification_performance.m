function plot_free_energy_patch_size_classification_performance(varargin)
% PLOT_FREE_ENERGY_PATCH_SIZE_CLASSIFICATION_PERFORMANCE  Plot bias and accuracy for 
% whole-brain patch size simulations
%
% Use as plot_free_energy_patch_size_classification_performance()
%
%   plot_free_energy_patch_size_classification_performance(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

stem_dir = '/data/pt_np-helbling/layer_opm_sim/';
save_dir = '/data/pt_np-helbling/layer_opm_sim/results_figures_spmdev/';

space = 35;
axis = 1;
invfoi = [10 30];

SNRs = [-50,-40,-30,-20,-10,-5];

% parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', <SURF_DIR>,'methind',1);  % define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames = {'EBB','MSP'};

% mesh file names
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

allmeshes = {white_mesh, pial_mesh};
Nmesh = length(allmeshes);

% combinations of simulation and reconstruction patch sizes
sim_patch_sizes = [5 5 10 10];
recon_patch_sizes = [5 10 5 10];

% allocate space for percentage correct and biased variables
f_correct_thresholded = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);
f_correct_unthresholded = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);
f_correct_significant = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);
f_pial_unthresholded = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);
f_pial_thresholded = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);
f_pial_significant = zeros(length(methodnames),length(sim_patch_sizes),...
    length(SNRs),params.nsims*Nmesh);

for methind = 1:length(methodnames) % loop across source reconstruction approaches
    method = methodnames{methind};
    for s = 1:length(SNRs) % loop across SNR levels
        SNR = SNRs(s);
        for idx = 1:length(sim_patch_sizes)
            sim_patch_size = sim_patch_sizes(idx);
            recon_patch_size = recon_patch_sizes(idx);

            % load whole-brain results for current patch size combination
            if sim_patch_size == 5&& recon_patch_size == 5
                data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_ds_spmdev/',space,axis),...
                    sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
                    invfoi(1),invfoi(2),SNR,params.dipole_moment));
            elseif sim_patch_size == 5&& recon_patch_size == 10
                data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_recon_patchsize_10_ds_spmdev/',space,axis),...
                    sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
                    invfoi(1),invfoi(2),SNR,params.dipole_moment));
            elseif sim_patch_size == 10 && recon_patch_size == 5
                data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_sim_patchsize_10_ds_spmdev/',space,axis),...
                    sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
                    invfoi(1),invfoi(2),SNR,params.dipole_moment));
            elseif sim_patch_size == 10 && recon_patch_size == 10
                data_file = fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_recon_sim_patchsize_10_ds_spmdev/',space,axis),...
                    sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
                    invfoi(1),invfoi(2),SNR,params.dipole_moment));
            end
            load(data_file,'allcrossF');


            for simmeshind = 1:Nmesh
                % get classification biases for the whole-brain analysis
                pialF = squeeze(allcrossF(simmeshind,1:params.nsims,2,methind));
                whiteF = squeeze(allcrossF(simmeshind,1:params.nsims,1,methind));
                pialWhiteF = pialF-whiteF;
                
                f_pial_thresholded(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = pialWhiteF>3;
                f_pial_unthresholded(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = pialWhiteF>0;
                f_pial_significant(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = abs(pialWhiteF)>3;

                % get number of correct classifications 
                trueF = squeeze(allcrossF(simmeshind,1:params.nsims,simmeshind,methind));
                otherF = squeeze(allcrossF(simmeshind,1:params.nsims,2-simmeshind+1,methind));
                trueOtherF = trueF-otherF;
                f_correct_thresholded(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = trueOtherF>3;
                f_correct_unthresholded(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = trueOtherF>0;
                f_correct_significant(methind,idx,s,...
                    (simmeshind-1)*params.nsims+1:simmeshind*params.nsims) = abs(trueOtherF)>3;
            end
        end
    end


    % plot percentage correct across SNRs for all patch size combination 
    styles = {'-','-','--','--'};
    colors = {'b','r','r','b'};
    addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
    
    figure
    hold all;
    for i = 1:length(sim_patch_sizes)
        perc_correct_unthresholded_stderr = squeeze(std(f_correct_unthresholded(methind,i,:,:),[],4))./sqrt(Nmesh*params.nsims).*100;
        perc_correct_unthresholded = squeeze(mean(f_correct_unthresholded(methind,i,:,:),4)).*100.0;
        perc_correct_significant = squeeze(mean(f_correct_significant(methind,i,:,:),4));       

        % print percentage correct and its statistical significance
        for s = 1:length(SNRs)
            SNR = SNRs(s);
            pout = myBinomTest(sum(squeeze(f_correct_unthresholded(methind,i,s,:))),length(f_correct_unthresholded(methind,i,s,:)),0.5,'two');
            fprintf('SNR %d, Sim = %d, Recon = %d, correct = %.4f, p = %.5f\n', SNR,sim_patch_sizes(i),recon_patch_sizes(i), perc_correct_unthresholded(s)/100, pout);
        end

        % set colour
        switch colors{i}
            case 'b'
                cm_all = brewermap(20, 'BuPu');
                cm = colormap(cm_all(8:end-3,:));

            case 'r'
                cm_all = brewermap(20, 'OrRd');
                cm = colormap(cm_all(8:end-3,:));
        end
        % do the actual plotting
        plot_fading_line(SNRs, perc_correct_unthresholded, ...
            perc_correct_unthresholded_stderr, perc_correct_significant, cm, styles{i});
    end

    % annotate the figure
    xlabel('SNR (dB)','FontSize',14)
    ylabel('% Classified Correct','FontSize',14)
    set(gca,'Xtick',SNRs,'XLim',[-55 0])
    ylim([-5 105]); 
    set(gca,'Ytick',[0, 0.2, 0.4, 0.6, 0.8,1]*100,'YLim',[-5 105])
    set(gcf,'Colormap',cm)
    c = colorbar;
    c.Title.String = '% signif';
    c.Title.FontSize = 18;
    c.Ticks = [0, 0.2, 0.4, 0.6, 0.8,1];
    c.TickLabels = {'0','20','40','60','80','100'};
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',21)
    set(gca,'linewidth',2.4)
    set(gcf,'color','w')
    title(method,'FontSize',24,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
    
    % save figure on classification accuracy across congruent and incongruent patch size combinations
    saveas(gcf,sprintf('%spatch_sizes_space_%d_axis_%d_ds_%s_perc_correct_OrRd_colorbar.jpg', save_dir, space, axis,method),'jpg')
    saveas(gcf,sprintf('%spatch_sizes_space_%d_axis_%d_ds_%s_perc_correct_OrRd_colorbar.svg', save_dir, space, axis,method),'svg')

    % plot percentage bias to pial surface across SNRs for all patch size combinations 
    figure
    hold all;
    for i = 1:length(sim_patch_sizes)
        perc_pial_unthresholded_stderr = squeeze(std(f_pial_unthresholded(methind,i,:,:),[],4))./sqrt(Nmesh*params.nsims).*100;
        perc_pial_unthresholded = squeeze(mean(f_pial_unthresholded(methind,i,:,:),4)).*100.0;
        perc_pial_significant = squeeze(mean(f_pial_significant(methind,i,:,:),4));      

        % print percentage biased to pial surface and its statistical significance
        for s = 1:length(SNRs)
            SNR = SNRs(s);
            pout = myBinomTest(sum(squeeze(f_pial_unthresholded(methind,i,s,:))),length(f_pial_unthresholded(methind,i,s,:)),0.5,'two');
            fprintf('SNR %d, Sim = %d, Recon = %d, pial = %.4f, p = %.5f\n', SNR, sim_patch_sizes(i),recon_patch_sizes(i), perc_pial_unthresholded(s)/100, pout);
        end

        % set colour
        switch colors{i}
            case 'b'
                cm_all = brewermap(20, 'BuPu');
                cm = colormap(cm_all(8:end-3,:));

            case 'r'
                cm_all = brewermap(20, 'OrRd');
                cm = colormap(cm_all(8:end-3,:));
        end
        % do the actual plotting
        plot_fading_line(SNRs, perc_pial_unthresholded, ...
            perc_pial_unthresholded_stderr, perc_pial_significant, cm, styles{i});
    end

    % annotate figure
    xlabel('SNR (dB)','FontSize',14)
    ylabel('% Classified Pial','FontSize',14)
    ylim([0 105]);
    set(gca,'Xtick',SNRs,'XLim',[-55 0])
    set(gca,'Ytick',[0, 0.2, 0.4, 0.6, 0.8,1]*100,'YLim',[-5 105])
    set(gcf,'Colormap',cm)
    c = colorbar;
        c.Title.String = '% signif';
        c.Title.FontSize = 18;
        c.Ticks = [0, 0.2, 0.4, 0.6, 0.8,1];
        c.TickLabels = {'0','20','40','60','80','100'};
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',21)
        set(gca,'linewidth',2.4)
               set(gcf,'color','w')
    title(method,'FontSize',24,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
         
    % save figure on classification bias across congruent and incongruent patch size combinations
    saveas(gcf,sprintf('%spatch_sizes_space_%d_axis_%d_ds_%s_perc_pial.jpg', save_dir, space, axis, method),'jpg')
    saveas(gcf,sprintf('%spatch_sizes_space_%d_axis_%d_ds_%s_perc_pial.svg', save_dir, space, axis, method),'svg')


% do stats based on exact McNemar tests
    disp('Comparing classifiers');
    for s = 1:length(SNRs)
        fprintf('SNR = %df, Sim patch size 5mm, Reconstruct 5mm-Reconstruct 10mm\n',SNRs(s));
        both_wrong = sum(~squeeze(f_correct_thresholded(methind,1,s,:)) & ~squeeze(f_correct_thresholded(methind,2,s,:)));
        right_wrong_not_over = sum(~squeeze(f_correct_thresholded(methind,1,s,:)) & squeeze(f_correct_thresholded(methind,2,s,:)));
        over_wrong_not_right = sum(squeeze(f_correct_thresholded(methind,1,s,:)) & ~squeeze(f_correct_thresholded(methind,2,s,:)));
        both_right = sum(squeeze(f_correct_thresholded(methind,1,s,:)) & squeeze(f_correct_thresholded(methind,2,s,:)));
        McNemarextest([both_wrong,right_wrong_not_over,over_wrong_not_right,both_right],2,0.05);

        fprintf('SNR = %df, Sim patch size 10mm, Reconstruct 10mm-Reconstruct 5mm\n',SNRs(s));
        both_wrong = sum(~squeeze(f_correct_thresholded(methind,4,s,:)) & ~squeeze(f_correct_thresholded(methind,3,s,:)));
        right_wrong_not_under = sum(~squeeze(f_correct_thresholded(methind,4,s,:)) & squeeze(f_correct_thresholded(methind,3,s,:)));
        under_wrong_not_right = sum(squeeze(f_correct_thresholded(methind,4,s,:)) & ~squeeze(f_correct_thresholded(methind,3,s,:)));
        both_right = sum(squeeze(f_correct_thresholded(methind,4,s,:)) & squeeze(f_correct_thresholded(methind,3,s,:)));
        McNemarextest([both_wrong,right_wrong_not_under,under_wrong_not_right,both_right],2,0.05);
    end

    disp('Comparing classifiers: pial');
    for s = 1:length(SNRs)
        fprintf('SNR = %df, Sim patch size 5mm, Reconstruct 5mm-Reconstruct 10mm\n',SNRs(s));
        both_wrong = sum(~squeeze(f_pial_thresholded(methind,1,s,:)) & ~squeeze(f_pial_thresholded(methind,2,s,:)));
        right_wrong_not_over = sum(~squeeze(f_pial_thresholded(methind,1,s,:)) & squeeze(f_pial_thresholded(methind,2,s,:)));
        over_wrong_not_right = sum(squeeze(f_pial_thresholded(methind,1,s,:)) & ~squeeze(f_pial_thresholded(methind,2,s,:)));
        both_right = sum(squeeze(f_pial_thresholded(methind,1,s,:)) & squeeze(f_pial_thresholded(methind,2,s,:)));
        McNemarextest([both_wrong,right_wrong_not_over,over_wrong_not_right,both_right],2,0.05);

        fprintf('SNR = %df, Sim patch size 10mm, Reconstruct 10mm-Reconstruct 5mm\n',SNRs(s));
        both_wrong = sum(~squeeze(f_pial_thresholded(methind,4,s,:)) & ~squeeze(f_pial_thresholded(methind,3,s,:)));
        right_wrong_not_under = sum(~squeeze(f_pial_thresholded(methind,4,s,:)) & squeeze(f_pial_thresholded(methind,3,s,:)));
        under_wrong_not_right = sum(squeeze(f_pial_thresholded(methind,4,s,:)) & ~squeeze(f_pial_thresholded(methind,3,s,:)));
        both_right = sum(squeeze(f_pial_thresholded(methind,4,s,:)) & squeeze(f_pial_thresholded(methind,3,s,:)));
        McNemarextest([both_wrong,right_wrong_not_under,under_wrong_not_right,both_right],2,0.05);
    end
end
