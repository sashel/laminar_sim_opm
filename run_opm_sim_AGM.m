function [] = run_opm_sim_AGM(stem_dir, space, axis, SNRs)
% RUN_OPM_SIM_AGM Run laminar simulations using an alternative generative
% model (AGM)for source reconstruction.
% With project folder STEM_DIR, SPACE the inter-sensor distance(s) of the OPM array,
% AXIS number of measurement axes and signal-to-noise ratios SNRs

% Use as run_opm_sim_AGM('/data/pt_np-helbling/layer_opm_sim/',55,1,[-5,-10,-20,-30,-40,-50])
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes and the forth the SNRs

for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the t-statistic ROI-based analysis
            simlayer_roi_AGM(SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d_ds.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_AGM/',space(i),axis(j))))

            fprintf('roi:/results_opm_sim_space_%d_axis_%d_AGM\n',space(i),axis(j))

            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy_AGM(SNRs(k), 'rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d_ds.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_AGM/',space(i),axis(j))))

            fprintf('free:/results_opm_sim_space_%d_axis_%d_AGM\n',space(i),axis(j))

            plot_wholebrain_roi_sim_results_ds11(SNRs(k),'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_AGM/',space(i),axis(j))),'prefix', sprintf('space_%d_axis_%d_ds',space(i),axis(j)),'recompute',true)

        end
    end
end
end

