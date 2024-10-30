function [] = run_opm_sim_BEM(stem_dir, space, axis, SNRs, recompute_BEM)
% RUN_OPM_SIM_BEM Run laminar simulations for the project folder STEM_DIR, for an OPM array with SPACE
% With project folder STEM_DIR, SPACE the inter-sensor distance(s) of the OPM array,
% AXIS number of measurement axes, signal-to-noise ratios SNRs and
% RECOMPUTE_BEM a logical flag indicating whether to recompute the BEM


% Use as run_opm_sim_BEM('/data/pt_np-helbling/layer_opm_sim/',55,1,-10)
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes and the forth the SNRs

for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy_BEM(SNRs(k),recompute_BEM,'rawfile',fullfile(stem_dir,sprintf('opm_sim_data_BEM/sim_opm_custom_space_%d_axis_%d_BEM.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_BEM/',space(i),axis(j))))
            fprintf('free:/results_opm_sim_space_%d_axis_%d_BEM\n',space(i),axis(j))

            % run laminar simulations for the t-statistic ROI analysis
            simlayer_roi_BEM(SNRs(k),recompute_BEM,'rawfile',fullfile(stem_dir,sprintf('opm_sim_data_BEM/sim_opm_custom_space_%d_axis_%d_BEM.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_BEM/',space(i),axis(j))))
            fprintf('roi:/results_opm_sim_space_%d_axis_%d_BEM\n',space(i),axis(j))

            plot_wholebrain_roi_sim_results_sashel(SNRs(k),'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_BEM/',space(i),axis(j))),'prefix', sprintf('space_%d_axis_%d',space(i),axis(j)),'recompute',true)
        end
    end
end
end
