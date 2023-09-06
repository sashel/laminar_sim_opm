function [] = run_opm_sim(stem_dir, space, axis, SNRs)
% RUN_OPM_SIM Run laminar simulations for the project folder STEM_DIR,for an OPM array with SPACE
% inter-sensor distance(s) and AXIS number of measurement axes across the
% signal-to-noise ratios in SNRs

% Use as run_opm_sim('/data/pt_np-helbling/layer_opm_sim/',1,[-5,-10,-20,-30,-40])
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes and the forth the SNRs

for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the t-statistic ROI analysis
            simlayer_roi(SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data/sim_opm_custom_space_%d_axis_%d.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_test/',space(i),axis(j))))
            fprintf('roi:/results_opm_sim_space_%d_axis_%d\n',space(i),axis(j))

            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy_opm(SNRs(k), 'rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d/',space(i),axis(j))))
            fprintf('free:/results_opm_sim_space_%d_axis_%d\n',space(i),axis(j))
           
            plot_wholebrain_roi_sim_results(SNRs(k),'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d/',space(i),axis(j))),'prefix', sprintf('space_%d_axis_%d',space(i),axis(j)),'recompute',true)
        end
    end
end
end
