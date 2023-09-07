function [] = run_opm_sim_offset(stem_dir, space, axis, SNRs, offsets)
% RUN_OPM_SIM_OFFSET Run laminar simulations for the data in project folder STEM_DIR, an OPM array 
% with SPACE inter-sensor distance(s) and AXIS number of measurement axes across the
% signal-to-noise ratios in SNRs for the scalp-sensor offsets specified in OFFSETS
%
% Use as run_opm_sim_offset('/data/pt_np-helbling/layer_opm_sim/',55,1,[-5,-10,-20,-30,-40],20)
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes, the forth the SNRs and the
% fifth the scalp-sensor offsets in mm

for i = space
    for j = axis
        for k = SNRs
            for m = offsets
                simlayer_roi_opm_sim([10 30], SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data/sim_opm_space_%d_axis_%d_offset_%d_ds.mat',space(i),axis(j),offsets(m))),...
                    'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_offset_%d/',space(i),axis(j),offsets(m))))
                fprintf('roi:/results_opm_sim_space_%d_axis_%d_offset_%d\n',space(i),axis(j),offsets(m))

                simlayer_free_energy_opm_sim([10,30], SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data/sim_opm_space_%d_axis_%d_offset_%d_ds.mat',space(i),axis(j),offset(m))),...
                    'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_offset_%d/',space(i),axis(j),offsets(m))))
                fprintf('free:/results_opm_sim_space_%d_axis_%d_offset_%d\n',space(i),axis(j),offsets(m))

                plot_wholebrain_roi_opm_sim_results(SNRs(k), 'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_offset_%d/',space(i),axis(j),offsets(m))),'prefix', sprintf('space_%d_axis_%d_offset_%d',space(i),axis(j),offsets(m)),'recompute',true)
            end
        end
    end
end
