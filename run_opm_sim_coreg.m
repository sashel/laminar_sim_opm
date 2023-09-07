function [] = run_opm_sim_coreg(stem_dir, space, axis, SNRs, coreg_err)
% RUN_OPM_SIM_COREG Run laminar simulations for the data in project folder STEM_DIR,
% an OPM array with SPACE inter-sensor distance(s) and AXIS number of measurement axes 
% across the signal-to-noise ratios in SNRs for the coregistration errors in mm
% specified in COREG_ERR
%
% Use as run_opm_sim_coreg('/data/pt_np-helbling/layer_opm_sim/',55,1,[-5,-10,-20,-30,-40],1:4)
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes, the forth the SNRs in dB and
% the fifth the coregistrtaion errors in mm

for i = space
    for j = axis
        for k = SNRs
            for m = coreg_err

                simlayer_roi_opm_sim_coreg([],[], [10 30], SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data/sim_opm_space_%d_axis_%d.mat',space(i),axis(j))),...
                    'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_coreg_%d/',space(i),axis(j),coreg_err(m))),'coreg_err',coreg_err(m))
                fprintf('roi:/results_opm_sim_space_%d_axis_%d_coreg_%d\n',space(i),axis(j),coreg_err(m))

                simlayer_free_energy_opm_sim_coreg([],[], [10 30], SNRs(k),'rawfile',fullfile(stem_dir,sprintf('opm_sim_data/sim_opm_space_%d_axis_%d.mat',space(i),axis(j))),...
                    'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_coreg_%d/',space(i),axis(j),coreg_err(m))),'coreg_err',coreg_err(m))
                fprintf('free:/results_opm_sim_space_%d_axis_%d_coreg_%d\n',space(i),axis(j),coreg_err(m))

                plot_wholebrain_roi_opm_sim_results([], [], SNRs(k), 'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_coreg_%d/',space(i),axis(j),coreg_err(m))),'prefix', sprintf('space_%d_axis_%d_coreg_%d',space(i),axis(j),coreg_err(m)),'recompute',true)

            end
        end
    end
end

