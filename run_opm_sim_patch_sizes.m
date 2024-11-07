function [] = run_opm_sim_patch_sizes(stem_dir, space, axis, SNRs)
% RUN_OPM_SIM_PATCH_SIZES Run laminar simulations with congruent and incongruent patch sizes
% using the whole-brain free energy analysis
% for the project folder STEM_DIR,fan OPM array with SPACE
% inter-sensor distance(s) and AXIS number of measurement axes across the
% signal-to-noise ratios in SNRs

% Use as run_opm_sim_patch_sizes('/data/pt_np-helbling/layer_opm_sim/',35,1,[-5,-10,-20,-30,-40])
% where the first argument is the project folder, the second the inter-sensor distance in mm,
% the third the number of measurement axes and the forth the SNRs
% Note that here we do not run simulations with congruent patches of 5 mm for simulated and reconstructed 
% sources as these simulations were performed using run_opm_sim

% sim patch: 10 mm, recon patch 5mm
for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy(SNRs(k), 'sim_patch_size', 10, rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_sim_patchsize_10/',space(i),axis(j))))
            fprintf('free:/results_opm_sim_space_%d_axis_%d_sim_patchsize_10\n',space(i),axis(j))
        end
    end
end

% recon patch: 10 mm, sim patch 5mm
for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy(SNRs(k), 'reconstruct_patch_size', 10, rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_recon_patchsize_10/',space(i),axis(j))))
            fprintf('free:/results_opm_sim_space_%d_axis_%d_recon_patchsize_10\n',space(i),axis(j))
        end
    end
end

% sim patch: 10 mm, recon patch 10 mm
for i = 1:length(space)
    for j = 1:length(axis)
        for k = 1:length(SNRs)
            % run laminar simulations for the whole-brain free energy analysis
            simlayer_free_energy(SNRs(k), 'sim_patch_size', 10, 'reconstruct_patch_size', 10,rawfile',fullfile(stem_dir,sprintf('opm_sim_data_spmdev/sim_opm_custom_space_%d_axis_%d.mat',space(i),axis(j))),...
                'out_path',fullfile(stem_dir,sprintf('/results_opm_sim_space_%d_axis_%d_sim_recon_patchsize_10/',space(i),axis(j))))
            fprintf('free:/results_opm_sim_space_%d_axis_%d_sim_recon_patchsize_10\n',space(i),axis(j))
        end
    end
endend
