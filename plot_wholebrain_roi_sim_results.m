function plot_wholebrain_roi_sim_results(snr, varargin)
% PLOT_WHOLEBRAIN_ROI_SIM_RESULTS  Plot free energy and t stastics for all
% simulations
%
% Use as
%   plot_wholebrain_roi_sim_results(-20), where the argument is the SNR (db).
% 
%   plot_wholebrain_roi_sim_results(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * surf_dir - directory containing subject surfaces
%    * sim_patch_size - 0 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 0 (default) or interger - reconstruction patch size

% parse inputs
addpath('/data/pt_02058/megdata/speech/analysis/software/brewermap')
cm = colormap(brewermap(8, 'Set2'));

defaults = struct('surf_dir','<FS_DIR>' , 'mri_dir', '<T1_DIR>',...
    'out_path', '', 'prefix', '', 'dipole_moment', 10, 'sim_patch_size', 0,...
    'reconstruct_patch_size', 0, 'nsims', 60, 'recompute',true,'freq',[10 30]);  % define default values

params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

if isempty(params.out_path)
    params.out_path = '/data/pt_np-helbling/layer_opm_sim/results_opm_sim_space_35_axis_1';
end

methodnames={'EBB','MSP'}; % alternatively {'IID','COH'}
method_picks = [1,2];
Nmeth = length(method_picks);

% define original and downsampled white matter surface
orig_white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','white.hc_PDw2.surf.gii');
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');

orig_pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','pial.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

allmeshes=strvcat(white_mesh,pial_mesh);
Nmesh=size(allmeshes,1);

wholeBrainPialWhiteF=zeros(Nmeth,Nmesh,params.nsims);
% load free energy results
freq = params.freq;

if params.sim_patch_size==0 && params.reconstruct_patch_size==0
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
        freq(1),freq(2),snr,params.dipole_moment);
    data_file = fullfile(params.out_path, fname);
else
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',...
        freq(1),freq(2),snr,params.dipole_moment, params.sim_patch_size,...
        params.reconstruct_patch_size);
    data_file = fullfile(params.out_path, fname);
end

roiWhitePialT=zeros(Nmeth,Nmesh*params.nsims);
data_dir=fullfile(params.out_path, sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1),...
    freq(2), snr, params.dipole_moment));

for methind=1:Nmeth   
    method=methodnames{method_picks(methind)};
    roiWhitePialT(methind,:)=get_wmpial_t(data_dir, method, params.nsims, deblank(allmeshes(2,:)), ...
        allmeshes(1,:), orig_pial_mesh, orig_white_mesh, 'recompute_trials',params.recompute,'delete_giftis',false);
end

figure('Position',[1 1 800 800]);

mesh_idx=zeros(Nmeth,Nmesh);
mesh_idx(1,1)=3;
mesh_idx(1,2)=1;
mesh_idx(2,1)=7;
mesh_idx(2,2)=5;
mesh_idx(3,1)=11;
mesh_idx(3,2)=9;
mesh_idx(4,1)=15;
mesh_idx(4,2)=13;
for methind=1:Nmeth       
    for simmeshind=1:Nmesh    
        subplot(Nmeth*Nmesh,2,mesh_idx(methind,simmeshind));
        [~,file,~]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        hold on
        for simind=1:params.nsims
            color = cm(7,:);
            if wholeBrainPialWhiteF(method_picks(methind),simmeshind,simind)<0
                color = cm(2,:); 
            end
            bar(simind,wholeBrainPialWhiteF(method_picks(methind),simmeshind,simind),'FaceColor',color,'EdgeColor','none');
        end
        xlim([0 params.nsims+1]);
        if (methind == Nmeth && simmeshind == Nmesh-1)
            xlabel('Simulation')
        end
        ylabel('\Delta F');      
        title(sprintf('%s - %s',methodnames{method_picks(methind)},simmeshname));        
    end

end

mesh_idx=zeros(Nmeth,Nmesh);
mesh_idx(1,1)=4;
mesh_idx(1,2)=2;
mesh_idx(2,1)=8;
mesh_idx(2,2)=6;
mesh_idx(3,1)=12;
mesh_idx(3,2)=10;
mesh_idx(4,1)=16;
mesh_idx(4,2)=14;

% dof=514;
dof = 199;
alpha=1.0-(0.05/2);
t_thresh=tinv(alpha, dof);

for methind=1:Nmeth,    
% for methind = method_picks,      
    method = methodnames{methind};
    
    % For each simulated mesh
    for simmeshind=1:Nmesh,
        [~,file,~]=fileparts(allmeshes(simmeshind,:));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        subplot(Nmeth*Nmesh,2,mesh_idx(methind,simmeshind));
        hold on
        for simind=1:params.nsims
            color = cm(7,:);
            if roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind)<0
                color = cm(2,:);
            end
            bar(simind,roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind),'FaceColor',color,'EdgeColor','none');
        end
        hold on

        xlim([0 params.nsims+1]);
        if (methind == Nmeth && simmeshind == Nmesh-1)
            xlabel('Simulation')
        end
        ylabel('ROI t-value');
        title(sprintf('%s - %s',methodnames{methind},simmeshname)); 
        
    end
end
set(gcf,'color','w')
print(sprintf('%s/%s_SNR%d_results',params.out_path,params.prefix,snr),'-dpdf','-bestfit')
save(sprintf('%s/%s_SNR%d_results',params.out_path,params.prefix,snr),'wholeBrainPialWhiteF','roiWhitePialT')
