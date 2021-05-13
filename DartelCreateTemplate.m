function DartelCreateTemplate(cfg)

% output dir
if ~exist(cfg.dartel_dir,'dir'); mkdir(cfg.dartel_dir); end

% collect segmentations and put in DARTEL folder
nsubjects = length(cfg.subjects);

for s = 1:nsubjects
    gm = str2fullfile(fullfile(cfg.root,cfg.subjects{s},cfg.nifti_dir),'c1*.nii');
    [~,name] = fileparts(gm);
    copyfile(gm,fullfile(cfg.dartel_dir,[name '.nii']));
    wm = str2fullfile(fullfile(cfg.root,cfg.subjects{s},cfg.nifti_dir),'c2*.nii');
    [~,name] = fileparts(wm);
    copyfile(wm,fullfile(cfg.dartel_dir,[name '.nii']));
    clear gm wm
end

% RUN DARTEL: Create templates
gm = str2fullfile(cfg.dartel_dir,'c1*.nii');
wm = str2fullfile(cfg.dartel_dir,'c2*.nii');

dt_create_templates{1}.spm.tools.dartel.warp.images = {gm',...
    wm'};
dt_create_templates{1}.spm.tools.dartel.warp.settings.template = 'Template';
dt_create_templates{1}.spm.tools.dartel.warp.settings.rform = 0;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
dt_create_templates{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
dt_create_templates{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
dt_create_templates{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
dt_create_templates{1}.spm.tools.dartel.warp.settings.optim.its = 3;

spm_jobman('run',dt_create_templates);

