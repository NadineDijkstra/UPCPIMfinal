function MakeContrasts(cfg)
% function MakeContrasts(cfg)

spm_file = str2fullfile(fullfile(cfg.root,cfg.subjectID,cfg.model_dir),...
    'SPM.mat');

matlabbatch{1}.spm.stats.con.spmmat = {spm_file}';

% t contrasts
if isfield(cfg,'tcontrasts')
for con = 1:length(cfg.tcontrasts)
    
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name = cfg.tnames{con};
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = cfg.tcontrasts{con};
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'none';
end
count = con;
else 
    count = 0;
end

% f contrasts
if isfield(cfg,'fcontrasts')
for con = 1:length(cfg.fcontrasts)
    
    matlabbatch{1}.spm.stats.con.consess{count+con}.fcon.name = cfg.fnames{con};
    matlabbatch{1}.spm.stats.con.consess{count+con}.fcon.weights = cfg.fcontrasts{con};
    matlabbatch{1}.spm.stats.con.consess{count+con}.fcon.sessrep = 'none';
end
end

matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run',matlabbatch)