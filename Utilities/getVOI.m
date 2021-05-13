function getVOI(spm_file,subjectname,fCon,tCon,pthresh,VOIname,VOIcentre)
% function getVOI(spm_file,fCon,tCon,pthresh,VOIname,VOIcentre)
% spm_file     : the SPM.mat of the GLM       
% fCon         : number f-contrast of effects of interest
% tCon         : number t-contrast to base peak on
% pthresh      : uncorrected threshold for the peak voxel
% VOIname      : name of the VOI
% VOIcentre    : centre of group effect around which to find the peak


matlabbatch{1}.spm.util.voi.spmmat = {spm_file};
matlabbatch{1}.spm.util.voi.adjust = fCon; % which F-contrast to adjust for
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = VOIname;
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {spm_file};
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = tCon;
matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = pthresh;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = VOIcentre;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 16;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = [0 0 0];
matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8;
matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';

qsubfeval('runBatch',matlabbatch,'timreq',60*60*1,'memreq',1024^3*4,'options','-V','batchid',sprintf('getVOI %s',subjectname))
