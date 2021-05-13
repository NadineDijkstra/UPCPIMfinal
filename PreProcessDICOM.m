function PreProcessDICOM(cfg)
% function PreProcessDICOM(cfg)
%
% reads in DICOM files and performs Realignment, Coregistration and
% Segmentation in that order. 'cfg' is a configuration structure with the
% following fields:
% cfg.
%     dicom_dir = directory of the dicom files
%     nifti_dir = where the niftis should be written
%     run_nr = optional, which runs to include


[~,subjectname] = fileparts(fileparts(cfg.nifti_dir));

if ~exist(cfg.nifti_dir,'dir') % if the map doesn't exist, make it
    mkdir(cfg.nifti_dir);
end

if ~exist(fullfile(cfg.nifti_dir,'dicom_conversion_done'),'file') % check if converted dicoms do not already exist
    %% Dicom import
    dicoms = [];%str2fullfile(cfg.dicom_dir,'*.IMA');
    folders = dir(cfg.dicom_dir); folders = {folders.name};
    
    if iscell(cfg.dicom_dir)
        for s = 1:size(cfg.dicom_dir,2)
            if isfield(cfg,'run_nr') && ~isempty(cfg.run_nr{s})
                for r = 1:size(cfg.run_nr{s},2)
                    dicoms = [dicoms,str2fullfile([cfg.dicom_dir{s} '/' folders{strncmp(folders,sprintf('%03d',1),3)}],'*.IMA')];%sprintf('*SKYRA.%04d*.IMA',cfg.run_nr{s}(r)))];
                end
            else
                dicoms = [dicoms,str2fullfile(cfg.dicom_dir{s},'*.IMA')];
            end
            
        end
    else
        if isfield(cfg,'run_nr') && ~isempty(cfg.run_nr)
            for r = 1:size(cfg.run_nr,2)
                dicoms = [dicoms,str2fullfile([cfg.dicom_dir '/' folders{strncmp(folders,sprintf('%03d',cfg.run_nr(r)),3)}],'*.IMA')];%[dicoms,str2fullfile(cfg.dicom_dir,sprintf('*SKYRA.%04d*.IMA',cfg.run_nr(r)))];
            end
        else
            dicoms = [dicoms;str2fullfile(cfg.dicom_dir,sprintf('*%s*.IMA',subjectname))];
        end
    end
    
    dicom{1}.spm.util.dicom.data = dicoms';
    dicom{1}.spm.util.dicom.root = 'flat';
    dicom{1}.spm.util.dicom.outdir = {cfg.nifti_dir};
    dicom{1}.spm.util.dicom.convopts.format = 'nii';
    dicom{1}.spm.util.dicom.convopts.icedims = 0;
    
    fprintf('Running dicom import for subject %s\n',subjectname)    

    spm_jobman('run',dicom)
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'dicom_conversion_done')))
    clear dicoms
else
    fprintf('dicom conversion for subject %s already done\n',subjectname)
end

%% Align the structural to a T1 template
if ~exist(fullfile(cfg.nifti_dir,'align_to_t1_done'),'file') % check if converted dicoms do not already exist
    
    structural = str2fullfile(cfg.nifti_dir, 's3*192-01*.nii'); %
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'/vol/ccnlab1/naddij/Templates/MNI/mni_t1.nii'};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {structural};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    spm_jobman('run',matlabbatch)
        
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'align_to_t1_done')));
else
    fprintf('Aligning to T1 for subject %s already done\n',subjectname)
end    

%% Realignment

if ~exist(fullfile(cfg.nifti_dir,'realignment_done'),'file')
    
    niftis = str2fullfile(cfg.nifti_dir,'f3*.nii'); % functional nifti's
    
    realign{1}.spm.spatial.realign.estwrite.data = {niftis'};
    realign{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realign{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realign{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realign{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    realign{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    realign{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realign{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    realign{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realign{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realign{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realign{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realign{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    fprintf('running realignment for subject %s\n',subjectname)
    
    spm_jobman('run',realign)
    clear niftis
    
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'realignment_done')))
    
else
    fprintf('realignment for subject %s already done \n',subjectname)
end

%% Coregistration
if ~exist(fullfile(cfg.nifti_dir,'coregistration_done'),'file')
    
    mean_file = str2fullfile(cfg.nifti_dir,'meanf*.nii');
    mean_coreg_file = fullfile(cfg.nifti_dir,'mean_coreg.nii');
    
    % copy mean file to be used for coregistration estimation
    if ~exist(mean_coreg_file, 'file')
        copyfile(mean_file,mean_coreg_file);
    end
    
    real_scans = str2fullfile(cfg.nifti_dir,'rf*.nii');
    structural = str2fullfile(cfg.nifti_dir, 'rs*192-01.nii');%'s*192-01*.nii'); %
    
    coreg{1}.spm.spatial.coreg.estimate.ref = {structural};
    coreg{1}.spm.spatial.coreg.estimate.source = {mean_coreg_file};
    coreg{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    coreg{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    coreg{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    coreg{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    fprintf('running coregistration for subject %s\n',subjectname)
    
    spm_jobman('run',coreg)
    
    % SPM coregistration parameters
    flags.estimate.cost_fun = 'nmi';
    flags.estimate.sep = [4 2];
    flags.estimate.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    flags.estimate.fwhm = [7 7];
    
    % create transformation matrix
    coregParameters = spm_coreg(spm_vol(mean_coreg_file), ...
        spm_vol(mean_file), ...
        flags.estimate);
    EPItransformationMatrix = spm_matrix(coregParameters);
    
    % apply transformation matrix to functionals
    
    % add mean functional to the list
    real_scans{end+1} = mean_file;
    
    for iFunc = 1:size(real_scans,2)
        if mod(iFunc,100)==0
            fprintf('\tworking on file %d\n',iFunc)
        end
        iFuncSpace = spm_get_space(char(real_scans(iFunc)));
        spm_get_space(char(real_scans(iFunc)), EPItransformationMatrix\iFuncSpace);
    end
    
    unix(sprintf('touch %s',fullfile(cfg.nifti_dir,'coregistration_done')))
    clear real_scans
else
    fprintf('coregistration for subjects %s already done \n',subjectname)
    
end

%% Segmentation
if ~exist(fullfile(cfg.nifti_dir,'segmentation_done'),'file')
    
    structural = str2fullfile(cfg.nifti_dir, 'rs*192-01.nii');%'s*192-01*.nii'); %
    
    segm{1}.spm.spatial.preproc.channel.vols = {structural};
    segm{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    segm{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    segm{1}.spm.spatial.preproc.channel.write = [0 0];
    segm{1}.spm.spatial.preproc.tissue(1).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,1'};
    segm{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    segm{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    segm{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    segm{1}.spm.spatial.preproc.tissue(2).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,2'};
    segm{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    segm{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    segm{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    segm{1}.spm.spatial.preproc.tissue(3).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,3'};
    segm{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    segm{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    segm{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    segm{1}.spm.spatial.preproc.tissue(4).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,4'};
    segm{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    segm{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    segm{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    segm{1}.spm.spatial.preproc.tissue(5).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,5'};
    segm{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    segm{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    segm{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    segm{1}.spm.spatial.preproc.tissue(6).tpm = {'/vol/optdcc/spm12/tpm/TPM.nii,6'};
    segm{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    segm{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    segm{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    segm{1}.spm.spatial.preproc.warp.mrf = 1;
    segm{1}.spm.spatial.preproc.warp.cleanup = 1;
    segm{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    segm{1}.spm.spatial.preproc.warp.affreg = 'mni';
    segm{1}.spm.spatial.preproc.warp.fwhm = 0;
    segm{1}.spm.spatial.preproc.warp.samp = 3;
    segm{1}.spm.spatial.preproc.warp.write = [1 1];
    
    fprintf('running segmentation for subject %s\n',subjectname)
    spm_jobman('run',segm)   

    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'segmentation_done')))
else
    fprintf('segmentation of subject %s already done \n',subjectname)
    
end

%% Reslice key files
if ~exist(fullfile(cfg.nifti_dir,'reslicing_done'),'file')
    
    resliceParameters = struct(... 
        'prefix', 'r',...
        'mask', 1,...
        'interp', 4, ...
        'wrap', [0 0 0], ...
        'which', [2 1]);
    
    resFlags = struct(...
        'interp', resliceParameters.interp,... % interpolation type
        'wrap', resliceParameters.wrap,...     % wrapping info (ignore...)
        'mask', resliceParameters.mask,...     % masking (see spm_reslice)
        'which', 1,...                  % what images to reslice
        'mean', 0);                     % write mean image
    
    mean_file = str2fullfile(cfg.nifti_dir,'meanf*.nii');
    
    spm_reslice({mean_file, str2fullfile(cfg.nifti_dir,'c1*.nii')}, resFlags);
    spm_reslice({mean_file, str2fullfile(cfg.nifti_dir,'c2*.nii')}, resFlags);
    spm_reslice({mean_file, str2fullfile(cfg.nifti_dir,'c3*.nii')}, resFlags);
    
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'reslicing_done')))
else
    fprintf('reslicing of subject %s already done \n',subjectname)
    
end

