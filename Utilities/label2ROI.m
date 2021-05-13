function status = label2ROI(subject, ROIs, hemis, funcTemplate, regFile, outDir)

% label2ROI(subject, ROIs, hemis, funcTemplate, regFile, outDir)
%
% Converts all label files for ROIs <ROIs> and hemispheres <hemis> for
% subject <subject> to 3D NIFTI binary masks, using <regFile> and
% <funcTemplate> for alignment.

status = 0;

if nargin < 6, outDir = './'; end
if nargin < 5,
    out('\nERROR: Not enough input arguments.\n')
    status = 1;
    return
end

[s, subjects_dir] = unix('echo $SUBJECTS_DIR');
if s==1
    out('\nERROR: Variable SUBJECTS_DIR does not seem to be defined. Exiting.\n');
    return
end

subjects_dir(end) = []; %Because for some reason, there's always a return "character" at the end
if strcmp(subjects_dir(end), filesep), subjects_dir(end) = []; end %Remove filesep character at the end if necessary

for iHemi = 1:length(hemis)
    for iROI = 1:length(ROIs);
        labelfile = [subjects_dir, filesep, subject, '/label/', hemis{iHemi}, '.', ROIs{iROI}, '.label'];
        if exist(labelfile, 'file')
            unixcmd = ['mri_label2vol --label ', labelfile, ' --temp ', funcTemplate, ' --subject ', subject, ' --hemi ', hemis{iHemi}, ...
                ' --o ', outDir, filesep, hemis{iHemi}, '_', ROIs{iROI}, '.nii --proj frac 0 1 0.1 --fillthresh 0.3 --reg ', regFile];
            [s,~] = unix(unixcmd, '-echo');
            if s==1
                out(['\nERROR: An error occured while processing label ', labelfile, '. Exiting.\n']);
                status = 1;
                return
            end
        end
    end %iROI
end



end