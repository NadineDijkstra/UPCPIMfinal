function bb_fsl_merge(cfg)
% bb_fsl_merge(cfg)
%
% Merges images from input array into 4D Nifti. Preferable input cell array
% should contain full filenames e.g. {'/home/file1.nii','/home/file2.nii'}
% and the output filename outputFn should also contain a full path.
%
% cfg.inputfiles    = cell array with files to be merged
% cfg.outputfile    = output filename of merged nifti

% load variables from cfg
get_vars_from_struct(cfg)


%%

try
    concatImgStr = sprintf('%s ', inputfiles{:});
catch
    concatImgStr = inputfiles;
end

fslCommand = ['fslmerge -t ' outputfile ' ' concatImgStr];

unix(fslCommand);