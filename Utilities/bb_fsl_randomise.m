function bb_fsl_randomise(cfg)
% bb_fsl_randomise(cfg)
%
% Uses FSL machinary to perform non-parametric statistics.
%
% cfg.inputfiles        = cell array with image filenames (e.g. single subject
%                           betas or contrasts)
% cfg.outputroot        = string containing output file-rootname - provide
%                           without extension e.g. '/home/fsl_randomise_'
%
% cfg.design            = design matrix, [images-by-experimental variables]
%                           see: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM
% cfg.contrasts         = contrast matrix, [contrasts-by-experimental variables]
% cfg.contrastType      = 't' or 'f';
% cfg.exchangeability   = exchangeability-block column vector, [images-by-1]
% cfg.mask              = stats mask
% cfg.options           = randomise options
%                           Default: 'nopermute' (no permutation)
%                           Built-in option: 'tcfe' (Threshold-Free Cluster Enhancement
%                           with 5mm sd variance smoothing to T-stats)
%                           Other: find out by typing 'randomise' in
%                           terminal (example of other options:  '-T -v 5 -n 10000')

% load variables from cfg
get_vars_from_struct(cfg)


%%

if ~isfield(cfg,'options') || isempty(cfg.options)
    options = 'nopermute';
    fprintf('Using default %s option for fast second level t-stat image\n',options)
end

switch options
    case 'nopermute'
        options = '-1 -n 1';
    case 'tfce'
        options = '-1 -T -v 5 -n 10000';
end

[outputPath,~,~] = fileparts(outputroot);

if isfield(cfg,'design') && ~isempty(cfg.design)
    useDesign = true;
    % construct design matrix file
    matFn  = fullfile(outputPath,'design.mat');
    mat2fsltxt(design,matFn)
    
    if ~isfield(cfg,'contrasts') || isempty(cfg.design)
        error('Specify contrasts in accordance with design matrix')
    else
        % construct contrasts file
        if ~isfield(cfg,'contrastType') && ~isempty(cfg.contrastType)
            error('Specify contrasts type. t or f contrast?')
        else
            switch contrastType
                case 't'
                    useTcon = true;
                    useFcon = false;
                    conFn  = fullfile(outputPath,'design.con');
                case 'f'
                    useFcon = true;
                    useTcon = false;
                    conFn  = fullfile(outputPath,'design.fts');
                otherwise
                    error('Unknown contrasts type. t or f contrast?')
            end
        end
        
        mat2fsltxt(contrasts,conFn)
        
    end
    
else
    useDesign = false;
    useTcon = false;
    useFcon = false;
end

if isfield(cfg,'exchangeability') && ~isempty(cfg.exchangeability)
    useGrp = true;
    % construct exchangeability-block file
    grpFn  = fullfile(outputPath,'design.grp');
    mat2fsltxt(exchangeability,grpFn)
else
    useGrp = false;
end

if isfield(cfg,'mask') && ~isempty(cfg.mask)
%     maskImgFn = mask;
    statsMask = true;
else
    statsMask = false;
end

%% setup

% merge all images to 4D NIFTI
mergedDataFn    = [outputroot 'merged_data.nii.gz'];
cfg             = [];
cfg.inputfiles  = inputfiles;
cfg.outputfile  = mergedDataFn;
bb_fsl_merge(cfg);

% load images
[V, Y] = read_nii(mergedDataFn);

if statsMask
    copyfile(mask,[outputroot 'mask.nii.gz']);
    maskImgFn   = [outputroot 'mask.nii.gz'];
    V3d = V;
    V3d.dim(4) = [];
else
    % construct binary mask that only includes voxels that are not zeros across
    % all images
    maskImgFn   = [outputroot 'mask.nii.gz'];
    maskImg     = all(Y,4);
    V3d = V;
    V3d.dim(4) = [];
    write_nii(V3d,maskImg,maskImgFn);
end

% write average image to disk
meanImgFn   = [outputroot 'mean.nii.gz'];
meanImg     = mean(Y,4);
write_nii(V3d,meanImg,meanImgFn);
clear Y

% construct fsl command
fslCommand = ['randomise -i ' mergedDataFn ' -m ' maskImgFn ...
    ' -o ' outputroot ' ' options];

if useDesign
    fslCommand = [fslCommand ' -d ' matFn ' '];
end
if useTcon
    fslCommand = [fslCommand ' -t ' conFn ' '];
end
if useFcon
    fslCommand = [fslCommand ' -f ' conFn ' '];
end
if useGrp
    fslCommand = [fslCommand ' -e ' grpFn ' '];
end

% execute fsl command
unix(fslCommand);


%% cleanup
delete(mergedDataFn);
%delete(maskImgFn);

% if useDesign
%     delete(matFn)
%     delete(conFn)
% end
%
% if useGrp
%     delete(grpFn)
% end
