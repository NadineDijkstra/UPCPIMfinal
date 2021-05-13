function clusterTable = bb_cluster_table(cfg)
% clusterTable = bb_cluster_table(cfg)
%
% cfg.inputfile       = statisical image (e.g. t-map)
% cfg.pImage          = image containing (un)corrected p-values, usually in
%                           the 1-p format for easy thresholding
% cfg.threshold       = t-value or (1 minus) p-value (in this case also 
%                           specify cfg.pImage)
% cfg.df              = integer, degrees of freedom (e.g. N-1) for conversion
%                          T-to-P-toZ conversion (only for parametric
%                          stats!)
% cfg.atlas           = string, can be 'aal', 'afni' or a filename of a
%                           atlas file supported by ft_read_atlas
% cfg.searchRadius    = radius of sphere in which anatomical label may be
%                           searched
% cfg.outputfile      = tab delim txt file where the table is stored

% load variables from cfg
get_vars_from_struct(cfg)


if ~isfield(cfg,'atlas') || isempty(cfg.atlas)
    findLabels = false;
else
    findLabels = true;
end 

if ~isfield(cfg,'searchRadius') || isempty(cfg.searchRadius)
    searchRadius = 0;
end

if isfield(cfg,'outputfile') && ~isempty(cfg.outputfile) 
    writeOutput = true;
else
    writeOutput = false;
end

if isfield(cfg,'pImage') && ~isempty(cfg.pImage) 
    
else

end



%%

fid = fopen(clusterTable);
cellTable        = textscan(fid,'%s %s %s %s %s %s %s %s %s');
cellTable        = horzcat(cellTable{:});
cellTable(1:3,:) = [];
cellTable(:,7:9) = [];

tVals  = cellfun(@str2num,cellTable(:,3));

%[pVals zVals] = t2p2z(tVals,df);
pVals = tVals; zVals = tVals;

mniCoord = cellfun(@str2num,cellTable(:,4:6));


%% reorder table into following format:
% 1: cluster number
% 2: number of voxels
% 3: X coord mni
% 4: Y coord mni
% 5: Z coord mni
% 6: peakZ
% 7: peakT
% 8: peakP
% 9: anatomical labels

clusterTable = cellTable(:,[1 2 4 5 6 3]);

nRois = size(clusterTable,1);

if findLabels
    
    switch atlas
        case 'aal'
            atlasFn = 'D:\Documents\MATLAB\Analyses\FieldTrip\fieldtrip-20191125\fieldtrip-20191125\template\atlas\aal\ROI_MNI_V4.nii';
        case 'afni'
            atlasFn = 'D:\Documents\MATLAB\Analyses\FieldTrip\fieldtrip-20191125\fieldtrip-20191125\template\atlas\afni\TTatlas+tlrc.BRIK';
        otherwise
            atlasFn = atlas;
    end
    
    addpath('D:\Documents\MATLAB\Analyses\FieldTrip\fieldtrip-20191125\fieldtrip-20191125')
    ft_defaults;
    atlasStruct = ft_read_atlas(atlasFn);
    
    labels = atlas_label_from_coord(mniCoord,atlasStruct.tissue,atlasStruct.tissuelabel,atlasStruct.transform,searchRadius);
else
    labels = cell(nRois,1); % empty cells
end


for iRoi = 1:nRois
    
    clusterTable{iRoi,6} = sprintf('%.2f',zVals(iRoi));
    
    if ~isempty(labels{iRoi})
        
        nRegions = length(labels{iRoi}{1});
        clusterTable(iRoi,7:7+nRegions-1) = labels{iRoi}{1};
    else
        clusterTable{iRoi,7} = '';
    end
end

if length(clusterTable) < 10
    disp(clusterTable);
else
    disp(clusterTable(1:10,:))
end

%% write output to csvfile

if writeOutput
    
    fid = fopen(outputfile,'w');
    
    fprintf(fid,'clusterId\tnVox\tmniX\tmniY\tmniZ\tpeakZ\tpeakT\tpeakP\tlabel\n');
    
    for iRoi = 1:nRois
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',clusterTable{iRoi,:});
    end
    fclose(fid);
    
end

