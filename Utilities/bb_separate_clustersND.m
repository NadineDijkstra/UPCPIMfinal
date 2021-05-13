function bb_separate_clustersND(cfg)
% bb_separate_clustersND(cfg)
%
% This function distinguishes clusters within a statistics image and saves
% them in separate nii-files.
%
% cfg.inputfile       = statisical image (e.g. t-map)
% cfg.threshold       = t-value
% cfg.numVox          = minimum num voxels per cluster (default = 0)

% load variables from cfg
get_vars_from_struct(cfg)


%%

[root, name, ext] = fileparts(inputfile);
name = strtok(name,'.');
outputFn = fullfile(root, [name '_clusters.nii.gz']);

fslCmd = sprintf('cluster -i %s --thresh=%f -o %s --mm',inputfile,threshold,outputFn);

[~,clusterOutput] = unix(fslCmd);
fprintf(clusterOutput)

clusters = strsplit(clusterOutput(106:end-1));
clusters = reshape(clusters,[9,length(clusters)/9])';
numVoxels = str2double(clusters(:,2));
clusterIdx = str2double(clusters(:,1));
cIdx     = clusterIdx(numVoxels >= cfg.numVox);

% extract and save big enough clusters
[V, Y] = read_nii(outputFn);

Y(~ismember(Y,cIdx)) = 0;

for i = 1:length(cIdx)
    Y(Y==cIdx(i)) = i;
end

write_nii(V,Y,outputFn);



% for iCluster = 1:nClusters
%     
%     iCluster_Y = Y==iCluster;
%     if strcmp(ext,'.nii')
%         write_nii(V,iCluster_Y, fullfile(root, [name '_cluster' num2str(iCluster) ext]));
%     elseif strcmp(ext,'.gz')
%         write_nii(V,iCluster_Y, fullfile(root, [name '_cluster' num2str(iCluster) '.nii' ext]));
%     end
%     
% end


