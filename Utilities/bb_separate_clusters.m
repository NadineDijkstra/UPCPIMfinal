function bb_separate_clusters(cfg)
% bb_separate_clusters(cfg)
%
% This function distinguishes clusters within a statistics image and saves
% them in separate nii-files.
%
% cfg.inputfile       = statisical image (e.g. t-map)
% cfg.threshold       = t-value

% load variables from cfg
get_vars_from_struct(cfg)


%%

[root, name, ext] = fileparts(inputfile);
name = strtok(name,'.');
outputFn = fullfile(root, [name '_clusters.nii.gz']);

fslCmd = sprintf('cluster -i %s --thresh=%f -o %s --mm',inputfile,threshold,outputFn);

[~,clusterOutput] = unix(fslCmd);
fprintf(clusterOutput)

% extract and save clusters
[V, Y] = read_nii(outputFn);
nClusters = max(unique(Y));

for iCluster = 1:nClusters
    
    iCluster_Y = Y==iCluster;
    if strcmp(ext,'.nii')
        write_nii(V,iCluster_Y, fullfile(root, [name '_cluster' num2str(iCluster) ext]));
    elseif strcmp(ext,'.gz')
        write_nii(V,iCluster_Y, fullfile(root, [name '_cluster' num2str(iCluster) '.nii' ext]));
    end
    
end


