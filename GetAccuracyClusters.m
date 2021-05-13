function GetAccuracyClusters(cfg)

workingDir = fullfile(cfg.root,'GroupResults',cfg.dir);

% load the accuracy map
[V,accuracy] = read_nii(fullfile(workingDir,cfg.map));
%accuracy = squeeze(accuracy(:,:,:,cfg.pair)); 
V.dim = size(accuracy);

% load the cluster map
[~,sigClusters] = read_nii(fullfile(workingDir,cfg.clusters));
sigClusters = sigClusters > 0;

% mask
sigAcc = accuracy; sigAcc(~sigClusters) = 0;

% write
[~,contrast] = fileparts(cfg.map);
[~,contrast] = fileparts(contrast); contrast = contrast(4:end);
name = fullfile(workingDir,sprintf('SigAcc_%d_%s.nii',cfg.pair,contrast));

write_nii(V,sigAcc,name);
