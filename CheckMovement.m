function CheckMovement(cfg)
% function CheckMovement(cfg)

movementFile = str2fullfile(cfg.nifti_dir,'rp*.txt');

fid = fopen(movementFile{1});
out = textscan(fid,'%f %f %f %f %f %f');
movement = cell2mat(out); clear out

% plot the things
[~,subName] = fileparts(fileparts(cfg.nifti_dir));
figure;
subplot(2,1,1);
plot(movement(:,1:3)); ylabel('mm'); xlabel('scans')
legend('X','Y','Z');
title(subName);
subplot(2,1,2);
plot(movement(:,4:6)); ylabel('mm'); xlabel('scans')
legend('pitch','roll','yaw');