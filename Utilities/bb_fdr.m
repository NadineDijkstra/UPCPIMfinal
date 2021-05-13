function fdr_threshold = bb_fdr(cfg)

% bb_fdr(cfg)
%
% Uses FSL's fdr function to determine FDR-corrected p-threshold
%
% cfg.inputfile         = input file (p-map)
% cfg.qvalue            = FDR threshold (default = 0.05)
% cfg.mask              = brain mask

% load variables from cfg
get_vars_from_struct(cfg)

%%

fslCommand = ['fdr -i ' inputfile ' -m ' mask ' -q ' num2str(qvalue)];

[~,fdr_threshold] = unix(fslCommand);

% fdr_threshold is now a char array
[~, fdr_threshold] = strtok(fdr_threshold,'0');
fdr_threshold = str2double(fdr_threshold);