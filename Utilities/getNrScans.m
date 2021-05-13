function [nrscans] = getNrScans(real_niftis,nblocks)

nrscans = zeros(nblocks,1); % usually there is one run per block

end_code = strfind(real_niftis{1,1},'-0');
code = real_niftis{1,1}(1:end_code(1,1));

for run = 1:nblocks
    if run == 1
        nrscans(run,1) = size(find(strncmp([code real_niftis{1,1}(end_code(1,1)+1:end_code(1,1)+4)],real_niftis',size([code real_niftis{1,1}(end_code(1,1)+1:end_code(1,1)+4)],2))),1);
    else
        if sum(nrscans(1:run-1,1))+1 < size(real_niftis,2)
            nrscans(run,1) = size(find(strncmp([code real_niftis{1,sum(nrscans(1:run-1,1))+1}(end_code(1,1)+1:end_code(1,1)+4)],real_niftis',size([code real_niftis{1,sum(nrscans(1:run-1,1))+1}(end_code(1,1)+1:end_code(1,1)+4)],2))),1);
        else
            nrscans(run,1) = 0;
        end
    end
end

% fill 0 as number of scans for block if it was included in previous run
two_blocks = find(nrscans > 400); % if the number of scans was more that 400, then there were two blocks
nrscans(two_blocks+1:nblocks) = [0;nrscans(two_blocks+1:nblocks-1)];
nrscans  = nrscans(nrscans~=0);
