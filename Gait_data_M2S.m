% This script first generates micro-Doppler signatures Combining multiple raw
% files, then the sx matrix can be used to get the torso for legs/hands
% reflections.

clear; clc; close all; warning off;
data = ['/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Data_21 jan_Top/*.bin'];
mDout = '/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Data_21 jan_Top/';


files = dir(data);

seqPerRecord = 1; 

filenames2 = {files.name};
for z = 1:length(filenames2)
        temp{1,z} = filenames2{z}(1:end-10);
end
uniqs = unique(temp);

for j = 1:length(uniqs) % 12
        match = strfind(filenames2,uniqs{j}); % find matches
        idx = find(~cellfun(@isempty,match)); % find non-empty indices
%         if j >=12
%                 idx = max(idx);
%         end
        RDC = [];
        % concat RDCs with same names
        for r = 1:length(idx)
                fname = fullfile(files(idx(r)).folder,files(idx(r)).name);
                temp2 = RDC_extract_2243(fname);
                RDC = [RDC temp2];
        end
        % divide into sub RDCs
        numChirps = floor(size(RDC,2)/seqPerRecord);
        for r =1:seqPerRecord
                tic

                subRDC = RDC(:,(r-1)*numChirps+1:r*numChirps,:);
               sx=RDC_to_sx_2243( subRDC);
                torso=find_torso(sx);
                
        end
        final_torso=prf_idx2vel(torso);
        matname=[fname(1:end-4) '.mat'];
        save(matname,'final_torso');
end





 


