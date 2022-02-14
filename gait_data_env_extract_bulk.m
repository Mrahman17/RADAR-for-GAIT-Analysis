clear; clc; close all; warning off;
DATA_DIR='/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Summer 2021 Data/77GHz/02_Normal Walk/';
OUT_DIR = '/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Data_Martelli_10sept/';   
% if 2~=exist(OUT_DIR,'dir')
%                 mkdir(OUT_DIR);
%  end
pattern = strcat(DATA_DIR, '*.bin');    % file pattern
files = dir(pattern);

% w = waitbar(0);

I_MAX = numel(files); % # of files in "files"

for i =1%: I_MAX  % for the first 20 iteration
        tic
        msg = strcat(['Processing file ', int2str(i), ' of ', int2str(I_MAX)]);   % loading message
        %     waitbar(i/I_MAX, w, msg);
        disp(msg);
        fName = files(i).name;
        [foo1, name, foo2] = fileparts(fName);
        fIn = fullfile(files(i).folder, files(i).name);
        fOut = strcat(OUT_DIR, name, '.mat');
        
        sx = get_sx_AWR1642_bulk_BPM(fIn);
        %sx = get_sx_AWR2243_bulk_BPM(fIn);
        torso=find_torso(sx);
%         leg_track=limb_tracker(sx);
        torso_vel_env=prf_idx2vel(torso);
       
%        limb_vel_env=prf_idx2vel(leg_track);
%         limb_vel_env=-1*limb_vel_env;

       % save(fOut,'torso_vel_env','torso')
        toc
end