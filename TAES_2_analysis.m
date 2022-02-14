clc; clear; close all; warning off
datapath = '/mnt/HDD02/WGAN/GAIT DATA/Original/train/02/';

saveColor ='/mnt/HDD02/WGAN/GAIT DATA/analysis/ENV_09/';


pattern = strcat(datapath, '*.png');    % file pattern
% pattern = strcat(datapath, '*.png');    % file pattern
files = dir(pattern);

I_MAX = numel(files); % # of files in "files" 

all_env=cell(I_MAX,2);
for ii =2:I_MAX 
    fIn = strcat(datapath,files(ii).name);
    [upper_env, central_env, lower_env] = env_find(fIn);
    all_env{ii,1}=pix2vel(upper_env);
    all_env{ii,2}=pix2vel(central_env);  
end
close all;
jj=55;
figure; plot(-1*all_env{jj,2}); grid on; grid minor;
torso=-1*all_env{jj,2};
torso_cut=torso(116:355);

yd=diff(torso_cut);
[b,a]=butter(1,2/28);
torso_filt=filter(b,a,yd);
figure; plot(torso_filt);grid on; grid minor;
figure; plot(smooth(torso_filt));grid on; grid minor;