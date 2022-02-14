 
function vel_env=prf_idx2vel(prf_env)
        % this script assume original sx variable has 4096 rows.
        % This script only works when you extract the trorso and upper envelope
        % from the spectrogram matrix, not from the saved .png/jpg
        new_env=prf_env;
y_vel=linspace(-6.2338,6.2338,4096);
 y_act=new_env+1000; % add lower crop point to move it to actual
for kk=1:length(new_env) 
        y_new(kk)=y_vel(round(y_act(kk)));
end
vel_env=y_new; 

% interpolation
% leg_vel=-1*y_new; % y axis in origianl sx matrix was going bottom to top, but in MS matrix it's going top to bottom.
% %so to take care of that we need to multiply with -1.
% 
% Inerpolate, if we want to plot radar data with vicon, then interpolate radar
% data to match with vicon

        %  
% 
% figure; plot(leg_vel,'r'); hold on; plot(AP_walk,'--g');
