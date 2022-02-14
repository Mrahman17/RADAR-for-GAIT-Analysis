function intrp_env=interp_to_match_vicon(env,total_vicon_samples)
     %  This funciton will increase the number of samples in radar envelope to
     %  macth it's length with the total VICON samples
     % Inputs:
        % env= Radar envelop
        % total_vicon_samples= How many samples are there in vicon data
        
        x=1:17995;%length(env);
        v=env;
        xq=linspace(1,total_vicon_samples,total_vicon_samples);
        intrp_env=interp1(x,v,xq);
end
