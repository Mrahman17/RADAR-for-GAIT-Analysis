function leg_track=limb_tracker(sx)
trimY_L=1000;
trimY_U=3000;
sx2 = abs(flipud(fftshift(sx,1)));
sx4=flipud(sx2(trimY_L:trimY_U,:));
img2=20*log10(sx4./max(sx4(:)));
img2(img2<-50)=-100; % This is equivalent to apply caxis([-34 0]), we assign them to a minimal value of -100
G=img2;

% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));

MS=Gs;
MS(MS<=1)=0;
% figure;
% imagesc(MS)  


% Find upper, lower and central envelope, which also contains buterworth filter
%  to smooth the curve a bit. otherwise it's really messy.


[upper_env, central_env, lower_env] = env_find(MS);
%figure;  plot(lower_env); hold on; plot(upper_env);hold on; plot(central_env)

% Combine upper and lower envelope, look at the index number of 0 hertz , (in
% this case it is 1000), if central and upper env is above 0 hertz line, then,
% new env should come from upper env, if central and upper env is below 0 hertz
% line, then new env should come from lower env.

a=upper_env;
b= central_env;
c=lower_env;
new_env=zeros(1,length(a));
for idx=1:length(central_env)
        if b(idx)<1000 
                if  a(idx)< 1000
                        
                        new_env(idx)=a(idx);
                else
                        new_env(idx)=b(idx);
                end
        end
        if b(idx)>1000
                if c(idx)> 1000
                        new_env(idx)=c(idx);
                else
                        new_env(idx)=b(idx);
                end
        end
       
end

  %figure;  plot(lower_env); hold on; plot(upper_env);hold on; plot(central_env);  hold on; plot(new_env,'r', 'LineWidth',2);               
 figure; imagesc(MS); hold on; plot(new_env,'r', 'LineWidth',1.3)
 leg_track=new_env;