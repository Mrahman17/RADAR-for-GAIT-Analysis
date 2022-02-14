function vel = pix_to_vel(pix)

pix = squeeze(squeeze(pix));
crop_height = 343;
prf = 6400;
orig_im_height = 343;
crop_resize = 343;
limits = linspace(-prf/ 4, prf/4, prf+1);
freq_per_pix = length(limits)/orig_im_height;
idx = floor(pix /crop_resize*crop_height / crop_height*orig_im_height * freq_per_pix)+1;

idx(idx>length(limits)) = length(limits);

freq = limits(idx);
vel = -3*10^8*freq/2/(77*10^9); % lift off (-) for lower env


end
%((3e8)*(prf/4)/(2*77e9))=3.1169
