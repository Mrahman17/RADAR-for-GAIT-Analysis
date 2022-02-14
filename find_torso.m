
function torso=find_torso(sx)
trimY_L=1000;
trimY_U=3000;
sx2 = abs(flipud(fftshift(sx,1)));
sx4=flipud(sx2(trimY_L:trimY_U,:));
img=20*log10(sx4./max(sx4(:)));



for t = 1:size(img,2)  
 
     vt =( max( img(:,t)));
    mt=find(img(:,t)==vt);
    torso(t) = mt(1);
end
%torso(torso<6) = max(torso);
%  torso2=median(torso,30);
figure; colormap(jet(256))
imagesc(img);caxis([-40 0]);
frame = frame2im(getframe(gca));
imwrite(frame,'neimage.png');
hold on; plot(torso,'m','LineWIdth',1.5);


