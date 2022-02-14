function center = centroid_cfar(RDM_mask,fout) % finds centroid for each frame

center = zeros(2,size(RDM_mask,3));

for i = 1:size(RDM_mask,3)
   [row,col,v] = find(RDM_mask(:,:,i));
   if ~isempty(row)
       center(1,i) = round(sum(row.*v)/length(v));
       center(2,i) = round(sum(col.*v)/length(v));
   else
       continue
   end
end

figure('visible','off')
for i = 1:size(RDM_mask,3)
    if center(1,i) ~= 0
        hold off
        imshow(RDM_mask(:,:,i))
        hold on
        scatter(center(2,i),center(1,i),'o','filled','MarkerFaceColor','red','linewidth',1)
        drawnow;
        F(i) = getframe(gca);
    else
        F(i) = im2frame(repmat(RDM_mask(:,:,i),1,1,3));
    end
end

writerObj = VideoWriter(fout);
writerObj.FrameRate = 25;
open(writerObj);

for i=1:length(F)
%         frame = 256*uint8(F(:,:,i));
        frame = F(i);
        writeVideo(writerObj, frame);
end
close(writerObj);
        
end