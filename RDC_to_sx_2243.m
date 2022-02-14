
function sx=RDC_to_sx_2243(RDC)
  rp = fft(RDC(:,:,1));


    %% MTI v2
    [b,a]=butter(1, 0.01, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
%                                      [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
%                                      highpass, Wn is cutoff freq (half the sample rate)
    [m,n]=size(rp(:,:,1));
    rngpro=zeros(m,n);
    for k=1:size(rp,1)
        rngpro(k,:)=filter(b,a,rp(k,:,1));
    end
    %% STFT
    rBin = 60:160; %covid 18:30, front ingore= 7:nts/2, %lab 15:31 for front
    nfft = 2^12;window = 256;noverlap = 192;shift = window - noverlap;
%      sx = myspecgramnew(rngpro(rBin,:),window,nfft,shift);
     sx = myspecgramnew(sum(rngpro(rBin,:)),window,nfft,shift); % mti filter and IQ correction
%     sx2 = abs(flipud(fftshift(sx,1)));
%     %% Spectrogram
%     timeAxis = [1:NPpF*NoF]*SweepTime/NPpF ; % Time
%     freqAxis = linspace(-prf/2,prf/2,nfft); % Frequency Axis
%     fig = figure('visible','on');
%     colormap(jet(256));
%     set(gca,'units','normalized','outerposition',[0,0,1,1]);
%     doppSignMTI = imagesc(timeAxis,[-prf/2 prf/2],20*log10(abs(sx2/max(sx2(:)))));
% 
%     caxis([-50 0]) % 40
%     set(gca, 'YDir','normal')
%     set(gcf,'color','w');
% %     colorbar;
%     axis([0 timeAxis(end) -prf/4 prf/4])
% %     saveas(fig,[fOut(1:end-4) '.fig']);
%     set(gca,'xtick',[],'ytick',[])
%     frame = frame2im(getframe(gca));
%     imwrite(frame,[fOut(1:end-4) '.png']);
%     close all