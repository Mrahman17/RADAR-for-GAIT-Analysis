%% NOT a working script

fname='sean_trial_02_Raw_0.bin';
fid = fopen(fname,'r');
    % DCA1000 should read in two's complement data
    Data = fread(fid, 'int16');
    fclose(fid);
    %% Parameters
    fileSize = size(Data, 1);

    numADCBits = 16; % number of ADC bits per sample
    SweepTime = 40e-3; % Time for 1 frame
    NTS = 256; %256 Number of time samples per sweep
    numADCSamples = NTS;
    numTX = 2; % '1' for 1 TX, '2' for BPM
    NoC = 128;%128; % Number of chirp loops
    NPpF = numTX*NoC; % Number of pulses per frame
    numRX = 4;
    sampleFreq = 6.25e6; % 2e6 ADC Sampling frequency
    slope = 66.578e12; %29.982e12; % Mhz / us = e6/e-6 = e12
    
    fstart = 77e9; % Start Frequency
     fstop = fstart+4e9;%1.79892e9;%   Stop Frequency
    fc = (fstart+fstop)/2; % Center Frequency
    c = physconst('LightSpeed'); % Speed of light
     lambda = c/fc; % Lambda
    numLanes = 4; % do not change. number of lanes is always 4 even if only 1 lane is used. unused lanes
    % NoF = fileSize/2/NPpF/numRX/NTS; % Number of frames
    numChirps = ceil(fileSize/2/NTS/numRX);
    NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
    Np = numChirps;%floor(size(Data(:,1),1)/NTS); % #of pulses
    dT = SweepTime/NPpF; % 
    prf = 1/dT; %
    
    Bw =4e9; % Bandwidth
    
    timeAxis = linspace(0,SweepTime*NoF,numChirps);%[1:NPpF*NoF]*SweepTime/NPpF ; % Time
    duration = max(timeAxis);
    
    idletime = 100e-6;
    adcStartTime = 6e-6;
    rampEndTime = 60e-6;
    
    isReal = 0; % set to 1 if real only data, 0 if complex dataare populated with 0 %% read file and convert to signed number

    %% Data reshape
    % if 12 or 14 bits ADC per sample compensate for sign extension
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        Data(Data > l_max) = Data(Data > l_max) - 2^numADCBits;
    end

    
%% organize data by LVDS lane
% for real only data
if isReal
% reshape data based on one samples per LVDS lane
Data = reshape(Data, numLanes, []);
%for complex data
else
% reshape and combine real and imaginary parts of complex number
LVDS = zeros(1, fileSize/2);

Data = reshape(Data, numLanes*2, []);
Data = Data([1,2,3,4],:) + sqrt(-1)*Data([5,6,7,8],:);
Data1=Data(1,:);
Data2=Data(2,:);
Data3=Data(3,:);
Data4=Data(4,:);
LVDS(1:4:end) = Data1(1:1:end);
LVDS(2:4:end) = Data2(1:1:end);
LVDS(3:4:end) = Data3(1:1:end);
LVDS(4:4:end) = Data4(1:1:end);
end
% check array size (if any frames were dropped, pad zeros)
if length(LVDS) ~= numADCSamples*numRX*numChirps
        numpad =  numADCSamples*numRX*numChirps - length(LVDS); % num of zeros to be padded
        LVDS = padarray(LVDS, [0 numpad],'post');
end

clear Data Data1 Data2 Data3 Data4
  LVDS = reshape(LVDS, numADCSamples*numRX, double(numChirps));
    %% If BPM, i.e. numTX = 2, see MIMO Radar sec. 4.2
    BPMidx = [1:2:numChirps-1];
    if numTX == 2
       LVDS_TX0 = 1/2 * (LVDS(:,BPMidx)+LVDS(:,BPMidx+1));
       LVDS_TX1 = 1/2 * (LVDS(:,BPMidx)-LVDS(:,BPMidx+1));
       LVDS0 = kron(LVDS_TX0,ones(1,2));
       LVDS1 = kron(LVDS_TX1,ones(1,2));
       LVDS = zeros(NTS*numRX*numTX,numChirps);
       LVDS(1:end/2,:) = LVDS0;
       LVDS(end/2+1:end,:) = LVDS1;
    end
  
    clear LVDS0 LVDS1 LVDSTX0 LVDSTX1 
    %% Organize data per RX
    DataN = zeros(numChirps*numADCSamples, numRX*numTX);
    for i = 1:numRX*numTX
        DataN(:,i) = reshape(LVDS((i-1)*numADCSamples+1:i*numADCSamples,:),[],1);
    end

     %% No IQ Correction
    rawData = reshape(DataN,NTS,numChirps, numRX*numTX); 
    
        %% Range-Velocity Map
        
        Rmax = sampleFreq*c/(2*slope);
        Tc = idletime+adcStartTime+rampEndTime;
        Tf = SweepTime;
        velmax = lambda/(Tc*4); % Unambiguous max velocity
        DFmax = velmax/(c/fc/2);
        rResol = c/(2*Bw);
        vResol = lambda/(2*Tf);
        % define frame size
        PN = NTS; %10 equally time spaced matricex: 10X500=5000
        RANGE_FFT_SIZE = NTS;
        DOPPLER_FFT_SIZE = PN*2; %*2
        
        
        RNGD2_GRID = linspace(0, Rmax, RANGE_FFT_SIZE);
        DOPP_GRID = linspace(DFmax, -DFmax, DOPPLER_FFT_SIZE);
        
        V_GRID = (c/fc/2)*DOPP_GRID;
        
        RCData = rawData(:,:,1);
        fps = 25;%1/SweepTime;
        n_frames = duration*fps;
        shft = floor(size(RCData,2)/n_frames);
        
      %%  CA-CFAR params
        numGuard = 4;
        numTrain = numGuard*2;
        P_fa = 1e-5; % Prob of false alarm
        SNR_OFFSET = -5; % -10
        %     cfar_bins = ones(2,n_frames);
        figure('Visible','off')%,
        % set(gcf,  'units', 'normalized','position', [0.2 0.2 0.4 0.6])
        
        for k = 1:n_frames
                
                RData_frame = RCData(:, 1+(k-1)*shft:k*shft);
                RData_frame = bsxfun(@minus, RData_frame, mean(RData_frame,2));   % subtract stationary objects
                G_frame = fftshift(fft2(RData_frame, RANGE_FFT_SIZE,DOPPLER_FFT_SIZE),2); % 2 adjust  shift
                RDM_dB = 10*log10(abs(G_frame)./max(abs(G_frame(:))));
                %         time_Counter = (k/n_frames)*duration;
                [RDM_mask, cfar_ranges, cfar_dopps] = ca_cfar_vid(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                center = centroid_cfar(RDM_mask,'centroid_RDC') ;
                
                
                    if center(1,k) ~= 0
                            hold off
                            imshow(RDM_mask(:,:,i))
                            hold on
                            scatter(center(2,i),center(1,i),'o','filled','MarkerFaceColor','red','linewidth',1)
                            drawnow;
                            F(k) = getframe(gca);
                    else
                            F(k) = im2frame(repmat(RDM_mask(:,:,i),1,1,3));
                    end
        end
                
        
        figure('visible','off')


writerObj = VideoWriter(fout);
writerObj.FrameRate = 25;
open(writerObj);

for i=1:length(F)
%         frame = 256*uint8(F(:,:,i));
        frame = F(i);
        writeVideo(writerObj, frame);
end
close(writerObj);
        

        
        
  %{              
                imagesc(V_GRID,RNGD2_GRID,RDM_dB);
                %         xlabel('Radial Velocity (m/s)','FontSize',13, 'FontName','Times')
                %         ylabel('Range (meter)','FontSize',13, 'FontName','Times')
                %         title({'Range-Velocity Map';num2str(time_Counter,'%.2f')},'FontSize',13, 'FontName','Times')
                %         colorbar
                set(gca, 'CLim',[-13,0]); % [-35,0],
                colormap(jet) % jet
                %         caxis([90 130]) % 90 130
                %         axis xy;
                axis([-velmax/numTX velmax/numTX 0 4])
                %         set(gcf, 'Position',  [100, 100, size(G_frame,1), size(G_frame,2)])
                
                drawnow
                F(k) = getframe(gca); % gcf returns the current figure handle
                
                %       colormap(gray)
                %       F2(k) =  getframe(gca);
                
        end
        
        %     fGray = [fNameOut(1:end-4) '_gray.avi'];
        fNameOut='RDVideo';
        writerObj = VideoWriter(fNameOut);
        writerObj.FrameRate = fps;
        open(writerObj);
        
        %     writerObj2 = VideoWriter(fGray);
        %     writerObj2.FrameRate = fps;
        %     open(writerObj2);
        
        for i=1:length(F)
                % convert the image to a frame
                frame = F(i) ;
                writeVideo(writerObj, frame);
                %         frame2 = F2(i);
                %         writeVideo(writerObj2, frame2);
        end
        close(writerObj);
        %     close(writerObj2);
        close all
                %}