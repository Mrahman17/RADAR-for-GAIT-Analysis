function rawData = RDC_extract_2243(fname)
    % read .bin file
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
    
    numLanes = 4; % do not change. number of lanes is always 4 even if only 1 lane is used. unused lanes
    % NoF = fileSize/2/NPpF/numRX/NTS; % Number of frames
    numChirps = ceil(fileSize/2/NTS/numRX);
    NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
    Np = numChirps;%floor(size(Data(:,1),1)/NTS); % #of pulses
    dT = SweepTime/NPpF; % 
    prf = 1/dT; %
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
tmp=0;
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
        if size(LVDS,2) > size(LVDS0,2)
               LVDS(:,end) = [];
               tmp = 1; % frame drop
               DataN = zeros((numChirps-1)*numADCSamples, numRX*numTX);
       else
               DataN = zeros(numChirps*numADCSamples, numRX*numTX);
       end
      
       LVDS(1:end/2,:) = LVDS0;
       LVDS(end/2+1:end,:) = LVDS1;
    end
  
    clear LVDS0 LVDS1 LVDSTX0 LVDSTX1 
    %% Organize data per RX
    %DataN = zeros(numChirps*numADCSamples, numRX*numTX);
    for i = 1:numRX*numTX
        DataN(:,i) = reshape(LVDS((i-1)*numADCSamples+1:i*numADCSamples,:),[],1);
    end
    if tmp
             rawData = reshape(DataN,NTS,numChirps-1, numRX*numTX); 
    else
            rawData = reshape(DataN,NTS,numChirps, numRX*numTX); 
    end
     %% No IQ Correction
    %rawData = reshape(DataN,NTS,numChirps, numRX*numTX); 
    rp = fft(rawData(:,:,1));
    clear DataN