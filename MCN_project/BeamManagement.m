
function [ BeamManagementRSRQ , BeamManagementRSRP] = BeamManagement(x,y,z)
carrier = nrCarrierConfig;
% Maximum transmission bandwidth configuration for 50 MHz carrier with 30 kHz subcarrier spacing
carrier.NSizeGrid = 133;
carrier.SubcarrierSpacing = 30;
carrier.NSlot = 0;
carrier.NFrame = 1;

numNZPRes = 12;
csirs = nrCSIRSConfig;
csirs.CSIRSType = repmat({'nzp'},1,numNZPRes);
csirs.CSIRSPeriod = 'on';
csirs.Density = repmat({'one'},1,numNZPRes);
csirs.RowNumber = repmat(2,1,numNZPRes);
csirs.SymbolLocations = {0,1,2,3,4,5,6,7,8,9,10,11};
csirs.SubcarrierLocations = repmat({0},1,numNZPRes);
csirs.NumRB = 25;

% Validate CSI-RS antenna ports
validateCSIRSPorts(csirs);

% Get the binary vector to represent the presence of each CSI-RS resource
% in a specified slot
csirsTransmitted = getActiveCSIRSRes(carrier,csirs);

powerCSIRS = 0;

csirsSym = nrCSIRS(carrier,csirs,'OutputResourceFormat','cell');

csirsInd = nrCSIRSIndices(carrier,csirs,'OutputResourceFormat','cell');

% Set the carrier frequency
fc = 28e9;
freqRange = validateFc(fc);

% Set the propagation speed
c = physconst('LightSpeed');
% Calculate wavelength
lambda = c/fc;

txArySize = [8 8];
rxArySize = [2 2];

nTx = prod(txArySize);
nRx = prod(rxArySize);



% Configure antenna array positions
txArrayPos = [0;0;0];
rxArrayPos = [x;y;z];

% Calculate the free space path loss
toRxRange = rangeangle(txArrayPos,rxArrayPos);
spLoss = fspl(toRxRange,lambda);

% Initialize the flags to choose between URA and ULA
isTxRectArray = false;
isRxRectArray = false;

% Enable isTxRectArray if both the number of rows and columns of transmit
% antenna array are greater than one
if ~any(txArySize == 1)
    isTxRectArray = true;
end
% Enable isRxRectArray if both the number of rows and columns of receive
% antenna array are greater than one
if ~any(rxArySize == 1)
    isRxRectArray = true;
end

% Configure the transmit and receive antenna elements
txAntenna = phased.IsotropicAntennaElement('BackBaffled',true);  % To avoid transmission beyond +/- 90
                                                                 % degrees from the broadside, baffle
                                                                 % the back of the transmit antenna
                                                                 % element by setting the BackBaffled
                                                                 % property to true
rxAntenna = phased.IsotropicAntennaElement('BackBaffled',false); % To receive the signal from 360 degrees,
                                                                 % set the BackBaffled property to false

% Configure transmit antenna array
if isTxRectArray
    % Create a URA System object for signal transmission
    txArray = phased.URA('Element',txAntenna,'Size',txArySize,'ElementSpacing',lambda/2);
else
    % Create a ULA System object for signal transmission
    txArray = phased.ULA('Element',txAntenna,'NumElements',nTx,'ElementSpacing',lambda/2);
end

% Configure receive antenna array
if isRxRectArray
    % Create a URA System object for signal reception
    rxArray = phased.URA('Element',rxAntenna,'Size',rxArySize,'ElementSpacing',lambda/2);
else
    % Create a ULA System object for signal reception
    rxArray = phased.ULA('Element',rxAntenna,'NumElements',nRx,'ElementSpacing',lambda/2);
end

fixedScatMode = true;
rng(42);
if fixedScatMode
    % Fixed single scatterer location
    numScat = 1;
    scatPos = [60;10;15];
else
    % Generate scatterers at random positions
    numScat = 10; %#ok<UNRCH> 
    azRange = -180:180;
    randAzOrder = randperm(length(azRange));
    elRange = -90:90;
    randElOrder = randperm(length(elRange));
    azAngInSph = deg2rad(azRange(randAzOrder(1:numScat)));
    elAngInSph = deg2rad(elRange(randElOrder(1:numScat)));
    r = 20;
    
    % Transform spherical coordinates to Cartesian coordinates
    [x,y,z] = sph2cart(azAngInSph,elAngInSph,r);
    scatPos = [x;y;z] + (txArrayPos + rxArrayPos)/2;
end

txArrayStv = phased.SteeringVector('SensorArray',txArray,'PropagationSpeed',c);

[~,scatAng] = rangeangle(scatPos(:,1),txArrayPos); % Pointing towards the first scatterer position, in case of multiple scatterers

azTxBeamWidth = 60; % In degrees
elTxBeamWidth = 30; % In degrees

ssbTxAng = getInitialBeamDir(scatAng,azTxBeamWidth,elTxBeamWidth);

% Get the number of transmit beams based on the number of active CSI-RS resources in a slot
numBeams = sum(csirsTransmitted);

% Get the azimuthal sweep range based on the SSB transmit beam direction
% and its beamwidth in azimuth plane
azSweepRange = [ssbTxAng(1) - azTxBeamWidth/2 ssbTxAng(1) + azTxBeamWidth/2];

% Get the elevation sweep range based on the SSB transmit beam direction
% and its beamwidth in elevation plane
elSweepRange = [ssbTxAng(2) - elTxBeamWidth/2 ssbTxAng(2) + elTxBeamWidth/2];

% Get the azimuth and elevation angle pairs for all NZP-CSI-RS transmit beams
azBW = beamwidth(txArray,fc,'Cut','Azimuth');
elBW = beamwidth(txArray,fc,'Cut','Elevation');
csirsBeamAng = hGetBeamSweepAngles(numBeams,azSweepRange,elSweepRange,azBW,elBW);

wT = zeros(nTx,numBeams);
for beamIdx = 1:numBeams
    tempW = txArrayStv(fc,csirsBeamAng(:,beamIdx));
    wT(:,beamIdx) = tempW;
end

% Number of CSI-RS antenna ports
ports = csirs.NumCSIRSPorts(1);
% Initialize the beamformed grid
bfGrid = nrResourceGrid(carrier,nTx);
% Get the active NZP-CSI-RS resource indices
activeRes = find(logical(csirsTransmitted));
for resIdx = 1:numNZPRes
    % Initialize the carrier resource grid for one slot and map NZP-CSI-RS symbols onto
    % the grid
    txSlotGrid = nrResourceGrid(carrier,ports);
    txSlotGrid(csirsInd{resIdx}) = db2mag(powerCSIRS)*csirsSym{resIdx};
    reshapedSymb = reshape(txSlotGrid,[],ports);
    
    % Get the transmit beam index
    beamIdx = find(activeRes == resIdx);
    
    % Apply the digital beamforming
    if ~isempty(beamIdx)
        bfSymb = reshapedSymb * wT(:,beamIdx)';
        bfGrid = bfGrid + reshape(bfSymb,size(bfGrid));
    end
end

% Perform OFDM modulation
[tbfWaveform,ofdmInfo] = nrOFDMModulate(carrier,bfGrid);

% Normalize the beamformed time-domain waveform over the number of transmit
% antennas
tbfWaveform = tbfWaveform/sqrt(nTx);

chan = phased.ScatteringMIMOChannel;
chan.PropagationSpeed = c;
chan.CarrierFrequency = fc;
chan.Polarization = 'none';
chan.SpecifyAtmosphere = false;
chan.SampleRate = ofdmInfo.SampleRate;
chan.SimulateDirectPath = false;
chan.ChannelResponseOutputPort = true;

% Configure transmit array parameters
chan.TransmitArray = txArray;
chan.TransmitArrayMotionSource = 'property';
chan.TransmitArrayPosition = txArrayPos;

% Configure receive array parameters
chan.ReceiveArray = rxArray;
chan.ReceiveArrayMotionSource = 'property';
chan.ReceiveArrayPosition = rxArrayPos;

% Configure scatterers
chan.ScattererSpecificationSource = 'property';
chan.ScattererPosition = scatPos;
chan.ScattererCoefficient = ones(1,numScat);

% Get the maximum channel delay by transmitting random signal
[~,~,tau] = chan(complex(randn(chan.SampleRate*1e-3,nTx), ...
    randn(chan.SampleRate*1e-3,nTx)));
maxChDelay = ceil(max(tau)*chan.SampleRate);

% Append zeros to the transmit waveform to account for channel delay
tbfWaveform = [tbfWaveform; zeros(maxChDelay,nTx)];
% Pass the waveform through the channel
fadWave = chan(tbfWaveform);

% Configure the receive gain
rxGain = 10.^((spLoss)/20); % Gain in linear scale
% Apply the gain
fadWaveG = fadWave*rxGain;

% Configure the SNR in dB
SNRdB = 20;
SNR = 10^(SNRdB/10); % SNR in linear scale
% Calculate the standard deviation for AWGN
N0 = 1/sqrt(2.0*nRx*double(ofdmInfo.Nfft)*SNR);

% Generate AWGN
noise = N0*complex(randn(size(fadWaveG)),randn(size(fadWaveG)));
% Apply AWGN to the waveform
rxWaveform = fadWaveG + noise;

% Generate reference symbols and indices
refSym = nrCSIRS(carrier,csirs);
refInd = nrCSIRSIndices(carrier,csirs);

% Estimate timing offset
offset = nrTimingEstimate(carrier,rxWaveform,refInd,refSym);
if offset > maxChDelay
    offset = 0;
end

% Correct timing offset
syncTdWaveform = rxWaveform(1+offset:end,:);

rxGrid = nrOFDMDemodulate(carrier,syncTdWaveform);

rxArrayStv = phased.SteeringVector('SensorArray',rxArray,'PropagationSpeed',c);

[~,scatRxAng] = rangeangle(scatPos(:,1),rxArrayPos); % Pointing towards the first scatterer position, in case of multiple scatterers

azRxBeamWidth = 30; % In degrees
elRxBeamWidth = 30; % In degrees

rxAng = getInitialBeamDir(scatRxAng,azRxBeamWidth,elRxBeamWidth);

wR = rxArrayStv(fc,rxAng);

temp = rxGrid;
if strcmpi(freqRange,'FR1')
    % Beamforming without combining
    rbfGrid = reshape(reshape(temp,[],nRx).*wR',size(temp,1),size(temp,2),[]);
else % 'FR2'
    % Beamforming with combining
    rbfGrid = reshape(reshape(temp,[],nRx)*conj(wR),size(temp,1),size(temp,2),[]);
end

sceneParams.TxArray = txArray;
sceneParams.RxArray = rxArray;
sceneParams.TxArrayPos = txArrayPos;
sceneParams.RxArrayPos = rxArrayPos;
sceneParams.ScatterersPos = scatPos;
sceneParams.Lambda = lambda;
sceneParams.ArrayScaling = 100;   % To enlarge antenna arrays in the plot
sceneParams.MaxTxBeamLength = 45; % Maximum length of transmit beams in the plot
sceneParams.MaxRxBeamLength = 25; % Maximum length of receive beam in the plot

hPlotSpatialMIMOScene(sceneParams,wT,wR);
axis tight;
view([74 29]);

% Perform RSRP measurements
meas = hCSIRSMeasurements(carrier,csirs,rbfGrid);

% Display the measurement quantities for all CSI-RS resources in dBm
RSRPdBm = meas.RSRPdBm;
disp(['RSRP measurements of all CSI-RS resources (in dBm):' 13 num2str(RSRPdBm)]);

% Get the transmit beam index with maximum RSRP value
[~,maxRSRPIdx] = max(RSRPdBm(logical(csirsTransmitted)));

% Get the CSI-RS resource index with maximum RSRP value
[~,maxRSRPResIdx] = max(RSRPdBm);

% Get the steering weights corresponding to refined transmit beam
if numBeams == 0
    disp('Refinement has not happened, as NZP-CSI-RS is not transmitted')
else
    refBeamWts = wT(:,maxRSRPIdx);
    csirsAzBeamWidth = beamwidth(txArray,fc,'PropagationSpeed',c,'Weights',refBeamWts,'CutAngle',csirsBeamAng(2,maxRSRPIdx));
    csirsElBeamWidth = beamwidth(txArray,fc,'PropagationSpeed',c,'Weights',refBeamWts,'Cut','Elevation','CutAngle',csirsBeamAng(1,maxRSRPIdx));
    disp("Beam For Position "+x+","+y+","+z+" is "+maxRSRPIdx);
%     disp(['From initial beam acquisition:' 13 '    Beamwidth of initial SSB beam in azimuth plane is: '...
%         num2str(azTxBeamWidth) ' degrees' 13 ...
%         '    Beamwidth of initial SSB beam in elevation plane is: '...
%         num2str(elTxBeamWidth) ' degrees' 13 13 ...
%         'With transmit-end beam refinement:' 13 '    Refined transmit beam ('...
%         num2str(maxRSRPIdx) ') corresponds to CSI-RS resource '...
%         num2str(maxRSRPResIdx) ' is selected in the direction ['...
%         num2str(csirsBeamAng(1,maxRSRPIdx)) ';' num2str(csirsBeamAng(2,maxRSRPIdx))...
%         ']' 13 '    Beamwidth of refined transmit beam in azimuth plane is: '...
%         num2str(csirsAzBeamWidth) ' degrees' 13 ...
%         '    Beamwidth of refined transmit beam in elevation plane is: '...
%         num2str(csirsElBeamWidth) ' degrees']);
end

BeamManagementRSRQ = meas.RSRQ(maxRSRPResIdx);
BeamManagementRSRP = meas.RSRPdBm(maxRSRPResIdx);

end















function validateCSIRSPorts(csirs)
%   validateCSIRSPorts validates the CSI-RS antenna ports, given the
%   CSI-RS configuration object CSIRS.

    numPorts = csirs.NumCSIRSPorts;
    if any(numPorts > 1)
        error('nr5g:PortsGreaterThan1','CSI-RS resources must be configured for single-port for RSRP measurements.');
    end
end



function csirsTransmitted = getActiveCSIRSRes(carrier,csirs)
%   getActiveCSIRSRes returns a binary vector indicating the presence of
%   all CSI-RS resources in a specified slot, given the carrier
%   configuration object CARRIER and CSI-RS configuration object CSIRS.

    % Extract the following properties of carrier
    NSlotA        = carrier.NSlot;          % Absolute slot number
    NFrameA       = carrier.NFrame;         % Absolute frame number
    SlotsPerFrame = carrier.SlotsPerFrame;  % Number of slots per frame
    
    % Calculate the appropriate frame number (0...1023) based on the
    % absolute slot number
    NFrameR = mod(NFrameA + fix(NSlotA/SlotsPerFrame),1024);
    % Relative slot number (0...slotsPerFrame-1)
    NSlotR = mod(NSlotA,SlotsPerFrame);
    
    % Loop over the number of CSI-RS resources
    numCSIRSRes = numel(csirs.CSIRSType);
    csirsTransmitted = zeros(1,numCSIRSRes);
    csirs_struct = validateConfig(csirs);
    for resIdx = 1:numCSIRSRes
        % Extract the CSI-RS slot periodicity and offset
        if isnumeric(csirs_struct.CSIRSPeriod{resIdx})
            Tcsi_rs = csirs_struct.CSIRSPeriod{resIdx}(1);
            Toffset = csirs_struct.CSIRSPeriod{resIdx}(2);
        else
            if strcmpi(csirs_struct.CSIRSPeriod{resIdx},'on')
                Tcsi_rs = 1;
            else
                Tcsi_rs = 0;
            end
            Toffset = 0;
        end
        
        % Check for the presence of CSI-RS, based on slot periodicity and offset
        if (Tcsi_rs ~= 0) && (mod(SlotsPerFrame*NFrameR + NSlotR - Toffset, Tcsi_rs) == 0)
            csirsTransmitted(resIdx) = 1;
        end
    end
end

function freqRange = validateFc(fc)
%   validateFc validates the carrier frequency FC and returns the frequency
%   range as either 'FR1' or 'FR2'.

    if fc >= 410e6 && fc <= 7.125e9
        freqRange = 'FR1';
    elseif fc >= 24.25e9 && fc <= 52.6e9
        freqRange = 'FR2';
    else
        error('nr5g:invalidFreq',['Selected carrier frequency is outside '...
            'FR1 (410 MHz to 7.125 GHz) and FR2 (24.25 GHz to 52.6 GHz).']);
    end
end

function beamDir = getInitialBeamDir(scatAng,azBeamWidth,elBeamWidth)
%   getInitialBeamDir returns the initial beam direction BEAMDIR, given the
%   angle of scatterer position with respect to transmit or receive antenna
%   array SCATANG, beamwidth of transmit or receive beam in azimuth plane
%   AZBEAMWIDTH, and beamwidth of transmit or receive beam in elevation
%   plane ELBEAMWIDTH.

    % Azimuth angle boundaries of all transmit/receive beams
    azSSBSweep = -180:azBeamWidth:180;
    % Elevation angle boundaries of all transmit/receive beams
    elSSBSweep = -90:elBeamWidth:90;
    
    % Get the azimuth angle of transmit/receive beam
    azIdx1 = find(azSSBSweep <= scatAng(1),1,'last');
    azIdx2 = find(azSSBSweep >= scatAng(1),1,'first');
    azAng = (azSSBSweep(azIdx1) + azSSBSweep(azIdx2))/2;
    
    % Get the elevation angle of transmit/receive beam
    elIdx1 = find(elSSBSweep <= scatAng(2),1,'last');
    elIdx2 = find(elSSBSweep >= scatAng(2),1,'first');
    elAng = (elSSBSweep(elIdx1) + elSSBSweep(elIdx2))/2;
    
    % Form the azimuth and elevation angle pair (in the form of [az;el])
    % for transmit/receive beam
    beamDir = [azAng;elAng];
end

function beamAng = hGetBeamSweepAngles(numBeams,azSweepRange,elSweepRange,azBW,elBW,varargin)
%hGetBeamSweepAngles Return azimuth and elevation angle pairs for all beams
%
%   BEAMANG = hGetBeamSweepAngles(NUMBEAMS,AZSWEEPRANGE,ELSWEEPRANGE, ...
%   AZBW,ELBW) returns azimuth and elevation angle pairs BEAMANG, given the
%   number of beams NUMBEAMS, angular sweep range in azimuth plane
%   AZSWEEPRANGE, angular sweep range in elevation plane ELSWEEPRANGE,
%   beamwidth of antenna array in azimuth plane AZBW, and beamwidth of
%   antenna array in elevation plane ELBW.
%   The beams are equally spaced within the sweeps in both azimuth and
%   elevation directions.
%
%   BEAMANG = hGetBeamSweepAngles(...,ELSWEEP) also allows a boolean flag
%   ELSWEEP to enable or disable the elevation sweep.

%   Copyright 2020 The MathWorks, Inc.

    narginchk(5,6);
    if nargin == 5
        elSweep = true;
    else
        elSweep = varargin{1};
    end
    if ~elSweep
        elSweepRange = [0 0];
    end

    if numBeams == 0
        beamAng = [];
    else
        if numel(unique(azSweepRange)) == 2 && numel(unique(elSweepRange)) == 2
            % Both azimuth and elevation sweep
            effNumAzBeams = ceil((azSweepRange(2) - azSweepRange(1))/azBW);
            effNumElBeams = ceil((elSweepRange(2) - elSweepRange(1))/elBW);
            effAzElBeamsRatio = effNumAzBeams/effNumElBeams;

            temp = numBeams./(numBeams:-1:1);
            divVals = temp(temp == floor(temp))';
            pairsAzEl = [divVals flipud(divVals)];
            [~,index] = min(abs((pairsAzEl(:,1)./pairsAzEl(:,2)) - effAzElBeamsRatio));

            beamsInAz = pairsAzEl(index,1);
            if beamsInAz == 1
                azDir = (azSweepRange(1) + azSweepRange(2))/2;
            else
                azBeamSeparation = (azSweepRange(2) - azSweepRange(1))/beamsInAz;
                azDir = azSweepRange(1) + azBeamSeparation:azBeamSeparation:azSweepRange(2);
            end
            
            beamsInEl = pairsAzEl(index,2);
            if beamsInEl == 1
                elDir = (elSweepRange(1) + elSweepRange(2))/2;
            else
                elBeamSeparation = (elSweepRange(2) - elSweepRange(1))/beamsInEl;
                elDir = elSweepRange(1) + elBeamSeparation:elBeamSeparation:elSweepRange(2);
            end
            azDir = repmat(azDir,1,beamsInEl);
            elDir = reshape(repmat(elDir,beamsInAz,1),[],1)';
            
        elseif numel(unique(azSweepRange)) == 2 && numel(unique(elSweepRange)) == 1 
            % Azimuth sweep only
            beamSeparation = (azSweepRange(2) - azSweepRange(1))/numBeams;
            azDir = azSweepRange(1) + beamSeparation:beamSeparation:azSweepRange(2);
            elDir = elSweepRange(1)*ones(1,numBeams);

        elseif numel(unique(azSweepRange)) == 1 && numel(unique(elSweepRange)) == 2 
            % Elevation sweep only
            beamSeparation = (elSweepRange(2) - elSweepRange(1))/numBeams;
            azDir = azSweepRange(1)*ones(1,numBeams);
            elDir = elSweepRange(1) + beamSeparation:beamSeparation:elSweepRange(2);

        else
            % No sweep
            azDir = azSweepRange(1)*ones(1,numBeams);
            elDir = elSweepRange(1)*ones(1,numBeams);
        end
        
        % Form azimuth and elevation angle pairs for all beams
        beamAng = [azDir;elDir];
    end
end
