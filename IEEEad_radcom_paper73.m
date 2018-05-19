clc
clear all
close all

%% Create dmg object and 1 frame

dmg = wlanDMGConfig;
dmg.MCS = 24;
dmg.TrainingLength = 0; 
dmg.PSDULength = 100; % 100 bytes of data
% seed = rng(0);
psdu = randi([0 1],dmg.PSDULength*8,1);
tx = wlanWaveformGenerator(psdu,dmg); % transmit signal
K = length(tx); % number of symbols per frame
s = tx/rms(tx); % normalized transmit signal
Es = sum(abs(tx).^2)/K; % transmit signal energy per symbol
% ind = wlanFieldIndices(dmg);

%% Create TX planar array on SOURCE VEHICLE

c = physconst('LightSpeed');
fc = 60.48e9;
lambda = c/fc;
d = lambda/2;
txArray = phased.URA([2 8],'ElementSpacing',d);
% viewArray(txArray)

%% Calculate TX steering vector and TX beamforming vector

txSteerVec = phased.SteeringVector('SensorArray',txArray);
AoD = [90;0]; % AoD is preambleDataAngle
aTX = txSteerVec(fc,AoD); % aTX do not change within CPI
fTX = conj(aTX); % preambleDataAMV aka TX Beamforming weights
% plotResponse(txArray,fc,c,'Format','Polar','RespCut','3D','Weights',fTX);

%% Create RX planar array on RECIPIENT VEHICLE

rxArray = txArray;
rxSteerVec = txSteerVec; % since we use the same array for TX and RX

%% Calculate RX steering vector

AoA = AoD; % since we only consider LOS
aRXcom = rxSteerVec(fc,AoA);

%% Calculate RX beamforming vector 

beamformer = phased.PhaseShiftBeamformer('SensorArray',rxArray,'OperatingFrequency',fc,'PropagationSpeed',c,'DirectionSource','Property','Direction',AoA);
beamformer.WeightsOutputPort = true;
rxElementSig = collectPlaneWave(rxArray,tx,AoA,fc,c); % Contains received signal at each element in rxArray
[rxSig,fRXcom] = beamformer(rxElementSig); 
fRXrad = conj(fRXcom);

%% Set constant parameters

m = 1; % frame index
W = 1.76e9; % signaling bandwidth
Ts = 1/W; % symbol period = 1/W
T = 4.2e-3; % CPI duration
M = floor(T/(K*Ts)); % no. of frames in one CPI
H_nlos = randn(16)+1i*randn(16); % NLOS modelled as G.R.V
beta0 = exp(1i*9348);
v0 = 20; % relative velocity of receipient vehicle
doppler_rad = (2*v0)/lambda; % round trip doppler shift
RCS0 = 10^(20/10);
PL1 = 2.0; % path loss exponent
PL2 = 2.5;
doppler_com = v0/lambda; % doppler shift
DopplerPhaseShift = exp(1i*2*pi*doppler_com*m*K*Ts);
alpha = 1; % does not affect SNR
Jcom = 9.74; % Rician K factor ranges from 6.35 to 15.1 dB. Mean 9.89dB.
boltz = physconst('Boltzmann');
temperature = 290; % kelvins
NF = 6; % in dB
noise_sigma = boltz*temperature*(W/2)*10^(NF/10); % noise power = kTBN
distance = [1:200]';

%% Calculate SNR and SCNR for different distances

for k = 1:length(distance)
    % Set varying parameters
    p0 = distance(k,1);
    
    % Calculate Effective Radar Channel Coeff
    tau0 = (2*p0)/c;
    G0_1 = ((lambda/(4*pi*p0^2))^PL1)*(RCS0/(4*pi));
    G0_2 = ((lambda/(4*pi*p0^2))^PL2)*(RCS0/(4*pi));
    h0_1 = sqrt(G0_1)*beta0*fRXrad'*conj(aRXcom)*aTX'*fTX;
    h0_2 = sqrt(G0_2)*beta0*fRXrad'*conj(aRXcom)*aTX'*fTX;

    % Calculate Effective Comm Channel Coeff
    Gcom1 = (lambda/(4*pi*p0))^PL1; % large scale channel gain
    Gcom2 = (lambda/(4*pi*p0))^PL2;
    H_los = alpha*DopplerPhaseShift*aRXcom*aTX';
    Hcom = (sqrt(Jcom/(Jcom+1)))*H_los + (sqrt(1/(Jcom+1)))*H_nlos;
    Heff1 = (sqrt(Gcom1))*fRXcom'*Hcom*fTX; % effective channel gain for one frame
    Heff2 = (sqrt(Gcom2))*fRXcom'*Hcom*fTX;
    
    %%% Calculate SNR and SCNR of received signal
    SNRdb1(k,1) = 10*log10((Es*(abs(Heff1))^2)/noise_sigma);
    SNRdb2(k,1) = 10*log10((Es*(abs(Heff2))^2)/noise_sigma);
    SCNRdb1(k,1) = 10*log10((Es*(abs(h0_1))^2)/(noise_sigma));
    SCNRdb2(k,1) = 10*log10((Es*(abs(h0_2))^2)/(noise_sigma));
end

%% Plot results

figure;
plot(distance,SNRdb1,'b-.')
axis([0 200 -60 80])
grid ON;
hold
plot(distance,SNRdb2,'b-')
plot(distance,SCNRdb1,'r-.')
plot(distance,SCNRdb2,'r-')
legend('Comm PL 2.0','Comm PL 2.5','Rad PL 2.0','Rad PL 2.5')
hold


