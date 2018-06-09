clc
clear
close all

%%% CREATE DMG OBJECT & ONE TRANSMIT FRAME

dmg = wlanDMGConfig;
dmg.MCS = 2; % Choose modulation technique
dmg.TrainingLength = 0; % No beam training
dmg.PSDULength = 100; % 100 bytes of data
psdu = randi([0 1],dmg.PSDULength*8,1); % Create communication bits
tx = wlanWaveformGenerator(psdu,dmg); % Create 1 transmit frame
K = length(tx); % number of symbols per frame
s = tx/rms(tx); % normalize transmit signal
Es = sum(abs(tx).^2)/K; % transmit signal energy per symbol
ind = wlanFieldIndices(dmg); % indices of 1 frame

%%% CREATE TX PLANAR ARRAY ON SOURCE VEHICLE

c = physconst('LightSpeed');
fc = 60.48e9;
lambda = c/fc;
d = lambda/2; % antenna array element spacing
txArray = phased.URA([2 8],'ElementSpacing',d);

%%% CALCULATE TX STEERING VECTOR & TX BEAMFORMING VECTOR

txSteerVec = phased.SteeringVector('SensorArray',txArray);
AoD = [90;0]; % AoD is preambleDataAngle - the angle at which you radiate the transmit signal
aTX = txSteerVec(fc,AoD); % Transmit steering vector (phase shifts experienced by a plane wave incident at AoD) - does not change within CPI 
fTX = conj(aTX); % Transmit beamforming weights (phase shifts used to make the incoming plane wave add up constructively)

%%% CREATE RX PLANAR ARRAY ON RECIPIENT VEHICLE

rxArray = txArray;
rxSteerVec = txSteerVec; % since we use the same array for TX and RX

%%% CALCULATE RX STEERING VECTOR

AoA = AoD; % since we only consider LOS & assume perfect beam alignment
aRXcom = rxSteerVec(fc,AoA); % Receive steering vector (for communication, at target vehicle)

%%% CALCULATE RX BEAMFORMING VECTOR

beamformer = phased.PhaseShiftBeamformer('SensorArray',rxArray,'OperatingFrequency',fc,'PropagationSpeed',c,'DirectionSource','Property','Direction',AoA);
beamformer.WeightsOutputPort = true;
txAnt = phased.Transmitter('PeakPower',0.1,'Gain',0); 
tx = txAnt(tx); % Transmit signal with increased power and antenna gain
rxElementSig = collectPlaneWave(rxArray,tx,AoA,fc,c); % Contains received signal at each element in rxArray
NF = 6;
rxAnt = phased.ReceiverPreamp('Gain',20,'LossFactor',0','NoiseMethod','Noise temperature','SampleRate',(1.76e9)/2,'NoiseFigure',NF); % 20, 41.5818
rxElementSig1 = rxAnt(rxElementSig); % Receive signal with antenna gain and loss due to NF
[rxSig,fRXcom] = beamformer(rxElementSig1); % Apply beamforming gain to received signal
fRXrad = conj(fRXcom); % Receive beamforming weights (for radar, at source vehicle)

%%% SET CONSTANT PARAMETERS

m = 1; % frame index
W = 1.76e9; % signaling bandwidth
Ts = 1/W; % symbol period = 1/W
T = 4.2e-3; % CPI duration
M = floor(T/(K*Ts)); % no. of frames in one CPI
H_nlos = randn(16)+1i*randn(16); % NLOS modelled as G.R.V
beta0 = exp(1i*3000000);
v0 = 224.94; % relative velocity of receipient vehicle
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

% Calculate SNR and SCNR for different distances

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

% Plot results

figure(1);
plot(distance,SNRdb1,'b-.')
axis([0 200 -60 80])
grid ON;
hold
plot(distance,SNRdb2,'b-')
plot(distance,SCNRdb1,'r-.')
plot(distance,SCNRdb2,'r-')
legend('Comm PL 2.0','Comm PL 2.5','Rad PL 2.0','Rad PL 2.5')
hold

%%% SIMULATE RECEIVED SIGNAL AT SOURCE VEHICLE

% Respecify channel gain

p0 = 100.21;
G0_1 = ((lambda/(4*pi*p0^2))^PL1)*(RCS0/(4*pi));
h0_1 = sqrt(G0_1)*beta0*fRXrad'*conj(aRXcom)*aTX'*fTX;

% Assume that RX and TX pulse shaping done perfectly. 
% Rx signal is the noisy, channel-attenuated, doppler shifted version of TX signal.
% Consider only 1 frame.
% Frequency offset will cause phase noise to the baseband received tx. So you can include it as exp(-1i*2*pi*freqOffset*k*Ts) if needed.

rx = zeros(length(tx),1);
for k = 1:length(tx)
    rx(k,1) = sqrt(Es)*h0_1*exp(1i*2*pi*doppler_rad*k*Ts)*rxSig(k,1);
end
rx = [zeros(round(((2*p0)/c)/Ts),1); rx]; % artificially include delay

% FRAME SYNCHRONIZATION USING STF

txPreamble = tx(ind.DMGSTF(1,1):ind.DMGCE(1,2),1); % extract TX preamble
Ktr = length(txPreamble); % length of training sequence/preamble
rxPreamble = rx(1:Ktr,1);

txSTF = txPreamble(ind.DMGSTF(1,1):ind.DMGSTF(1,2));
Kstf = length(txSTF);
rxSTF = rxPreamble(1:Kstf,1);

P = 128; % window length
ND = 128; % distance between points selected at each iteration
R1 = zeros(2,1);
R12 = zeros(30*P,1);
R11 = zeros(30*P,1);

for l = 1:30*P
        R13 = 0;
    for nn = 0:P-1
        n = nn+1;
        if (l-n)>0 && (l-n-ND)>0
            R11(l,1) = R11(l,1) + rx(l-n)*conj(rx(l-n-ND));
        end
        if (l-n)>0
            R12(l,1) = R12(l,1) + (abs(rx(l-n)))^2;
        end
        if (l-n-ND)>0
            R13 = R13 + (abs(conj(rx(l-n-ND))))^2;
        end
    end
    R1(l,1) = R11(l,1)./(sqrt(R12(l,1))*sqrt(R13));
end

getOut = 0;
threshold = 0.25;

for k = 1:length(R1)
    if abs(R1(k,1))>threshold
       frameStart = k;
       count = count + 1;
    if count == 256
           getOut = 1;
    end
    else
        count = 0;
    end
    if getOut == 1
        frameStart = k-392;
        if frameStart < 0
            frameStart = 1;
        end
        break;
    end
end

l01 = frameStart;
figure(2)
plot(1:length(R1),abs(R1))
hold
plot(1:length(R1),threshold*ones(length(R1),1),'r-')
plotFrameStart = zeros(length(R1),1);
plotFrameStart(frameStart,1) = 1;
plot(1:length(R1),plotFrameStart,'g-','LineWidth',2);
hold

%%% ESTIMATE RANGE OF TARGET

R2 = zeros(4353,1);
for l = 1:28*P
    for nn = 0:Ktr-1
        n = nn+1;
        if n>0 && (l+n)>0
            R2(l,1) = R2(l,1) + txPreamble(n)*conj(rx(l+n));
        end
    end
    R2(l,1) = (1/Ktr)*R2(l,1);
end

[maximum,integerSymDelay] = max((abs(R2)).^2); % delay (in number of symbols)
l04 = integerSymDelay;
p0_est = integerSymDelay*Ts*0.5*c; 

%%% PROCESS IMAG PART OF R11 TO GET DOPPLER ESTIMATE (& FREQUENCY OFFSET ESTIMATE)

doppler_rad_est = angle(sum(R11(140:1500,1)))/(2*pi*128*Ts);
v0_est_coarse = floor(0.5*lambda*doppler_rad_est);

% figure
% plot(1:length(R11),angle(R11)) 
% Only the angle of R11 over STF (1:2176) is of concern. The slight DC
% offset from 0 here represents the phase shift due to frequency offset and 
% doppler shift due to velocity of target. If noise power is high, you'll
% also see a variation about the DC offset value.
% If there is residual frequency offset that cannot be removed, this estimate will not give accurate velocity estimate. Hence the need for
% Morelli/Mengali algorithms.
% Even without frequency offset, this does not achieve required 0.1m/s accuracy. But Morellai/Mengalli algorithms can achieve it.

%%% FINE DELAY ESTIMATE & CHANNEL GAIN ESTIMATE

havg = zeros(Ktr+round(((2*p0)/c)/Ts),1);
txCEF = tx(ind.DMGCE(1,1):ind.DMGCE(1,2),1);
a256u = conj(txCEF(1:256,1));
b256u = conj(txCEF(257:512,1));
a256v = conj(txCEF(513:768,1));
b256v = conj(txCEF(769:1024,1));

for l = 1:length(havg)
    for nn = 0:255
        n = nn+1;
        if (l+n)>0 && (n+l+256)>0
            havg(l,1) = havg(l,1) + rx(l+n)*a256u(n)+ rx(l+n+256)*b256u(n)+rx(l+n+512)*a256v(n) + rx(l+n+256+512)*b256v(n);
        end
    end
    havg(l,1) = (1/(2*256)/2)*havg(l,1);
end

havgmod(:,1) = (abs(havg)).^2;
[peak,fineDelayEst] = max(havgmod);
h0_1_est = havg(fineDelayEst); % channel estimate will be wrong since Doppler shift is left uncorrected
l03 = fineDelayEst-Kstf; % but if doppler shift is corrected, fineDelayEst cannot be estimated
chan_est = havg(fineDelayEst-128:fineDelayEst+127);

fprintf('\nActual target range: %.2f',p0)
fprintf('\nEstimated target range: %.2f',p0_est)
fprintf('\nActual target velocity: %.2f',v0)
fprintf('\nEstimated target velocity: %.2f',v0_est_coarse)
