% estEDR
% Estimate frequency location of posterior dominant activity from "alpha prominence"
% IZ, 26.11.18

% Analyse relative power in EEG channel. Choose intervals with rel.Alpha >
% 50% and not extreme global power.
% On mean spectrum from these intervals: Find peak.
% Input: channel-vector and EEG.srate. Output: PDR estimate. Peak
% amplitude. Seconds of data used.
% Band definitions fixed: bands variable
% Uses wavelet transform. S parameter linearly changing from 4 to 12.
% Frequency resolution 0.1 Hz.
% if more than one peak - take the tallest one

function [PDR,PDRpeak,secondsUsed] = estPDR(data,sr)

if nargin<2
   disp('Insufficient inputs')
   disp('Channel data and sampling rate')
   return
end

% band definitions
bands = [0 4; 4.1 7.9; 8 12.9; 13 20]; % delta, theta, alpha, beta

% Morlet wavelet transform
EEGtime = 0+1/sr:(1/sr):(length(data)/sr);
wavtime = -2:1/sr:2;
nData = length(EEGtime);
nKern = length(wavtime);
nConv = nData + nKern -1;
halfwav = floor(length(wavtime)/2)+1;

frex = 0.1:.1:20; % frequencies
nFrex = length(frex); % number of frequencies
s = linspace(4,12,nFrex) ./ (2*pi .* frex); % cycles between 4-12

% init
tf = zeros(nFrex, nData); % frex-by-samples-by-numberOfChannels

sigX = fft(data,nConv);

for ii = 1:nFrex
    cmw = exp(2*1i*pi*frex(ii)*wavtime - (wavtime .^2)/(2*s(ii)));
    cmwX = fft(cmw,nConv);
    cmwX = cmwX ./ max(cmwX);
    as = ifft(sigX .* cmwX);
    tf(ii,:) = abs(as(halfwav:end-halfwav+1)).^2; % power
end

% bands in frequency vector
deltaFreq = dsearchn(frex',bands(1,1)) : dsearchn(frex',bands(1,2));
thetaFreq = dsearchn(frex',bands(2,1)) : dsearchn(frex',bands(2,2));
alphaFreq = dsearchn(frex',bands(3,1)) : dsearchn(frex',bands(3,2));
betaFreq = dsearchn(frex',bands(4,1)) : dsearchn(frex',bands(4,2));

% spectral power ratio.
relPwr = zeros(4,nData); % rows 1-4, Delta, Theta, Alpha, Beta
globalPwr = zeros(1,nData); % init
for ii = 1:length(globalPwr)
D = trapz(tf(deltaFreq,ii));
T = trapz(tf(thetaFreq,ii));
A = trapz(tf(alphaFreq,ii));
B = trapz(tf(betaFreq,ii));
Total = D+T+A+B;
relPwr(1,ii) = round(D/Total*100,2);
relPwr(2,ii) = round(T/Total*100,2);
relPwr(3,ii) = round(A/Total*100,2);
relPwr(4,ii) = round(B/Total*100,2);
globalPwr(ii) = Total;
end

% sd
sdGPWR = std(globalPwr); % standard deviation of global pwr

% Relative Alpha Requirement
relA_Threshold = 50; % in percent

% alpha above 50 % and global pwr not extreme
[~, req] = find(globalPwr < sdGPWR*3 & relPwr(3,:)>relA_Threshold); % requirements fulfilled test
if length(req) < 60 % insufficient part of data fulfilling requirements
    PDR = NaN; PDRpeak = NaN; secondsUsed = NaN;
    return
end

% calculate "alpha prominence spectrum"
aSpec = mean(tf(:,req),2); % spectrum from selected parts of TF
%Spec = mean(tf,2); % mean spectrum

% find peak
[pk, loc] = findpeaks(aSpec,'SortStr','descend','MinPeakProminence',1.5);
if isempty(loc) % if prominence req too severe
    [pk, loc] = findpeaks(aSpec,'SortStr','descend');
end
PDR = frex(loc(1));
PDRpeak = pk(1);
secondsUsed = length(req)/sr;
disp(['PDR peak location: ',num2str(round(PDR,1)),' Hz.'])
end
