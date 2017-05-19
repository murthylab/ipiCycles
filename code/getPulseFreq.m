function [Fpeak, amp, F] = getPulseFreq(pulses, threshold, Fs)
% [Fpeak, amp, F] = getPulseFreq(pulses, threshold=1/exp(1), Fs=10000)

if ~exist('threshold','var') || isempty(threshold)
   threshold = 1/exp(1);
end
if ~exist('Fs','var') || isempty(Fs)
   Fs = 10000;
end

fftLen = 2000;
tmp = abs(fft(pulses, fftLen))';
amp = tmp(:,1:fftLen/2);
F = (1:ceil(fftLen/2)) * Fs/fftLen;
Fpeak = centerOfMass(F, amp, threshold)';
