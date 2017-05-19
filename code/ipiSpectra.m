function [spec, peak, a] = ipiSpectra(t,ipi, OFAC, HIFAC, alpha, F)
% Computes spectra, detects spectral peaks and generates interpolated
% spectrum for IPIs.
% [spec, peak, a] = ipiSpectra(t,ipi, OFAC, HIFAC, alpha, F)
%
% PARAMETERS
%  t     - vector of pulse times in seconds
%  ipi   - vector of pulse intervals in seconds
%  OFAC  - see lomb.m
%  HIFAC - see lomb.m
%  alpha - significance level for spectral peaks 
%  F     - vector of frequencies for interpolation
%
% RETURNS
%  spec  - F - frequencies, P - spectral power, p - probability
%          (see lomb.m for details)
%  peaks - spectral peaks using findpeaks 
%          (amp - power, loc - frequency, prob - probability, significant - frequency with significant peaks;
%  a     - a.F - frequency, a.spec - power (interpolate spectra so they are all on a common x-axis
%          across flies)

% lomb scargle periodogram
[spec.F, spec.P, spec.p] = lomb(t,ipi,OFAC, HIFAC);   % estimate power spectrum

% find peaks and their significance in the raw spectra
[peak.amp, peak.loc] = findpeaks(spec.P);             % find spectral peaks
peak.prob = spec.p(peak.loc);                         % determine their p-value
peak.significant = spec.F(peak.loc(peak.prob<alpha)); % get frequencies of significant peaks

% interpolate spectra to period-axis common to all flies
a.F = F;
a.spec = interp1(spec.F, spec.P, a.F);
