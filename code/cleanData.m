function [t, ipi] = cleanData(t, ipi, ipiCutoffLow, ipiCutoffHigh)
% [t, ipi] = cleanData(t, ipi, ipiCutoffLow, ipiCutoffHigh)
% Removes NaN values, too short and too long IPIs.
% 
% PARAMS
%  t        - vector of pulse times in seconds
%  ipi      - vector pulse intervals in seconds
%  ipiCutoffLow  - smallest ipi value to retain
%  ipiCutoffHigh - greates ipi value to retain
%
% RETURNS
%  t, ipi - cleaned-up inputs

% vectors are padded with nan - remove them here
badIdx = isnan(ipi) | isnan(t);
t(badIdx) = [];
ipi(badIdx) = [];

% remove short ipis - these are likely segmentation errors
t(ipi<ipiCutoffLow) = [];
ipi(ipi<ipiCutoffLow) = [];

% removes long ipis
t(ipi>ipiCutoffHigh) = [];
ipi(ipi>ipiCutoffHigh) = [];
         