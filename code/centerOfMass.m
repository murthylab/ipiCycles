function com = centerOfMass(x,y, threshold)
% computes center of mass
% USAGE: centerOfMass(x,y, threshold)

if nargin==2
   threshold =  0.5;
end
y = normalizeMax(y);          % norm to max=1
y = limit(y-threshold,0,inf); % get rid of everything subthreshold
com = x*normalizeSum(y)';     % get weighted sum of super thres portion
