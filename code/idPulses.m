function [confMat, eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesCell, tolerance)
% Identifies the pulse events across different annotations and computes
% confusion matrices 
% [confMat, eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesCell, tolerance)
%
% PARAMS
% pulseTimesCell - cell array of event times corresponding to different
%                  annotations
% tolerance      - tolerance within which to assign different times to the
%                  same event
%
% RETURNS
% confMat, eventMat, pulseId, pulseTimes, pulseGroup

if ~exist('tolerance','var') || isempty(tolerance)
	tolerance = 5/1000;% 5ms
end

% pool and sort all pulses
pulseTimes = [];
pulseGroup = [];
for pl = 1:length(pulseTimesCell)
   pulseTimes = [pulseTimes; pulseTimesCell{pl}];
   pulseGroup = [pulseGroup; pl*ones(size(pulseTimesCell{pl}))];
end
sortIdx = argsort(pulseTimes);
pulseTimes = pulseTimes(sortIdx);
pulseGroup = pulseGroup(sortIdx);

% go through pulses and assign those within tolerance to same event
cnt = 1;
pulseId = nan(size(pulseTimes));
pulseId(1) = 1;
for pul = 2:length(pulseTimes)
   if pulseTimes(pul)-tolerance>pulseTimes(pul-1)
      cnt = cnt+1;
   end
   pulseId(pul) = cnt;
end

% build binary event matrix for each set of pulse times
ids = max(pulseId);
eventMat = zeros(ids, max(pulseGroup));
for g = 1:max(pulseGroup)
   eventMat(pulseId(pulseGroup==g), g) = 1;
end
%% generate confusion matrices
confMat = zeros(2, 2, 1, 1);
for g1 = 1:max(pulseGroup)
   for g2 = 1:max(pulseGroup)
      [tmp, gn] = confusionmat(eventMat(:,g1), eventMat(:,g2));
      if length(tmp)==1
          confMat(1,1,g1,g2) = tmp;
      else
          confMat(:,:,g1,g2) = tmp;
      end
   end
end