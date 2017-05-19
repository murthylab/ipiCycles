cc()
addpath('code')
% files with pulsetimes/intervals of different datasets and annotations
fileNames = {'Stern2014_FSSCoen2014.mat';...
   'Stern2014_FSSStern2014.mat';...
   'Stern2014_KyriacouManual2017.mat';...
   'Stern2017_FSSCoen2014.mat'}

%%
ipiCutoffLow = 15/1000;    % [ms] lower cutoff for ipis - anything shorter is likely segmentation error
ipiCutoffHigh = [55 75]/1000;   % [ms] upper cutoff for ipis - anything longer will be removed from analysis
OFAC = 20;        % Lomb Scargle oversampling factor (see help lomb)
HIFAC = 0.4;      % Lomb Scargle (see help lomb)
alpha = 0.05;     % significance level for peak detection
F = 1./(20:.1:150);
dK = load(['data/' 'Stern2014_KyriacouManual2017.mat']); %  manual annotation by Kyriacou et al. (2017)
cutLabel = {'full','cut'}; % cut recordings to the part that was manually annotated manually or use full recording

for fil = 1:length(fileNames) % for each annotation - manual Kyriacou et al. 2017, automatic with parameters from Stern et al. (2014) and Coen et al. (2014)
   disp(fileNames{fil})
   d = load(['data/' fileNames{fil}]); % load results from all flies for the current segmentation
   clear a spec peak
   for ic = 1:length(ipiCutoffHigh) % cut off IPIs at 55 or 75 ms
      for ct = 1:2  % cut recordings to the part that was manually annotated manually or use full recording
         if ct==2 && fil>2 % do not cut the manually annotated data set or the new data set for which no reference annotation exists
            warning('no cutting')
            continue
         end
         for fly = 1:size(d.ipi,2) % for each fly
            
            try
               disp(d.flyNames{fly})
               ipi = d.ipi(:,fly);        % ipis in seconds
               t = d.pulseTimes(:,fly);   % pulses times in seconds
               
               [t, ipi] = cleanData(t, ipi, ipiCutoffLow, ipiCutoffHigh(ic));
               
               if ct == 2 % cut recording to include only the song manually annotated by Kyriacou
                  idx = find(endsWith(dK.flyNames, d.flyNames{fly}), 1, 'first'); % find current fly in man$
                  if ~isempty(idx)
                     tStart = nanmax(0,ceil(nanmax(0, nanmin(dK.pulseTimes(:,idx)))));
                     tEnd = floor(nanmax(dK.pulseTimes(:,idx))+1);
                     
                     toCut = t<tStart | t>tEnd;
                     t(toCut) = [];
                     ipi(toCut) = [];
                  else
                     t = [];
                     ipi = [];
                     warning('no manual annotation found - skipping')
                     continue
                  end
               end
               
               %% calculate spectra
               [thisSpec, thisPeak, thisA] = ipiSpectra(t,ipi, OFAC, HIFAC, alpha, F);
               % reassemble data
               spec.F{fly} = thisSpec.F; spec.P{fly} = thisSpec.P; spec.p{fly} = thisSpec.p;
               peak.amp{fly} = thisPeak.amp; peak.loc{fly} = thisPeak.loc;
               peak.prob{fly} = thisPeak.prob; peak.significant{fly} = thisPeak.significant;
               a.F = thisA.F; a.spec(:,fly) = thisA.spec; a.flyNames{fly} = d.flyNames{fly};
            catch ME
               disp(ME.getReport())
            end
         end
         saveFileName = sprintf('spectra/%s_spec_ipiCutoff%dms_%s', fileNames{fil}(1:end-4), ipiCutoffHigh(ic)*1000, cutLabel{ct});
         fprintf('saving to %s.\n', saveFileName)
         save(saveFileName, 'spec', 'peak', 'a')
      end
   end
end