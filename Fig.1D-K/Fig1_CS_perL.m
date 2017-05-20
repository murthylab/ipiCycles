% compares spectral peaks (irrepsecitve of their significance) of CantonS wild type and perL mutant flies for 
% different period ranges and w/ and w/o false discovery rate correction
addpath(genpath('../code'))
cc()
Fs = 10000;
strLabel = {'Kyriacou manual', 'Stern new'};

figure('Name', 'Fig. 1D-K - genotype-specific periodicity')
for ii = 1:2
   clearvars -except ii Fs strLabel
   cnt = 0;
   if ii == 1
      load('../spectra/Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_full.mat')
   else
      load('../spectra/Stern2017_FSSCoen2014_spec_ipiCutoff75ms_full.mat')
   end
   %%
   Fb = 1:10:length(a.F);
   G = grp2idx(contains(a.flyNames, 'CS'));
   Glabel = {'perL','CS'};
   for fly = 1:size(a.spec,2) % individual flies
      peakFreq(fly) = a.F(argmax(a.spec(:,fly)));
      for f1 = 1:length(Fb)
         for f2 = 1:length(Fb)
            if f2>f1
               peakFreqRange(fly,f1,f2) = a.F(f1-1+argmax(a.spec(Fb(f1):Fb(f2),fly)));
            else
               peakFreqRange(fly,f1,f2) = nan;
            end
         end
      end
   end
   
   %%
   mySubPlot(2,7,ii,1)
   cla
   hold on
   [M, QL, QU] = grpstats(1./peakFreq, G, {@mean, @quartile25, @quartile75});
   plot([1 2]+0.25, M, '.r', 'MarkerSize', 32)
   plot([1 2; 1 2]+0.25, [QL QU]', 'r')
   
   for g = 1:max(G)
      jitterXvalues(plot(g, 1./peakFreq(G==g),'ok'),0.05)
   end
   [hn1, pn1] = jbtest(1./peakFreq(G==1));
   [hn2, pn2] = jbtest(1./peakFreq(G==2));
   disp([hn1, pn1, hn2, pn2])
   ylabel('peak period [s]')
   [~, pttest] = ttest2( 1./peakFreq(G==1),  1./peakFreq(G==2), 'tail','right');
   title(sprintf('%1.1fs (perL):%1.1fs (CS),\nranksum p=%1.2f, right-tailed ttest p=%1.2f', ...
      grpstats(1./peakFreq, G), ...
      ranksum(1./peakFreq(G==1), 1./peakFreq(G==2)), ...
      pttest))
   set(gca, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel',Glabel, 'YLim', [20 150])
   
%    % frequency
%    mySubPlot(2,6,ii,2)
%    hold on
%    [M, QL, QU] = grpstats(peakFreq, G, {@mean, @quartile25, @quartile75});
%    plot([1 2]+0.25, M, '.r', 'MarkerSize', 32)
%    plot([1 2; 1 2]+0.25, [QL QU]', 'r')
%    for g = 1:max(G)
%       jitterXvalues(plot(g, peakFreq(G==g),'ok'),0.05)
%    end
%    [hn1, pn1] = jbtest(peakFreq(G==1));
%    [hn2, pn2] = jbtest(peakFreq(G==2));
%    disp([hn1, pn1, hn2, pn2])
%    
%    ylabel('peak frequency [Hz]')
%    [~, pttest] = ttest2( peakFreq(G==1),  peakFreq(G==2), 'tail','left');
%    title(sprintf('%1.3fHz (perL):%1.3fHz (CS),\nranksum p=%1.2f, left-t. ttest p=%1.2f', ...
%       grpstats(peakFreq, G), ...
%       ranksum(peakFreq(G==1), peakFreq(G==2)), ...
%       pttest))
%    set(gca, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel',Glabel, 'YLim', [0 0.06])
%    
   %% now look at dependence on freq range considered
   clear m e
   for g = 1:max(G)
      m(:,:,g) = nanmean(1./peakFreqRange(G==g,:,:),1);
      e(:,:,g) = sem(1./peakFreqRange(G==g,:,:),1);
   end
   %%
   clear p pL
   for f1 = 1:length(Fb)
      for f2 = f1+1:length(Fb)
         if f2>f1
            [~, p(f1,f2), ci, stats] = ttest2(1./peakFreqRange(G==1,f1,f2), 1./peakFreqRange(G==2,f1,f2), 'tail','right');
            t(f1,f2) = stats.tstat;
            [~, pL(f1,f2), ci, stats] = ttest2(1./peakFreqRange(G==1,f1,f2), 1./peakFreqRange(G==2,f1,f2), 'tail','left');
         else
            p(f1,f2) = nan;
            p(f1,f2) = nan;
         end
      end
   end
   %%
   meanRatio = log2(m(:,:,1)./m(:,:,2));
   alpha = 0.05;
   Ntests = sum(p(:)>0);
   
   pCorr = fdr(p,1);
   pRatio = log2(p/alpha);
   pRatio(isinf(pRatio)) = nan;
   pRatioCorr = log2(pCorr/alpha);
   pRatioCorr(isinf(pRatioCorr)) = nan;
   
   pLCorr = fdr(pL,1);
   pLRatio = log2(pL/alpha);
   pLRatio(isinf(pLRatio)) = nan;
   pLRatioCorr = log2(pLCorr/alpha);
   pLRatioCorr(isinf(pLRatioCorr)) = nan;

   %%
   cmap = flipud(cbrewer('Div','RdBu',64));
   cmap(1,:) = [1 1 1];
   
   hP(1) = mySubPlot(2,7,ii,3:4);
   imagesc(1./a.F(Fb), 1./a.F(Fb), pRatio')
   [X,Y] = meshgrid(1./a.F(Fb(1:end-1)), 1./a.F(Fb));
   v = [0, 0];
   hold on
   contour(X,Y, pRatio',v, 'k', 'LineWidth', 2)
   xlabel('lower bound on period [s]')
   ylabel('upper bound on period [s]')
   colormap(gca, cmap)
   
   hcb = colorbar('east');
   hcb.Position(3:4) = hcb.Position(3:4)/2;
   title(hcb, 'p-value')
   set(hcb, 'Ticks', -4:2:4, 'TickLabels', round((2.^(-4:2:4)*0.05)*1000)/1000)
   title({'p-values for right-tailed ttest' '(perL>CS), contour at q=0.05'})
   
   hP(2) = mySubPlot(2,7,ii,5:6);
   imagesc(1./a.F(Fb), 1./a.F(Fb), pRatioCorr')
   [X,Y] = meshgrid(1./a.F(Fb(1:end-1)), 1./a.F(Fb));
   v = [0, 0];
   hold on
   contour(X,Y, pRatioCorr',v, 'k', 'LineWidth', 2)
   colormap(gca, cmap)
   title({'q-values (BH method) for right-tailed ttest' '(perL>CS), contour at q=0.05'})
   xlabel('lower bound on period [s]')
   fprintf('significant cells\n   right-tailed %1d\n   left-tailed %1d\n',...
      sum(pRatio(:)>0 & pRatio(:)<0.05), ...
      sum(pLRatio(:)>0 & pLRatio(:)<0.05))
   
   axis(hP, 'square', 'xy')
   set(hP, 'CLim', [-4 4], 'XTick', 25:25:150, 'YTick', 25:25:150)
   
   %%
   % retain only the values for DeltaF>20 and vectorize
   p = triu(p,20); p = p(p(:)>0);
   pCorr = triu(pCorr,20);pCorr = pCorr(pCorr(:)>0);
   pL = triu(pL,20);pL = pL(pL(:)>0);
   pLCorr = triu(pLCorr,20);pLCorr = pLCorr(pLCorr(:)>0);
   % calculate fraction of significant frequency ranges
   fractionSignificant = [mean(p(:)>0  & p(:)<0.05),  mean(pCorr(:)>0  & pCorr(:)<0.05),...
                          mean(pL(:)>0 & pL(:)<0.05), mean(pLCorr(:)>0 & pLCorr(:)<0.05)];
   % plot
   mySubPlot(2,7,ii,7);
   bar(fractionSignificant, 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6])
   set(gca, 'YLim', [0 1], 'XLim', [0 5],...
      'XTick', 1:4, 'XTickLabel', {'perL>CS, p' 'q' 'perL<CS, p' 'q' }, 'XTickLabelRotation', 60);
   ylabel('fraction of significant frequency ranges with \DeltaF>20Hz (p<0.05 or q<0.05)')
   drawnow
end
clp()
