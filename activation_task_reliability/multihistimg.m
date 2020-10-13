function multihistimg(mat,numbins,rowlabels,skipfirst)
if nargin<2, numbins=100; end
if nargin<3, rowlabels=''; end
if nargin<4, skipfirst=1; end
clf;
[nh,bins]=hist(mat,numbins);
clf;
if skipfirst
  imagesc(nh(2:end,:)');
else
  imagesc(nh');
end
histvals=reshape(nh,prod(size(nh)),1);
set(gca,'YTick',1:length(rowlabels));
set(gca,'YTickLabel',rowlabels);
ticks=get(gca,'XTick');
tickvals=gsresample(bins,length(bins),length(ticks));
tickvals=round(tickvals.*100)./100;
set(gca,'xTickLabel',tickvals);
colormap(min(1,hot+.3));
%caxis([prctle(histvals,10) prctile(histvals,90)]);
colorbar;
