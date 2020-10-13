function gfh=gamfitwithheight(X,graphics,onset_delay,rise_decay_rate,height,fitcol)
% does an amoeba-style modified gamma variate fit to a vector X
% with parameters (onset-delay, rise-decay-rate, height)
% X=gammavar(1:30,2,7,.5); X(5)=.1; X(20)=.4; 
% e.g.,  gfh=gamfitwithheight(X)

if nargin<2, graphics=1; end
if nargin<3, onset_delay=2; end
if nargin<4, rise_decay_rate=7; end
if nargin<5, height=1; end
if nargin<6, fitcol='r'; end

global gfh;
gfh.X=X;
options=optimset('Display','off');
gfh.parms=fminsearch(@gammavarmatch,[onset_delay rise_decay_rate height],options);
gfh.mse=gammavarmatch(gfh.parms);
gfh.pred=gammavar(1:length(X),gfh.parms(1), gfh.parms(2), gfh.parms(3));
if size(X,1)~=size(gfh.pred,1)
  gfh.pred=gfh.pred';
end
gfh.r=r(X,gfh.pred);
gfh.parmlabels={'alpha - offset','beta - rise/decay time', 'height'};

if max(gfh.pred)==0
  gfh.X=-X;
  gfh.parms=fminsearch(@gammavarmatch,[onset_delay rise_decay_rate height]);
  gfh.mse=gammavarmatch(gfh.parms);
  gfh.pred=gammavar(1:length(X),gfh.parms(1), gfh.parms(2), gfh.parms(3));
  if size(X,1)~=size(gfh.pred,1)
    gfh.pred=gfh.pred';
  end
  gfh.r=r(X,gfh.pred);
  gfh.X=X;
  gfh.pred=-gfh.pred;
  gfh.parms(end)=-gfh.parms(end);
end

if graphics
  %figure(1); clf;
  plot(X);
  hold on
  plot(gfh.pred,fitcol);
  legend('Data','Predicted');
  title(sprintf('Gamma variate fit, r=%.2f',gfh.r));
end
