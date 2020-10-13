function [mr]=mregs(X,Y,outlier,getts,prettyprint)
% performs a multiple regression
% where X is a matrix of column vectors and Y is a column vector
% Useage: [mr]=mregs(X,Y[,outlier,getts,prettyprint]) 
% where mr is a structure with relevant values
%  including Rsq,B,B0,Ypred,stb,dfr,dff,F,p
%    where stb is standardized betas
% if getts is true, t and p values are also returned for each beta

if (nargin<3) | (isempty(outlier)) | (length(outlier)==1), outlier=zeros(size(Y)); end; 
if nargin<4, getts=0; end
if nargin<5, prettyprint=0; end

if getts
  [mr.Rsq, mr.B,mr.B0, mr.Ypred, mr.stb, mr.t, mr.tp]=mreg(X,Y,outlier,getts);
else
  [mr.Rsq,mr.B,mr.B0,mr.Ypred,mr.stb]=mreg(X,Y,outlier);
end
if size(X,2)==1
  mr.r=sqrt(mr.Rsq).*sign(mr.B(end));
end

ytouse=find(~outlier);
mr.meanPointsOff=mean(abs(Y(ytouse)-mr.Ypred));


Er=1; 
Ef=1-mr.Rsq;
mr.dfr=size(X,1)-length(find(outlier))        -1;
mr.dff=size(X,1)-length(find(outlier))-size(X,2)-1;
mr.F=((Er-Ef)./(mr.dfr-mr.dff))./(Ef./mr.dff);
mr.p=1-fcdf(mr.F,mr.dfr-mr.dff,mr.dff);

% ncols=size(X,2);
% if (getts & (ncols>1))
%   if ncols==2
%     rsqs(1)=mr.Rsq-mreg(X(:,2),Y,outlier);
%     rsqs(2)=mr.Rsq-mreg(X(:,1),Y,outlier);
%   else
%     for ct=1:ncols
%       mincol=1; if ct==1, mincol=2; end
%       maxcol=ncols; if ct==ncols, maxcol=ncols-1; end
%       Xs=[X(:,mincol:ct-1) X(:,ct+1:maxcol)];
%       rsqs(ct)=mr.Rsq-mreg(Xs,Y,outlier);
%     end
%   end
%   denom=sqrt((1-mr.Rsq)./mr.dff);
%   for ct=1:ncols
%     mr.t(ct)=rsqs(ct)./denom;
%     mr.tp(ct)=(1-tcdf(mr.t(ct),mr.dff))./2; % divide by 2 for 2-tailed
%   end
% end
if prettyprint
  if size(X,2)==1
    fprintf(1,'r=%.3f, Rsq = %.3f, F(%d,%d)=%.3f, p=%.3f, M(points off)=%.3f\n',mr.r, mr.Rsq, mr.dfr-mr.dff,mr.dfr,mr.F,mr.p,mr.meanPointsOff);
  else
    fprintf(1,'Rsq = %.3f, F(%d,%d)=%.3f, p=%.3f, M(points off)=%.3f\n',mr.Rsq, mr.dfr-mr.dff,mr.dfr,mr.F,mr.p,mr.meanPointsOff);
  end
  if (getts & prettyprint>1)
    for ct=1:length(mr.t)
      if ct>1
	fprintf(' B%d=%.3f, stB=%.3f, t(%d)=%.3f, p=%.3f\n',ct-1,mr.B(ct),mr.stb(ct-1),length(Y)-1,mr.t(ct), mr.tp(ct));
      else
	fprintf(' B%d=%.3f, stB=----, t(%d)=%.3f, p=%.3f\n',ct-1,mr.B(ct),length(Y)-1,mr.t(ct), mr.tp(ct));
      end
    end
  end
end

    
