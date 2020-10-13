function [Rsq,B,B0,Ypred,STB,ts,ps]=mreg(X,Y,outlier,getts)
% Usage: [Rsq,B,B0,Ypred,STB,[ts,ps]]=mreg(X,Y[,outlier,getts]) 
% Performs a multiple regression
%  X is a matrix of column vectors representing 
%    explantory (independent) variables
%  Y is a column vector representing the 
%    response (dependent) variable
%  outlier is an optional column vector the same size as Y 
%    in which entries with zero are used and non-zero entries
%    are censored
%  if getts is true t and p values are calculated for all betas
% by Greg Siegle

if (nargin<2)
  fprintf(1,'Useage: [Rsq,B,B0,Ypred]=mreg(X,Y)\n');
  return
end

if nargin<3, outlier=[]; end
if nargin<4, getts=0; end

if ~isempty(outlier) % restrict X and Y based on outlier
  indices=find(~outlier);
  nx=zeros(length(indices),size(X,2));
  ny=zeros(length(indices),1);
  for ind=1:length(indices)  
    nx(ind,:)=X(indices(ind),:);
    ny(ind)=Y(indices(ind));
  end
  X=nx;
  Y=ny;
end


mx=mean(X);
for col=1:size(X,2)
  Xn(:,col)=X(:,col)-mx(col);
end

Yn=Y-mean(Y);
B=inv(Xn'*Xn)*Xn'*Yn;
B0=mean(Y-X*B);
Ypred=X*B+B0;
Rsq=var(Ypred)/var(Yn);
for col=1:size(X,2)
  STB(:,col)=B(col)/(std(Y)/std(X(:,col)));
end
  

if ~getts
  ts=0; ps=0;
else
  B=0;
  bo = ones(length(X),1);
  Xnew =[bo,X];
  B=(inv(Xnew'*Xnew))*(Xnew'*Y);
  MSE=((Y'*Y)-(B'*Xnew'*Y))/(length(X)-size(X,2));
  sebcov=sqrt(MSE*(abs(inv(Xnew'*Xnew))));
  
  for ct=1:length(sebcov)
    seb(ct,:) = sebcov(ct,ct);
  end
  
  for ct=1:length(B)
    ts(ct,:)=B(ct,1)/seb(ct,1);
  end
  ts=abs(ts);
  
  df=length(X)-size(X,2)-1;
  
  for ct=1:length(ts)
    ps(ct,:)=2*(1-tcdf(ts(ct),df));
  end;
end



%stdX=(X-mean(X))./std(X);
%stdY=(Y-mean(Y))./std(Y);
%stdB=inv(stdX'*stdX)*stdX'*stdY;
%B0=mean(stdY-stdX*stdB);
%stdYpred=stdX*stdB+B0;
%Rsq=var(stdYpred)/var(stdY);
