function [mr]=mhreg(X1,Xa,Y,outlier,pprint,withpartials)
% performs a multiple hierarchical regression
% where X1 is a matrix of column vectors representing the reduced model
% Xa is a matrix of additional models to make the full model
% and Y is a column vector for the dependent variable
% Useage: [mr]=mhreg(X1,Xa,Y[,outlier,pprint]) 
% where mr is a structure with r2s and dr2s.
% e.g., aaa=rand(50,4); bbb=rand(50,4); ccc=rand(50,1)+bbb(:,1)+aaa(:,1);
% mr=mhreg(aaa,bbb,ccc)

if (nargin<3)
  fprintf(1,'Useage: mr=mreg(X1,X2,Y)\n');
  fprintf(1,'where mr is a structure with r2s and dr2s\n');
  fprintf(1,'Use mreg for a straight multiple regression');
  return
end

if (nargin<4) | (length(outlier)==1), outlier=zeros(size(Y)); end; 
if nargin<5, pprint=0; end
if nargin<6, withpartials=0; end

X2=[X1 Xa]; % the full model

%[mr.r1.Rsq,mr.r1.B,mr.r1.B0,mr.r1.Ypred]=mreg(X1,Y,outlier);
%[mr.r2.Rsq,mr.r2.B,mr.r2.B0,mr.r2.Ypred]=mreg(X2,Y,outlier);
mr.r1=mregs(X1,Y,outlier);
mr.r2=mregs(X2,Y,outlier);


mr.dr.Rsq=mr.r2.Rsq-mr.r1.Rsq;

Er=1-mr.r1.Rsq; Ef=1-mr.r2.Rsq;
mr.dr.dfr=size(X1,1)-length(find(outlier))-size(X1,2)-1;
mr.dr.dff=size(X2,1)-length(find(outlier))-size(X2,2)-1;
mr.dr.F=((Er-Ef)./(mr.dr.dfr-mr.dr.dff))./(Ef./mr.dr.dff);
mr.dr.p=1-fcdf(mr.dr.F,mr.dr.dfr-mr.dr.dff,mr.dr.dff);

if withpartials
  if (size(X1,2)==1) & (size(Xa,2)==1)
    mr.dr.sr=(r(Xa,Y)-r(Xa,X1).*r(Y,X1))./sqrt(1-r(Y,Xa).^2);
    mr.dr.pr=(r(Xa,Y)-r(Xa,X1).*r(Y,X1))./(sqrt(1-r(Y,Xa).^2).*sqrt(1-r(X1,Xa).^2));
  else
    % does not do signs gracefully
    mrxay=mregs(Xa,Y,outlier);   mrxax1=mregs(Xa,X1,outlier);  
    rxay=sqrt(mrxay.Rsq); rxax1=sqrt(mrxax1.Rsq);
    ryx1=sqrt(mr.r1.Rsq); 
    mr.dr.sr=(rxay-rxax1.*ryx1)./sqrt(1-ryx1.^2);
    mr.dr.pr=(rxay-rxax1.*ryx1)./(sqrt(1-ryx1.^2).*sqrt(1-rxax1.^2));
  end
end

if pprint
   fprintf('M1R^2: %.2f, M2R^2: %.2f, chaR^2: %.2f, chaF(%d,%d)=%.2f, p=%.2f\n',mr.r1.Rsq, mr.r2.Rsq,mr.dr.Rsq,mr.dr.dfr-mr.dr.dff,mr.dr.dff,mr.dr.F,mr.dr.p);
end

