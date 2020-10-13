function iccv = icc(mat,singlemeas)
% given a matrix, computes the interclass correlation
% among the COLUMNS of the matrix
% e.g., where matrix is subjects going down and time across

if nargin<2, singlemeas=1; end

if singlemeas
  k=size(mat,2); % # things across which we're calculating reliability
  N=size(mat,1); % # observations
  W=N;
  Pj=sum(mat');
  Ti=sum(mat);
  GrandSum=sum(sum(mat));
  ssb=sum((Pj.^2)./k)-(GrandSum.^2)./(k.*N); 
  msb=ssb./(N-1) ;
  matsq = mat.^2;
  ssw=sum(sum(matsq))-sum((Pj.^2)./k);
  msw=ssw./(N.*(k-1));
  iccv=(msb-msw)./(msb+(k-1).*msw);

% twodigitaccuracy but makes biased estimates using k vs. k-1 for msw.
%   k=size(mat,2); % # things across which we're calculating reliability
%   N=size(mat,1); % # observations
%   W=N;
%   Pj=sum(mat')
%   Ti=sum(mat);
%   GrandSum=sum(sum(mat))
%   ssb=sum((Pj.^2)./k)-(GrandSum.^2)./(k.*N); % correct
%   msb=ssb./(N) % correct
%   matsq = mat.^2;
%   ssw=sum(sum(matsq))-sum((Pj.^2)./k);
%   msw=ssw./(N.*(k))
%   iccv=(msb-msw)./(msb+(k).*msw);

% unused - was never quite right...
%  s=repmeas2s(mat);
%  iccv=s.icc_singlemeas;
else
  % compute the average measure inter-item reliability
  iccv=alpha(mat);
end
