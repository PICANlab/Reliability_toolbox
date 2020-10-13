function f=gammavar(t,alpha,beta,height)
% usage: f=gammavar(t,alpha,beta,height)
% basic form from:  http://pet.utu.fi/staff/vesoik/reports/tpcmod0010.pdf
% I added the height option for normalization
% where height =1 makes the height of the top point be 1 (i.e., unit
% height). height=0 is unconstrained.
% alpha is 
% beta affects the rise and decay. High beta is slow rise, slow
% decay.
% alpha is the offset. high alpha starts later.
% e.g., plot(gammavar(0.01:.05:10))
%  plot(gammavar(1:30,2,7,.5))
if nargin<2, alpha=1; end
if nargin<3, beta=1; end
if nargin<4, height=1; end

f=(t.^alpha.*exp(-t./beta))./(beta.^(alpha+1).*gamma(alpha+1));
if height, f=height.*f./max(f); end


%clf; hold on;
%for alpha=1:5
%plot(gammavar(0.01:.05:10,alpha,1,1))
%end
% with t out to 30, alpha and beta can both be up to 3.


