function gmse=gammavarmatch(parms)
% matches data and parameters to a gamma variate 
% for use with gamfitwithheight(X)

alpha=parms(1);
beta=parms(2);
height=parms(3);

global gfh;
t=1:length(gfh.X);
f=gammavar(t,alpha,beta,height);
if size(f,1)~=size(gfh.X,1)
    f=f';
end
gmse=mse(f,gfh.X);

