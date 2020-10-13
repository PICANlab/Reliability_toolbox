function mnew=rescaleoutliers(m,threshold)
% usage: mnew=rescaleoutliers(m,threshold)
% given a matrix Windzorizes each COLUMN
% That is, outliers, defined as outside 
%  25th prctile - threshold*IQR
%  75th prctile + threshold*IQR
% are rescaled to the maximum in the direction they went 
% outside the valid interval.
% threshold defaults to 1.5
% NOTE:  this behaves weirdly if you have a matrix of timeseries
% for that use rescaleoutlierstimeseriesmatrix.
% by Greg Siegle


% test with m=rand(10,4); m(2,2)=17; m(8,3)=10;

if nargin<2
  threshold=1.5;
end

if ((size(m,1)==1) & (size(m,2)>1))
   fprintf('Needs columnar data. Reorienting\n');
   m=m';
end

if threshold==0
  return;
end


missing=(abs(m)==999);
% don't use missing in computing percentiles
colmed=median(m);
medmat=repmat(colmed,size(m,1),1);
mnomissing=m;
eachmissing=find(abs(m)==999);
mnomissing(eachmissing)=medmat(eachmissing);

medm=median(mnomissing);
iqrm=tsnaniqr(mnomissing);
q1=tsprctile(mnomissing,25);
q3=tsprctile(mnomissing,75);

lbm=m<repmat(q1-threshold.*iqrm,size(m,1),1);
ubm=m>repmat(q3+threshold.*iqrm,size(m,1),1);

outlier=(lbm | ubm);
mnew=(1-outlier).*m+outlier.*repmat(medm,size(m,1),1);
maxgood=max(mnew);
mingood=min(mnew);
mnew=(1-lbm).*m+lbm.*repmat(mingood,size(m,1),1);
mnew=(1-ubm).*mnew+ubm.*repmat(maxgood,size(m,1),1);
mnew=(1-missing).*mnew+missing.*-999;
