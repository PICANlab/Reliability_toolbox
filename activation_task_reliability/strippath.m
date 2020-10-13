function [sp,fname]=strippath(str)
% takes in a full path and returns the 
% path and file components

strrev=fliplr(str);
firstslash=first(find((strrev=='/') | (str=='\')));
if firstslash
  sp=fliplr(strrev(firstslash:end));
  fname=fliplr(strrev(1:firstslash-1));
else
  sp=[]; fname=str;
end
