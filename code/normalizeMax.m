function data = normalizeMax(data)
%NORMALIZE DATA, such that max=1, but preserves offset (min)
%   data = normalizeMax(data)
%
%ARGS
%   data
%RETURNS
%   normalized data
%
% created 06/09/19 Jan
if (size(data,1)==1 | size(data,2)==1)
   data = data/(nanmax(data)+eps);
else
   for cat = 1:size(data,1)
      data(cat,:) = data(cat,:)/(nanmax(data(cat,:))+eps);
   end
end
