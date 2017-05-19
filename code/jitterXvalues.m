function varargout = jitterXvalues(hh, sigma)

if ~exist('sigma','var') || isempty(sigma)
   sigma = 0.03;
end

if isgraphics(hh, 'axes')
   lineHandles = hh.Children;
elseif isgraphics(hh, 'line') 
   lineHandles = hh; 
end

for child = 1:length(lineHandles)
   XData = lineHandles(child).XData;
   lineHandles(child).XData = XData+randn(size(XData))*sigma;
end

if nargout>0
   varargout{1} = hh;
end
