function [a,z] = significance(d,t)

[~,~,~,ofac,~,f,~] = init(d,t);
nf = length(f);

%% Significance
M = 2*nf/ofac;
a = [0.001 0.01 0.05];%[0.005  0.1 0.5]
z = -log(1-(1-a).^(1/M));
end


%% init (initialize)
function [x,t,hifac,ofac,a,f,fig] = init(x,t,varargin)
if nargin < 6, a = [];    % set default value for a
else           a = sort(varargin{4});
    a = a(:)';
end
if nargin < 5, ofac = 4;  % set default value for ofac
else           ofac = varargin{3};
end
if nargin < 4, hifac = 1; % set default value for hifac
else           hifac = varargin{2};
end
if nargin < 3, fig = 0;   % set default value for hifac
else           fig = varargin{1};
end

if isempty(ofac),  ofac  = 4; end
if isempty(hifac), hifac = 1; end
if isempty(fig),   fig   = 0; end

if ~isvector(x) ||~isvector(t),
    error('%s: inputs x and t must be vectors',mfilename);
else
    x = x(:); t = t(:);
    nt = length(t);
    if length(x)~=nt
        error('%s: Inputs x and t must have the same length',mfilename);
    end
end

[t,ind] = unique(t);    % discard double entries and sort t
x = x(ind);
if length(x)~=nt, disp(sprintf('WARNING %s: Double entries have been eliminated',mfilename)); end

T = t(end) - t(1);
nf = round(0.5*ofac*hifac*nt);
f = (1:nf)'/(T*ofac);
end