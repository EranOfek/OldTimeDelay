function ResF=timedelayed_lc(T,f1,varargin)
% Given a light curve, generate a linear combination of time delayed light curves.
% Package: +TimeDelay

% Input  : - A vector of times.
%          - A vector of the source (first image) light curve.
%          * Pairs of ...,key,val,... The following keys are available:
%            'A0' - \alpha_0 
%            'A'  - Vector of \alpha_i
%            'Tau' - vector of \tau (skipping tau_1=0)
%            'x0'  - x_0. Default is 0.
%            'y0'  - y_0. Default is 0.
%            'x'   - Vector of x_i.
%            'y'   - vector of y_i
%            'eps_x'
%            'eps_F'
%            'gamma'
%            'InterpMethod'
% Example: T = (1:1:240).';
%          TS = TimeDelay.simulate_lc_from_ps(T,[2.5 1],zeros(size(T)));
%          ResA=TimeDelay.timedelayed_lc(T,TS(:,2))


InPar = inputParser;
addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 0.4]); %2./3]);
addOptional(InPar,'Tau',[36.9]);   % all positive!
addOptional(InPar,'x0',0);  
addOptional(InPar,'y0',0);  
addOptional(InPar,'x',[0 1]);  
addOptional(InPar,'y',[0 0]);  
addOptional(InPar,'eps_x',0.02);  
addOptional(InPar,'eps_F',0.05);  
addOptional(InPar,'gamma',2.5);
addOptional(InPar,'InterpMethod','pchip');
parse(InPar,varargin{:});
InPar = InPar.Results;


FT = @(Flux,Time,Omega)  sum(Flux(:).*exp(-1i.*Omega(:).'.*Time(:)),1).'; %./2;



% add first image time delay
InPar.Tau = [0, InPar.Tau];

Nim = numel(InPar.Tau);
T   = T(:);
f1  = f1(:);

% generate shifted individual LC
f = zeros(numel(T),Nim);
for Iim=1:1:Nim
    f(:,Iim) = interp1(T,f1,T-InPar.Tau(Iim),InPar.InterpMethod);
end

% truncate the edges
Tcut = ceil(max(InPar.Tau));
f  = f(Tcut:end,:);



% total flux
F_t = InPar.A0 + sum(f.*InPar.A(:).',2);
x_t = (InPar.A0.*InPar.x0 + sum(InPar.A(:).'.*InPar.x(:).'.*f,2)) ./ F_t;
y_t = (InPar.A0.*InPar.y0 + sum(InPar.A(:).'.*InPar.y(:).'.*f,2)) ./ F_t;



% truncate time series to overlap region
Flag = ~isnan(F_t);
T    = T(Flag);
F_t  = F_t(Flag);
x_t  = x_t(Flag);
y_t  = y_t(Flag);

Nt   = numel(T);
F_t    = F_t + InPar.eps_F.*ones(Nt,1);
x_t    = x_t + InPar.eps_x.*ones(Nt,1);
y_t    = y_t + InPar.eps_x.*ones(Nt,1);

ResF.T       = T;
ResF.f       = f;
ResF.F_t     = F_t;
ResF.x_t     = x_t;
ResF.y_t     = y_t;


