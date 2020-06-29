function [ResF,ResS]=photometric_astrometric_model(T,f1,varargin)
% Calculate the the total flux and centeroid from an individaul light curve
% Package: +TimeDelay
% Input  : - A vector of times.
%          - A vector of the source (first image) light curve.
%          * Pairs of ...,key,val,... The following keys are available:
%            'A0'
%            'A'
%            'Tau'
%            'x0'
%            'y0'
%            'x'
%            'y'
%            'eps_x'
%            'eps_F'
%            'gamma'
%            'InterpMethod'
% Example: T = (1:1:240).';
%          TS = Util.stat.rand_ps(T,[2.5 1],zeros(size(T)));
%          ResA=TimeDelay.photometric_astrometric_model(T,TS(:,2))


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
F = InPar.A0 + sum(f.*InPar.A(:).',2);
X = (InPar.A0.*InPar.x0 + sum(InPar.A(:).'.*InPar.x(:).'.*f,2)) ./ F;
Y = (InPar.A0.*InPar.y0 + sum(InPar.A(:).'.*InPar.y(:).'.*f,2)) ./ F;



% truncate time series to overlap region
Flag = ~isnan(F);
T    = T(Flag);
F    = F(Flag);
X    = X(Flag);
Y    = Y(Flag);

Nt   = numel(T);
F    = F + InPar.eps_F.*ones(Nt,1);
X    = X + InPar.eps_x.*ones(Nt,1);
Y    = Y + InPar.eps_x.*ones(Nt,1);

ResF.T       = T;
ResF.f       = f;
ResF.F       = F;
ResF.X       = X;
ResF.Y       = Y;


if nargout>1

    Gx = X.*F;
    Gy = Y.*F;

    DeltaT = min(diff(T));
    TimeSpan = range(T);
    FreqVec = (-0.5./DeltaT:1./TimeSpan:0.5./DeltaT).';
    Omega = 2.*pi.*FreqVec;

    F_w  = FT(F,T,Omega);
    Gx_w = FT(Gx,T,Omega);
    Gy_w = FT(Gy,T,Omega);

    Ax_w = Gx_w./F_w;
    Ay_w = Gy_w./F_w;


    ResS.FreqVec = FreqVec;
    ResS.Omega   = Omega;

    ResS.Gx      = Gx;
    ResS.Gy      = Gy;
    ResS.Gx_w    = Gx_w;
    ResS.Gy_w    = Gy_w;
    ResS.Ax_w    = Ax_w;
    ResS.Ay_w    = Ay_w;
end





