function ResF=timedelayed_eqspaced_lc(Nt,varargin)
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
% Example: ResF=TimeDelay.timedelayed_lc;


InPar = inputParser;
addOptional(InPar,'Dt',1);
addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 0.66]); %2./3]);
addOptional(InPar,'Tau',[14.7]);   % all positive!
addOptional(InPar,'x0',0);  
addOptional(InPar,'y0',0);  
addOptional(InPar,'x',[1 -1]);  
addOptional(InPar,'y',[0 0]);
addOptional(InPar,'StdF',0.3);    % Std amplitude relative to mean(F) 
addOptional(InPar,'eps_x_abs',0.05);  
addOptional(InPar,'eps_F_rel',0.025);
addOptional(InPar,'t_vec',(1:1:1024).');
addOptional(InPar,'gamma',3);
addOptional(InPar,'InterpMethod','pchip');
parse(InPar,varargin{:});
InPar = InPar.Results;


Freq = TimeDelay.fft_freq(Nt,InPar.Dt);
Freq = Freq.';
w    = 2.*pi.*Freq;

Amp  = abs(w).^(-0.5.*InPar.gamma);
f1_w = Amp.*randn(Nt,1); +1i.*Amp.*randn(Nt,1);
f_DC = 200;
f1_w(1) = Nt.*f_DC;
f1_t = ifft(f1_w,'symmetric');
f1_w = fft(f1_w);

phi_w = sum((InPar.A(1) + InPar.A(2:end).*exp(1i.*w.*InPar.Tau)).*f1_w,2);
phi_w(1) = phi_w(1) + Nt.*InPar.A0;
phi_t = ifft(phi_w);
chi_t = sum((InPar.A(1).*InPar.x(1) + InPar.A(2:end).*InPar.x(2:end).*exp(1i.*w.*InPar.Tau)).*f1_w,2)./phi_t;




if nargin==0
    % simulation mode
    T = InPar.t_vec;
    TS = TimeDelay.simulate_lc_from_ps(T,[InPar.gamma 1],zeros(size(T)));
    f1 = TS(:,2);
    
    f1 = f1 + std(f1)./InPar.StdF; % add DC to simulated LC
end

if nargin>0 && ~isempty(InPar.StdF)
    f1 = f1 + std(f1)./InPar.StdF;
end


%FT = @(Flux,Time,Omega)  sum(Flux(:).*exp(-1i.*Omega(:).'.*Time(:)),1).'; %./2;



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

eps_F_abs = InPar.eps_F_rel.*mean(F_t);

Nt   = numel(T);
F_t    = F_t +       eps_F_abs.*randn(Nt,1);
x_t    = x_t + InPar.eps_x_abs.*rand(Nt,1);
y_t    = y_t + InPar.eps_x_abs.*randn(Nt,1);

ResF.T       = T;
ResF.f       = f;
ResF.F_t     = F_t;
ResF.x_t     = x_t;
ResF.y_t     = y_t;
ResF.eps_F_abs = eps_F_abs;
ResF.eps_x_abs = InPar.eps_x_abs;


