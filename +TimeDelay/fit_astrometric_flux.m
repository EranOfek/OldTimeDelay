function [Res]=fit_astrometric_flux(t,F_t,x_t,y_t,sigma_F,sigma_x,varargin)
% Fit the 2 images astrometric-flux time delay model to observations
% Package: +TimeDelay
% Description:
% Input  : - t : vector of times.
%          - F_t: vector of total flux.
%          - x_t: vector of x position.
%          - y_t: vector of y position.
%          - sigma_F: Error in flux.
%          - sigma_x: Error in position.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Solver' - Either @Util.fit.fminsearch_my | @Util.fit.fminunc_my
%                       Default is @Util.fit.fminunc_my
%            'FitPar' - The list of parameters to fit.
%                       [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%                       If NaN, then will attempt to fit the parameter.
%                       Default is 
%                       [NaN NaN NaN  NaN NaN NaN  NaN NaN NaN  3]
%            'DefPar' - The list of initial guess to use, or the parameter
%                       value if not fitted. Default is
%                       [2 1   0.66  0   1   -1   0   0   0    3]
%            'Limits' - A two column matrix of [lower, upper] bounds on the
%                       parameters. Default is
%                       [0 5;0 2;0 2;  -1 1; -2.1 2.1; -2.1 2.1;  -1 1; -2.1 2.1; -2.1 2.1;   1.5 3.5]
%            'TwoD'   - A logical indicate if to perform 2-D fit.
%                       Default is true.
%            'VecInvTau' - A vector of 1/time_delay to attempt fitting.
%                       Default is (1./100:0.5./100:1./10)
%            'Min_w'   - Minimum w. Default is 2.*pi./100.
%            'Verbose' - Default is true.
% Output : - An output structure containing the following fields:
%            .Tau 
%            .LL_H0
%            .LL_H1
%            .BestPar_H0
%            .BestPar_H1
% Example: Res=TimeDelay.fit_astrometric_flux
%          ResF=TimeDelay.timedelayed_lc;
% Res=TimeDelay.fit_astrometric_flux(ResF.T,ResF.F_t,ResF.x_t,ResF.y_t,ResF.eps_F_abs,ResF.eps_x_abs);

% [ResLC]=TimeDelay.rand_lensed;

% Res = TimeDelay.fit_astrometric_flux(ResLC.T,ResLC.F_t,ResLC.x_t,ResLC.y_t,ResLC.sigma_F_hat,ResLC.sigma_x_hat);


NPAR2D = 11 -1;
NPAR1D = 8 -1;

InPar = inputParser;
addOptional(InPar,'Solver',@Util.fit.fminunc_my);
%addOptional(InPar,'FitPar',[NaN NaN NaN  NaN NaN NaN  NaN NaN NaN  3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%addOptional(InPar,'FitPar',[0 NaN NaN  NaN   NaN   NaN   0   0   0    3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]

%addOptional(InPar,'FitPar',[0 NaN   NaN   0   0.4   -0.2    0   0   0    3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%addOptional(InPar,'DefPar',[0 1     0.66   0     0.4   -0.2  0   0   0    3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%addOptional(InPar,'Limits',[0 3; 0 5;0 5;  -1 1; -2.1 2.1; -2.1 2.1;  -1 1; -2.1 2.1; -2.1 2.1;   1.5 3.5]); %  without Tau

addOptional(InPar,'FitPar',[0 NaN   NaN      0   1   -1        3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
addOptional(InPar,'DefPar',[0 1     0.66     0   1   -1        3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
addOptional(InPar,'Limits',[0 3; 1e-5 5;0 5;  -1 1; -2.1 2.1; -2.1 2.1;     1.5 3.5]); %  without Tau

%addOptional(InPar,'FitPar',[NaN NaN NaN  NaN  NaN NaN   3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%addOptional(InPar,'DefPar',[0 1   0.66  0   1   -1       3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
%addOptional(InPar,'Limits',[0 5;0 2;0 2;  -1 1; -2.1 2.1; -2.1 2.1;    1.5 3.5]); %  without Tau

addOptional(InPar,'TwoD',false);
addOptional(InPar,'VecInvTau',(1./100:0.5./240:1./10)); %1./(20:1:33)); %1./(10.7:1:18.7)); %(2./100:0.5./100:1./10));  % sign has meaning!
addOptional(InPar,'Min_w',2.*pi./100);
addOptional(InPar,'Verbose',true);
parse(InPar,varargin{:});
InPar = InPar.Results;

InPar.FitPar = [InPar.FitPar];
InPar.DefPar = [InPar.DefPar];


if nargin==0
    %params = jsondecode(fileread('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/sims_4.json'));
    params = jsondecode(fileread('/raid/eran/projects/Algorithms/AstrometricTimeDelay/OferData/sims_4.json'));

    t   = load('/raid/eran/projects/Algorithms/AstrometricTimeDelay/OferData/t.txt')';
    F_t = load('/raid/eran/projects/Algorithms/AstrometricTimeDelay/OferData/F_t.txt')';
    x_t = load('/raid/eran/projects/Algorithms/AstrometricTimeDelay/OferData/x_t.txt')';
    N = length(F_t);

    sigma_x = params.sigma_x_prop * abs(params.x_1 - params.x_2);
    sigma_F = params.sigma_F_prop * params.f_dc;

    
    freqs = TimeDelay.fft_freq(N, params.t_step);

    gamma = params.gamma;
    sigma_F_hat = sigma_F;
    sigma_x_hat = sigma_x;
    w = 2.*pi*freqs;
    DFT = fft(eye(N), N, 1) / sqrt(N);
    DFT_dagger = DFT';
    x0 = 0;
    x1 = params.x_1;
    x2 = params.x_2;
    Alpha0 = params.alpha_0;

    Tau = params.tau;
    Alpha1 = params.alpha_1;
    Alpha2 = params.alpha_2;
    F_w = fft(F_t) / sqrt(N);

    Gx_t = x_t.*F_t;
    Gx_w = fft(Gx_t) ./ sqrt(N);

    y_t = randn(size(x_t)).*sigma_x;
    sigma_y = sigma_x;
    sigma_y_hat = sigma_y;
    
    Gy_t = y_t.*F_t;
    Gy_w = fft(Gy_t) ./ sqrt(N);

    y0 = 0;
    y1 = 0;
    y2 = 0;
    LogZ = sum(log(F_t));
   
    
else
    % Input arguments: t,F_t,x_t,y_t,sigma_F,sigma_x
    
    N = length(F_t);
    t_step = unique(diff(t));
    if numel(t_step)>1
        error('Time series must be equally spaced');
    end
    
    freqs = TimeDelay.fft_freq(N, t_step);

    sigma_F_hat = sigma_F;
    sigma_x_hat = sigma_x;
    w = 2.*pi*freqs;
    
    %x0 = 0;
    %x1 = params.x_1;
    %x2 = params.x_2;
    %Alpha0 = params.alpha_0;

    %Tau = params.tau;
    %Alpha1 = params.alpha_1;
    %Alpha2 = params.alpha_2;
    
    %F_t = F_t - mean(F_t);
    
    F_w = fft(F_t) ./ sqrt(N);
    Gx_t = x_t.*F_t;
    Gx_w = fft(Gx_t) ./ sqrt(N);
    Gy_t = y_t.*F_t;
    Gy_w = fft(Gy_t) ./ sqrt(N);

    
    sigma_y = sigma_x;
    sigma_y_hat = sigma_y;
    
end

%% Main Fitter 

N = length(F_t);

DFT = fft(eye(N), N, 1) ./ sqrt(N);
DFT_dagger = DFT';
LogZ = sum(log(F_t));
Gamma_1_ = ((DFT * diag(F_t.^2)) * DFT_dagger) * sigma_x_hat^2;

% verify the size of FitPar and DefPar
if InPar.TwoD
    if ~(numel(InPar.FitPar)==NPAR2D && numel(InPar.DefPar)==NPAR2D && size(InPar.Limits,1)==NPAR2D)
        error('For 2D fitting Limits, FitPar and DefPar must contain %d elements',NPAR2D);
    end
else
    if ~(numel(InPar.FitPar)==NPAR1D && numel(InPar.DefPar)==NPAR1D && size(InPar.Limits,1)==NPAR1D)
        error('For 1D fitting Limits, FitPar and DefPar must contain %d elements',NPAR1D);
    end
end


    
FlagN       = isnan(InPar.FitPar);

% switch func2str(InPar.Solver)
%     case 'Util.fit.fmincon_my'
%         AddPars = {[],[],[],[],InPar.Limits(FlagN,1),InPar.Limits(FlagN,2)};
%     otherwise
%         AddPars = {};
% end
AddPars = {};


BestGuessH1 = InPar.DefPar(FlagN);
BestGuessH0 = InPar.DefPar(1:2);  % A0, A1
    
% fitting for each time delay

VecTau = 1./InPar.VecInvTau; 
Ntau   = numel(VecTau);



Res.Tau    = VecTau(:);
Res.LL_H0  = nan(Ntau,1);
Res.LL_H1  = nan(Ntau,1);
Res.BestPar_H0 = nan(Ntau,numel(BestGuessH0));
Res.BestPar_H1 = nan(Ntau,numel(BestGuessH1));
Res.ExitFlag_H1 = nan(Ntau,1);
Res.ExitFlag_H0 = nan(Ntau,1);

%Options = optimoptions ( 'fminunc','UseParallel',false,'FunctionTolerance',1e-6,'Display','off');



for Itau=1:1:Ntau
    if InPar.Verbose
        fprintf('Fitting time delay %d of %d   -  Tau=%f\n',Itau,Ntau,VecTau(Itau));
    end
    
    Limits = [VecTau(Itau), VecTau(Itau); InPar.Limits];
    
    FitParH1 = [VecTau(Itau), InPar.FitPar(1:end)];
    
    if InPar.TwoD
        FitParH0 = [VecTau(Itau), NaN, NaN, 0,  0 0 0  0 0 0 InPar.FitPar(end)];

        [Res.BestPar_H1(Itau,:),Res.LL_H1(Itau),Res.ExitFlag_H1(Itau)]=InPar.Solver({@TimeDelay.logl_x2d_given_F, ...
                                         FitParH1, Limits, w, Gx_w, Gy_w, F_t, F_w, ...
                                         sigma_F_hat, sigma_x_hat, ...
                                         DFT, DFT_dagger, LogZ, Gamma_1_},...
                                        BestGuessH1,AddPars{:});

            
        [Res.BestPar_H0(Itau,:),Res.LL_H0(Itau),Res.ExitFlag_H0(Itau)]=InPar.Solver({@TimeDelay.logl_x2d_given_F, ...
                                         FitParH0, Limits, w, Gx_w, Gy_w, F_t, F_w, ...
                                         sigma_F_hat, sigma_x_hat, ...
                                         DFT, DFT_dagger, LogZ, Gamma_1_},...
                                        BestGuessH0,AddPars{:});
    else
        % 1-D
        FitParH0 = [VecTau(Itau), NaN, NaN, 0,  0 0 0   InPar.FitPar(end)];
        [Res.BestPar_H1(Itau,:),Res.LL_H1(Itau)]=InPar.Solver({@TimeDelay.logl_x_given_F, ...
                                         FitParH1, Limits, w, Gx_w, F_t, F_w, ...
                                         sigma_F_hat, sigma_x_hat, ...
                                         InPar.Min_w,...
                                         DFT, DFT_dagger, LogZ, Gamma_1_},...
                                        BestGuessH1,AddPars{:});



        [Res.BestPar_H0(Itau,:),Res.LL_H0(Itau)]=InPar.Solver({@TimeDelay.logl_x_given_F, ...
                                         FitParH0, Limits, w, Gx_w, F_t, F_w, ...
                                         sigma_F_hat, sigma_x_hat, ...
                                         InPar.Min_w,...
                                         DFT, DFT_dagger, LogZ, Gamma_1_},...
                                        BestGuessH0,AddPars{:});

       
    end
                                    
    %
    %if InPar.Verbose
    %end
    %Res.FitParH1 = FitParH1;
    %Res.FitParH0 = FitParH0;
    %Res.BestGuessH1 = BestGuessH1;
    %Res.BestGuessH0 = BestGuessH0;
    %Res.Limits = Limits;
    
end







%%



%[LogL_xF,LogL_GF,LogL_F]=logl_x_given_F(Pars, FitPar, Limits, w, G_w, F_t, F_w, sigma_F_hat, sigma_x_hat, DFT, DFT_dagger, LogZ)

%    Tau    = Pars(1);
%    Alpha0 = Pars(2);
%    Alpha1 = Pars(3);
%    Alpha2 = Pars(4);
%    x0     = Pars(5);
%    x1     = Pars(6);
%    x2     = Pars(7);
%    gamma  = Pars(8);
   

%Limits    = [5 100; 0 5; 0 2; 0 2; -1 1; -2.1 2.1; -2.1 2.1; 1.5 3.5];


%%



if 1==0
BestGuess = [2 1 0.66 0 1 -1];
BestGuessH0 = [1];

VecTau = 1./InPar.VecInvTau; 
Ntau   = numel(VecTau);

LL   = nan(Ntau,1);
LLH0 = nan(Ntau,1);
for Itau=1:1:Ntau

    FitPar = [VecTau(Itau)  NaN NaN NaN    NaN NaN NaN 3];

    tic;
    [P,Fval]=Util.fit.fminsearch_my({@TimeDelay.logl_x_given_F, ...
                                     FitPar, Limits, w, Gx_w, F_t, F_w, ...
                                     sigma_F_hat, sigma_x_hat, ...
                                     DFT, DFT_dagger, LogZ, Gamma_1_},...
                                    BestGuess);

    toc;
    
    FitPar = [VecTau(Itau)  2 NaN 0 0   0 0 3];
    tic;
    [PH0,FvalH0]=Util.fit.fminsearch_my({@TimeDelay.logl_x_given_F, ...
                                     FitPar, Limits, w, Gx_w, F_t, F_w, ...
                                     sigma_F_hat, sigma_x_hat, ...
                                     DFT, DFT_dagger, LogZ, Gamma_1_},...
                                    BestGuessH0);

    %
    toc
    
    P
    LL(Itau) = Fval;
    LLH0(Itau) = FvalH0;
    
    
end

plot(VecTau,LL)
hold on;
plot(VecTau,LLH0)

Res.VecTau = VecTau;
Res.LL = LL;
Res.LLH0 = LLH0;

end



%%
if 1==0

    ParsActual   = [14.7 2 1   0.66 0 1 -1 3];

Limits    = [5 100; 0 5; 0 2; 0 2; -1 1; -2.1 2.1; -2.1 2.1; -1 1; -2.1 2.1; -2.1 2.1;  1.5 3.5];


BestGuess = [2 1 0.66 0 1 -1 0 0 0];
BestGuessH0 = [2 1];



VecTau = 14.7; %1./InPar.VecInvTau; 
Ntau   = numel(VecTau);

LL   = nan(Ntau,1);
LLH0 = nan(Ntau,1);
for Itau=1:1:Ntau

    FitPar = [VecTau(Itau)  NaN NaN NaN    NaN NaN NaN  NaN NaN NaN  3];
    FitParH0 = [VecTau(Itau)  NaN NaN 0   0 0 0  0 0 0  3];
    
    Limits(1,1:2) = VecTau(Itau);
    
    tic;
    [P,Fval]=Util.fit.fminsearch_my({@TimeDelay.logl_x2d_given_F, ...
                                     FitPar, Limits, w, Gx_w, Gy_w, F_t, F_w, ...
                                     sigma_F_hat, sigma_x_hat, ...
                                     DFT, DFT_dagger, LogZ, Gamma_1_},...
                                    BestGuess);

    toc;
    
    tic;
    [PH0,FvalH0]=Util.fit.fminsearch_my({@TimeDelay.logl_x2d_given_F, ...
                                     FitParH0, Limits, w, Gx_w, Gy_w, F_t, F_w, ...
                                     sigma_F_hat, sigma_x_hat, ...
                                     DFT, DFT_dagger, LogZ, Gamma_1_},...
                                    BestGuessH0);

    %
    toc
    
    P
    LL(Itau) = Fval;
    LLH0(Itau) = FvalH0;
    

    Res.Tau = VecTau;
    Res.LL_H0 = FvalH0;
    Res.LL_H1 = Fval;
    Res.BestPar_H0 = PH0;
    Res.BestPar_H1 = P;
    Res.FitParH1 = FitPar;
    Res.FitParH0 = FitParH0;
    Res.BestGuessH1 = BestGuess;
    Res.BestGuessH0 = BestGuessH0;
    Res.Limits      = Limits;
    
    
end


end

%%


if 1==0

FitPar = [14.7  2 NaN NaN  1 -1 3];

VecA1    = (0.5:0.05:2).';
VecA2dA1 = (0.3:0.05:1).';
Na1 = numel(VecA1);
Na2 = numel(VecA2dA1);
LogL_xF = nan(Na1,Na2);
for Ia1=1:1:Na1
    A1 = VecA1(Ia1);
    for Ia2=1:1:Na2
        A2 = A1.*VecA2dA1(Ia2);
       
        [LogL_xF(Ia1,Ia2),LogL_GF,LogL_F]=TimeDelay.logl_x_given_F([A1 A2], FitPar, Limits, w, G_w, F_t, F_w, sigma_F_hat, sigma_x_hat, DFT, DFT_dagger, LogZ);

    end
end
    
end

if 1==0
    
FitPar = [14.7  2 1 0.66  NaN NaN 3];
Vecx1    = (-2:0.1:2).';
Vecx2    = (-2:0.1:2).';
Nx1 = numel(Vecx1);
Nx2 = numel(Vecx2);
LogL_xF = nan(Nx1,Nx2);
for Ix1=1:1:Nx1
    x1 = Vecx1(Ix1);
    for Ix2=1:1:Nx2
        x2 = Vecx2(Ix2);
        
        [LogL_xF(Ix1,Ix2),LogL_GF,LogL_F]=TimeDelay.logl_x_given_F([x1 x2], FitPar, Limits, w, G_w, F_t, F_w, sigma_F_hat, sigma_x_hat, DFT, DFT_dagger, LogZ);
        
        
    end
end
    
end


if 1==0
    
FitPar = [14.7  NaN 1 0.66  1 -1 3];
VecA0   = (-2:0.1:5).';
Na0 = numel(VecA0);
LogL_xF = nan(Na0,1);
for Ia0=1:1:Na0
    A0 = VecA0(Ia0);
    
    [LogL_xF(Ia0),LogL_GF,LogL_F]=TimeDelay.logl_x_given_F([A0], FitPar, Limits, w, G_w, F_t, F_w, sigma_F_hat, sigma_x_hat, DFT, DFT_dagger, LogZ);
        
        
end
    
end



