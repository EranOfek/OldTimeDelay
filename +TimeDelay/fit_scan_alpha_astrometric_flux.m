function [Res]=fit_scan_alpha_astrometric_flux(t,F_t,x_t,y_t,sigma_F,sigma_x,varargin)
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
%            'Verbose' - Default is true.
% Output : - An output structure containing the following fields:
%            .Tau 
%            .LL_H0
%            .LL_H1
%            .BestPar_H0
%            .BestPar_H1
% Example: Res=TimeDelay.fit_scan_alpha_astrometric_flux
%          ResF=TimeDelay.timedelayed_lc;
% Res=TimeDelay.fit_scan_alpha_astrometric_flux(ResF.T, ResF.F_t,ResF.x_t,ResF.y_t,ResF.eps_F_abs,ResF.eps_x_abs);


NPAR2D = 11 -1;
NPAR1D = 8 -1;

InPar = inputParser;
addOptional(InPar,'Solver',@Util.fit.fminunc_my);
addOptional(InPar,'FitPar',[0 1   0.66  0   1   -1   0   0   0    3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
addOptional(InPar,'Limits',[0 20;0 10;0 10;  -1 1; -2.1 2.1; -2.1 2.1;  -1 1; -2.1 2.1; -2.1 2.1;   1.5 3.5]); %  without Tau
addOptional(InPar,'VecA1',(0.5:0.05:1.5));
addOptional(InPar,'VecA2dA1',(0.4:0.05:1));

addOptional(InPar,'Tau',14.7); 
addOptional(InPar,'Verbose',true);
parse(InPar,varargin{:});
InPar = InPar.Results;


if nargin==0
    params = jsondecode(fileread('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/sims_4.json'));

    t   = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/t.txt')';
    F_t = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/F_t.txt')';
    x_t = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/x_t.txt')';
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
    
    F_w = fft(F_t) / sqrt(N);
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
% if InPar.TwoD
%     if ~(numel(InPar.FitPar)==NPAR2D && numel(InPar.DefPar)==NPAR2D && size(InPar.Limits,1)==NPAR2D)
%         error('For 2D fitting Limits, FitPar and DefPar must contain %d elements',NPAR2D);
%     end
% else
%     if ~(numel(InPar.FitPar)==NPAR1D && numel(InPar.DefPar)==NPAR1D && size(InPar.Limits,1)==NPAR1D)
%         error('For 1D fitting Limits, FitPar and DefPar must contain %d elements',NPAR1D);
%     end
% end


    
FlagN       = isnan(InPar.FitPar);


Limits = [InPar.Tau, InPar.Tau; InPar.Limits];
FitParsH1 = [InPar.Tau, InPar.FitPar(1:end)];
FitParsH1(3:4) = NaN;

Na1 = numel(InPar.VecA1);
Na2 = numel(InPar.VecA2dA1);
%Res.LL_H0  = nan(Ntau,1);
Res.A1     = InPar.VecA1;
Res.A2dA1  = InPar.VecA2dA1;
Res.LL_H1  = nan(Na1,Na2);

Na1.*Na2

for Ia1=1:1:Na1
    for Ia2=1:1:Na2
        
        A1 = InPar.VecA1(Ia1);
        A2 = A1.*InPar.VecA2dA1(Ia2);
        
        [LogL_xyF,LogLx_GF,LogLy_GF,LogL_F] = TimeDelay.logl_x2d_given_F([A1 A2], FitParsH1, Limits, w, Gx_w, Gy_w, F_t, F_w,...
            sigma_F_hat, sigma_x_hat, DFT, DFT_dagger, LogZ, Gamma_1_);
        Res.LL_H1(Ia1,Ia2) = LogL_xyF;
        
    end
end

