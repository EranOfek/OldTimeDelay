function LogL=flux_delay_logl(Par,PS,ErrF,FitFlag,Limits,WinPS)
% Calculate the log-likelihood for two images flux time delay fit
% Package: TimeDelay
% Description: Calculate the log-likelihood for two images flux time delay
%              fit.
% Input  : - Vector of free paramaeters - e.g.,
%            [A1, A2, Tau, gamma], or [A2, Tau], or [A1, A2, gamma]
%            In order to fit some of the papameters use the FitFlag input
%            parameter.
%          - A two column vector of [Frequency, PowerSpectrum].
%          - A flux error scalar.
%          - Vector of deafault values parameters which will not be fitted
%            by the order [A1, A2, Tau, gamma].
%            If NaN, then parameter will be fitted.
%            Default is [NaN NaN NaN 2].
%          - Limits [min A1, max A1;
%                    min A2/A1, max A2/A1;
%                    min Tau, Max Tau;
%                    min gamma, max gamma].
%            Default is [0 Inf; -eps 1; 5 100; 1.5 3.5].
%          - Optional window spectrum vector for the same frequencies.
%            If empty, then don't use the window function.
%            Default is []. 
% Output : - Minus Log likelihhod.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: LogL=TimeDelay.flux_delay_logl(Par,PS,ErrF)
% Reliable: 
%--------------------------------------------------------------------------

ColFreq  = 1;
ColPower = 2;

DefGamma = 2;

if nargin<6
    WinPS = [];
    if nargin<5
        Limits = [0 Inf; -eps 1; 5 100; 1.5 3.5];
        if nargin<4
            FitFlag = [NaN NaN NaN DefGamma];
        end
    end
end

Npar = numel(Par);
if sum(isnan(FitFlag))~=Npar
    error('Numner of parameters to fit is not consistent with FitFlag input');
end

IndPar = cumsum(isnan([FitFlag]));
% select parameters to fit
if isnan(FitFlag(1))
    A1 = Par(IndPar(1));
else
    A1 = FitFlag(1);
end
if isnan(FitFlag(2))
    A2 = Par(IndPar(2));
else
    A2 = FitFlag(2);
end
if isnan(FitFlag(3))
    Tau = Par(IndPar(3));
else
    Tau = FitFlag(3);
end
if isnan(FitFlag(4))
    gamma = Par(IndPar(4));
else
    gamma = FitFlag(4);
end
%[A1,A2, Tau,gamma]

if (A1<Limits(1,1) || A1>Limits(1,2) || ...
            (A2./A1)<Limits(2,1) || (A2./A1)>Limits(2,2) || ...
            Tau<Limits(3,1) || Tau>Limits(3,2) || ...
            gamma<Limits(4,1) || gamma>Limits(4,2))
    LogL = Inf;
else
    Omega  = 2.*pi.*PS(:,ColFreq);
    SigmaF = abs(Omega).^(-gamma).*( A1.^2 + A2.^2 + 2.*A1.*A2.*cos(Omega.*Tau) ) + ErrF.^2;

    % convolve with window function (doesn't work)
    if ~isempty(WinPS)
        SigmaF = conv(SigmaF,WinPS,'same');
    end
    

    % log determinant of a diagonal matrix
    LogL = -0.5.*sum(log(2.*pi.*SigmaF)) - sum(PS(:,ColPower)./(2.*SigmaF));

    % return the -LogL
    LogL = -LogL;
    
    
end

% if (numel(Par)<3)
%     % H0
%     % only A1, [gamma] are fitted
%     
%     A1    = Par(1);
%     A2    = 0;
%     Tau   = 0;
%     
%     if (numel(Par)>1)
%         gamma = Par(2);
%     else
%         gamma = DefGamma;
%     end
% 
%     if (A1<Limits(1,1) || A1>Limits(1,2) || ...
%             (A2./A1)<Limits(2,1) || (A2./A1)>Limits(2,2) || ...
%             Tau<Limits(3,1) || Tau>Limits(3,2) || ...
%             gamma<Limits(4,1) || gamma>Limits(4,2))
%         LogL = Inf;
%     else
%         LogL = [];
%     end
%     
% else
%     % H1
%     % [A1, A2, Tau [gamma]] are fitted
%     
%     A1    = Par(1);
%     A2    = Par(2);
%     Tau   = Par(3);
% 
%     if (numel(Par)>3)
%         gamma = Par(4);
%     else
%         gamma = DefGamma;
%     end
%     if (A1<Limits(1,1) || A1>Limits(1,2) || ...
%             gamma<Limits(4,1) || gamma>Limits(4,2))
%         LogL = Inf;
%     else
%         LogL = [];
%     end
%     
% end
% 
%     
% if isempty(LogL)
%     Omega  = 2.*pi.*PS(:,ColFreq);
%     SigmaF = abs(Omega).^(-gamma).*( A1.^2 + A2.^2 + 2.*A1.*A2.*cos(Omega.*Tau) ) + ErrF.^2;
% 
% 
%     % log determinant of a diagonal matrix
%     LogL = -0.5.*sum(log(2.*pi.*SigmaF)) - sum(PS(:,ColPower)./(2.*SigmaF));
% 
%     % return the -LogL
%     LogL = -LogL;
% end