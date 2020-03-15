function LogL=flux_delay_loglH0(Par,PS,ErrF,gamma,Limits)
% Calculate the log-likelihood for a single image flux fit
% Package: TimeDelay
% Description: Calculate the log-likelihood for the null hypotesis of the
%              two images flux time delay fit.
% Input  : - Vector of free paramaeters [A1, [gamma]]
%            gamma is optional.
%          - A two column vector of [Frequency, PowerSpectrum].
%          - A flux error scalar.
%          - gamma. Default is 2.
%          - Limits [min A1, max A1;
%                    min A2/A1, max A2/A1;
%                    min Tau, Max Tau;
%                    min gamma, max gamma].
%            Default is [0 Inf; -eps 1; 5 100; 1.5 3.5].
% Output : - Minus Log likelihhod.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: LogL=TimeDelay.flux_delay_logl(Par,PS,ErrF)
% Reliable: 
%--------------------------------------------------------------------------

ColFreq  = 1;
ColPower = 2;

if nargin<5
    Limits = [0 Inf; -eps 1; 5 100; 1.5 3.5];
    if nargin<4
        gamma = 2;
    end
end


A1    = Par(1);
%A2    = Par(2);
%Tau   = Par(3);

if (numel(Par)>1)
    gamma = Par(2);
end

if (A1<Limits(1,1) || A1>Limits(1,2) || ...
        gamma<Limits(4,1) || gamma>Limits(4,2))
    LogL = Inf;
else

    Omega  = 2.*pi.*PS(:,ColFreq);
    SigmaF = abs(Omega).^(-gamma).*A1.^2  + ErrF.^2;


    % log determinant of a diagonal matrix
    LogL = -0.5.*sum(log(2.*pi.*SigmaF)) - sum(PS(:,ColPower)./(2.*SigmaF));

    % return the -LogL
    LogL = -LogL;
end