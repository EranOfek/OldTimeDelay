function [F_w,LC]=gradient_descent(Fun,Y,ParInit,FitFlag,AddPar)
% 
% Package: +TimeDelay
% Description: 
% Input  : - Function handle. The function is of the form:
%            Chi2=Fun(Par,Y,FitFlag,AddPar{:})
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N=240;
%          Omega=ifftshift(TimeDelay.fft_freq(N))./N .*2.*pi;
%          F_w=TimeDelay.rand_psd(Omega.^-2)
% Reliable: 
%--------------------------------------------------------------------------

% Tau is a vector of all Tau to scan simultanosuly
% Y is the PS
% Omega is the corresponding frequency
% A1, A2, gamma, ErrF

dA1    = 1e-4;
dA2    = 1e-4;
dgamma = 1e-4;

Tau = Tau(:).';

Converge = false;
Step = ones(size(Tau)).*1e-4;
A    = [A1;A2].*ones(size(Tau));
while ~Converge
    
    LogL     = logLikelihood(Y,Omega,ErrF,Tau,A(1,:),             A(2,:),             gamma);
    LogL_dA1 = logLikelihood(Y,Omega,ErrF,Tau,A(1,:)+dA1.*A(1,:), A(2,:),             gamma);
    LogL_dA2 = logLikelihood(Y,Omega,ErrF,Tau,A(1,:),             A(2,:)+dA2.*A(2,:), gamma);
    Grad = [LogL-LogL_dA1; LogL-LogL_dA2];
    Anew = A - Step.*Grad;
    
    LogLn     = logLikelihood(Y,Omega,ErrF,Tau,Anew(1,:),             Anew(2,:),             gamma);
    LogLn_dA1 = logLikelihood(Y,Omega,ErrF,Tau,Anew(1,:)+dA1.*A(1,:), Anew(2,:),             gamma);
    LogLn_dA2 = logLikelihood(Y,Omega,ErrF,Tau,Anew(1,:),             Anew(2,:)+dA2.*A(2,:), gamma);
    GradN = [LogLn-LogLn_dA1; LogLn-LogLn_dA2];
    Anew = A - Step.*Grad;
    
    
    Step = (Anew - A).'.* (Grad - GradOld  )
   
    
    
end

end

function LogL=logLikelihood(Y,Omega,ErrF,Tau,A1,A2,gamma)
    %
    SigmaF = abs(Omega).^(-gamma).*( A1.^2 + A2.^2 + 2.*A1.*A2.*cos(Omega.*Tau) ) + ErrF.^2;
    % log determinant of a diagonal matrix
    LogL = -0.5.*sum(log(2.*pi.*SigmaF)) - sum(Y./(2.*SigmaF));
    % return the -LogL
    LogL = -LogL;  % LogL vector for all Tau, and specific A1, A2, gamma
end