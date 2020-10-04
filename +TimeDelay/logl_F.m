function [LogLikeF,Sigma_F,Sigma_phi]=logl_F(w,F_w,Tau,Alpha1,Alpha2,gamma,sigma_F_hat)
% Calculate the log-likelihood of F, Sigma_F, and Sigma_phi
% Package: +TimeDelay
% Input  : - Vector of angular frequencies (as generated by fft),
%            so 0 frequency is at position 1.
%          - fft(F_t)./sqrt(N_t);
%          - Time delay.
%          - Alpha1
%          - Alpha2
%          - gamma
%          - sigma_F_hat
% Output : - An array of Log likelihood of F for each
%            Tau,Alpha1,Alpha2,gamma
%          - Sigma_F, column per parameter Tau/Alpha/gamma
%          - Sigma_phi, column per parameter Tau/Alpha/gamma
% Example: [LogLikeF,Sigma_F,Sigma_phi]=TimeDelay.logl_F(w,F_w,14,1,0.65,2,0.01)

w_     = w;
w_(1)  = 1.0;
% frequency is column vector
w_     = w_(:);
% free parameters should be row vector
Tau    = Tau(:).';
Alpha1 = Alpha1(:).';
Alpha2 = Alpha2(:).';
gamma  = gamma(:).';

% sigma_phi is a matrix of (frequrncy/columns, parameters/row):
Sigma_phi = (Alpha1.^2 + Alpha2.^2 + 2.*Alpha1.*Alpha2.*cos(w.*Tau)) .* (abs(w_).^(-gamma));
Sigma_F   = Sigma_phi + sigma_F_hat.^2;

Power_F = (abs(F_w).^2);
 
LogLikeF = -0.5.*sum(Power_F(2:end) ./ Sigma_F(2:end,:)) - 0.5.*sum(log(2.*pi.*Sigma_F(2:end,:)));
