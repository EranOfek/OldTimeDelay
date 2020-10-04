function [LogLike,]=logLike_x_given_F(Par,ParFixed,Limits,F_t,x_t,X_i,sigma_F,sigma_x)
%

% Par:
% [A0 A1 A2 Tau gamma]

sigma_F_hat = sigma_F;
sigma_x_hat = sigma_x;

if isempty(FitFlag)
    %         [A0  A1  A2  Tau gamma]
    FitFlag = [NaN NaN NaN NaN NaN];
end

IndPar = cumsum(isnan([ParFixed]));

% select parameters to fit
ParInd = 1;
if isnan(ParFixed(ParInd))
    Alpha_0 = Par(IndPar(ParInd));
else
    Alpha_0 = ParFixed(ParInd);
end
Alpha = zeros(1,2);
ParInd = 2;
if isnan(ParFixed(ParInd))
    Alpha(1) = Par(IndPar(ParInd));
else
    Alpha(1) = ParFixed(ParInd);
end
ParInd = 3;
if isnan(ParFixed(ParInd))
    Alpha(2) = Par(IndPar(ParInd));
else
    Alpha(2) = ParFixed(ParInd);
end
ParInd = 4;
if isnan(ParFixed(ParInd))
    Tau(1) = Par(IndPar(ParInd));
else
    Tau(1) = ParFixed(ParInd);
end
ParInd = 5;
if isnan(ParFixed(ParInd))
    gamma = Par(IndPar(ParInd));
else
    gamma = ParFixed(ParInd);
end


%%
N_t  = numel(F_t);
D_t  = F_t(2) - F_t(1);

F_w  = fft(F_t) ./ sqrt(N_t);
Freq = TimeDelay.fft_freq(N_t, D_t);

w          = 2.*pi.*Freq;
DFT        = fft(eye(N_t), N_t, 1) / sqrt(N_t);
DFT_dagger = DFT';

% Sigma_phi
w_ = w;
w_(1) = 1.0;
Sigma_phi_ = (Alpha(1).^2 + Alpha(2).^2 + 2.*Alpha(1).*Alpha(2).*cos(w.*Tau)) ./ (abs(w_).^gamma);
Sigma_F_ = Sigma_phi_ + sigma_F_hat^2;

% mu_epsilon_F_given_F
Shrink = sigma_F_hat.^2 ./ (sigma_F_hat.^2 + Sigma_phi_);
mu_epsilon_F_given_F_ = Shrink .* F_w;





%x_of_F_w_ = x_of_F_w(F_w, tau, alpha_1, alpha_2);
%n = length(w);
F_w_tilde = F_w;
F_w_tilde(1) = F_w_tilde(1) - Alpha_0.*N_t;  %length(w);
Diag_A_hat   = A_hat(Tau, Alpha, X_i, w);
F_t_ = real(ifft(F_w)) * sqrt(N_t);
x_of_F_w = (real(ifft(Diag_A_hat.*F_w_tilde)) ./ F_t_) * sqrt(N_y);

% X_hat_of_F_w
X_t = diag(x_of_F_w(F_w, tau, alpha_1, alpha_2));
    X_hat_ = ((DFT * X_t) * DFT_dagger) / N_t;  % TODO: check factor len(w)**-1


G_t = x_t.*F_t;
G_w = fft(G_t) ./ sqrt(N_t);

%mu_G_given_F_ = mu_G_given_F(F_w, tau, alpha_1, alpha_2);

   %mu_epsilon_F_given_F_ = mu_epsilon_F_given_F(F_w, tau, alpha_1, alpha_2);
   
    %X_hat_ = X_hat_of_F_w(F_w, tau, alpha_1, alpha_2);
    
    
    mu_ = transpose(Diag_A_hat.*F_w) + (X_hat_ - diag(Diag_A_hat)) * transpose(mu_epsilon_F_given_F_);
    mu_G_given_F_ = transpose(mu_);
    
%A_hat_F_w_ = A_hat(tau, alpha_1, alpha_2) .* F_w;
A_hat_F_w_ = Diag_A_hat.*F_w;


%[Gamma_G_given_F_, Gamma_1_, Gamma_2_] = Gamma_G_given_F(F_t, F_w, tau, alpha_1, alpha_2);
%X_hat_ = X_hat_of_F_w(F_w, tau, alpha_1, alpha_2);
%A_hat_ = A_hat(tau, alpha_1, alpha_2);
X_minus_A = X_hat_ - diag(Diag_A_hat);
    
    %Gamma_epsilon_F_given_F_ = Gamma_epsilon_F_given_F(tau, alpha_1, alpha_2);
    Gamma_epsilon_F_given_F_ = real(diag((sigma_F_hat.^-2 + Sigma_phi_.^-1).^-1));
    
    Gamma_1_ = ((DFT * diag(F_t.^2)) * DFT_dagger) * sigma_x_hat^2;
    Gamma_2_ = (X_minus_A * Gamma_epsilon_F_given_F_) * X_minus_A';
    Gamma_G_given_F_ = Gamma_1_ + Gamma_2_;

    
LogZ = sum(log(F_t));

    
% thing that depands on the specific model
    
    
for i=1:length(tau_)
    display(i);
    for j=1:length(alpha_2_)
        ll_F(i, j) = log_like_F(F_w, tau_(i), alpha_1, alpha_2_(j));
        
        
        
    Sigma_F_ = Sigma_F(tau, alpha_1, alpha_2);
    Sigma_F_ = Sigma_F_(2:end);
    
    power_F_ = (abs(F_w).^2); power_F_ = power_F_(2:end);
    ll = -0.5*sum(power_F_ ./ Sigma_F_) - 0.5*sum(log(2*pi*Sigma_F_));
        ll_G_given_F(i, j) = log_like_G_given_F(G_w, F_t, F_w, tau_(i), alpha_1, alpha_2_(j));
        
    end
end

ll_x_given_F = ll_F + (-1./LogZ).*ll_G_given_F;
    
    