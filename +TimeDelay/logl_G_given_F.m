function [LogL_GF,LogL_F]=logl_G_given_F(w, G_w, F_t, F_w, Tau, Alpha0, Alpha1, Alpha2,x1,x2,gamma,sigma_F_hat,DFT,DFT_dagger,LogZ)
%
% Example: [LogL_GF,LogL_F]=TimeDelay.logl_G_given_F

if nargin==0
    params = jsondecode(fileread('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/sims_4.json'));

    t   = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/t.txt')';
    F_t = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/F_t.txt')';
    x_t = load('/home/eran/matlab/TimeDelay/Ofer/qtd_sim_4/output_txt/sims_4/x_t.txt')';

    sigma_x = params.sigma_x_prop * abs(params.x_1 - params.x_2);
    sigma_F = params.sigma_F_prop * params.f_dc;

    N = length(t);
    freqs = TimeDelay.fft_freq(N, params.t_step);

    gamma = params.gamma;
    sigma_F_hat = sigma_F;
    sigma_x_hat = sigma_x;
    w = 2.*pi*freqs;
    DFT = fft(eye(N), N, 1) / sqrt(N);
    DFT_dagger = DFT';
    x1 = params.x_1;
    x2 = params.x_2;
    Alpha0 = params.alpha_0;

    Tau = params.tau;
    Alpha1 = params.alpha_1;
    Alpha2 = params.alpha_2;
    F_w = fft(F_t) / sqrt(N);

    G_t = x_t.*F_t;
    G_w = fft(G_t) ./ sqrt(N);

    LogZ = sum(log(F_t));
end

Tau = 18
Alpha1 = 2
Alpha2 = 0.7



%mu_G_given_F_ = mu_G_given_F(F_w, tau, alpha_1, alpha_2);


%mu_epsilon_F_given_F_ = mu_epsilon_F_given_F(F_w, tau, alpha_1, alpha_2);

%Sigma_phi_ = Sigma_phi(tau, alpha_1, alpha_2);

N = numel(w);

if nargin<13
    DFT = fft(eye(N), N, 1) ./ sqrt(N);
    if nargin<14
        DFT_dagger = DFT';
        if nargin<15
            LogZ = sum(log(F_t));
        end
    end
end


% parameters should be row vectors
% Tau    = Tau(:).';
% Alpha1 = Alpha1(:).';
% Alpha2 = Alpha2(:).';
% gamma  = gamma(:).';

% w should be a column vector
G_w = G_w(:).';
F_t = F_t(:).';
F_w = F_w(:).';

w_ = w;
w_(1) = 1.0;


Sigma_phi_ = (Alpha1.^2 + Alpha2.^2 + 2.*Alpha1.*Alpha2.*cos(w.*Tau)) .* (abs(w_).^(-gamma));
Sigma_F    = Sigma_phi_ + sigma_F_hat.^2;
Power_F    = (abs(F_w).^2);
LogL_F     = -0.5.*sum(Power_F(2:end) ./ Sigma_F(:,2:end)) - 0.5.*sum(log(2.*pi.*Sigma_F(:,2:end)));

Shrink = sigma_F_hat.^2 ./ (sigma_F_hat.^2 + Sigma_phi_);
mu_epsilon_F_given_F_ = Shrink .* F_w;

%X_hat_ = X_hat_of_F_w(F_w, tau, alpha_1, alpha_2);
%x_t_ = x_of_F_w(F_w, tau, alpha_1, alpha_2));

F_w_tilde = F_w;
F_w_tilde(1) = F_w_tilde(1) - Alpha0.*N;

%A_hat_ = A_hat(tau, alpha_1, alpha_2);
num_ = Alpha1.*x1 + Alpha2.*x2.*exp(1j.*w.*Tau);
denom_ = Alpha1 + Alpha2.*exp(1j.*w.*Tau);
A_hat_ = num_ ./ denom_;

F_t_ = real(ifft(F_w)) .* sqrt(N);
x_t_ = (real(ifft(A_hat_.*F_w_tilde)) ./ F_t_) .* sqrt(N);

X_t = diag(x_t_);
X_hat_ = ((DFT * X_t) * DFT_dagger) ./N;  % TODO: check factor len(w)**-1

%A_hat_ = A_hat(tau, alpha_1, alpha_2);


%mu_ = transpose(A_hat_.*F_w) + (X_hat_ - diag(A_hat_)) * transpose(mu_epsilon_F_given_F_);
mu_ = transpose(A_hat_.*F_w) + (X_hat_ - diag(A_hat_)) * transpose(mu_epsilon_F_given_F_);

mu_G_given_F_ = transpose(mu_);

G_minus_mu = G_w - mu_G_given_F_; 
G_minus_mu = G_minus_mu(2:end);
%[Gamma_G_given_F_, ~, ~] = Gamma_G_given_F(F_t, F_w, tau, alpha_1, alpha_2);
%X_hat_ = X_hat_of_F_w(F_w, tau, alpha_1, alpha_2);

%A_hat_ = A_hat(tau, alpha_1, alpha_2);
X_minus_A = X_hat_ - diag(A_hat_);
%Gamma_epsilon_F_given_F_ = Gamma_epsilon_F_given_F(tau, alpha_1, alpha_2);
%Sigma_phi_ = Sigma_phi(tau, alpha_1, alpha_2);
Gamma_epsilon_F_given_F_ = real(diag((sigma_F_hat.^-2 + Sigma_phi_.^-1).^-1));


Gamma_1_ = ((DFT * diag(F_t.^2)) * DFT_dagger) * sigma_x_hat^2;
Gamma_2_ = (X_minus_A * Gamma_epsilon_F_given_F_) * X_minus_A';
Gamma_G_given_F_ = Gamma_1_ + Gamma_2_;


Gamma_G_given_F_ = Gamma_G_given_F_(2:end, 2:end);
inv_Gamma_G_given_F = inv(Gamma_G_given_F_);
logdet_ = TimeDelay.logdet(pi.*Gamma_G_given_F_);
mahal = conj(G_minus_mu) * (inv_Gamma_G_given_F * transpose(G_minus_mu));
LogL_GF = -real(logdet_) - real(mahal);

