function LogL=ast_delay_logl(InPar,Data,ErrF,FitFlag,Limits,WinPS)
% Calculate the log-likelihood for two images flux time delay fit
% Package: TimeDelay
% Description: Calculate the log-likelihood for two images flux time delay
%              fit.
% Input  : - Vector of free paramaeters - e.g.,
%            [gamma, A0, A1 x1, y1, Tau2, A2, x2, y2,...]
%            In order to fit some of the papameters use the FitFlag input
%            parameter.
%          - A structure containing the observational data:
%            The following fields should be provided:
%            .Time  - Vector of times.
%            .F_t   - Vector of total flux
%            .X_t   - Vector of center of mass X position.
%            .Y_t   - Vector of center of mass Y position.
%            .Freq  - Vector of frequencies
%            .Omega - Vector of 2.*pi.*Frequency
%            .F_w   - power spectrum at frequencies (freq 0 is at the first
%            position).
%          - A two column vector of [Frequency, PowerSpectrum].
%          - A three element vector [sigma_F, sigma_x, sigma_y].
%          - Vector of deafault values parameters which will not be fitted
%            by the order 
%            [gamma, A0, A1 x1, y1, Tau2, A2, x2, y2,...]
%            If NaN, then parameter will be fitted.
%            Default is [NaN NaN NaN NaN NaN NaN NaN NaN NaN].
%          - two column matrix of limits [min, max],
%            for each of the parameters.
%            Default is [2 3;0 Inf;0 Inf; -2 2;-2 2; 5 100; 0 Inf; -2 2; -2 2].
%          - Optional window spectrum vector for the same frequencies.
%            If empty, then don't use the window function.
%            Default is []. 
% Output : - Minus Log likelihhod.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Par = [3 2 1 1 0 14.7 0.66 -1 0];
%          Data.F_t = Util.IO.load2('F_t.txt');
%          Data.F_w = fft(Data.F_t);
%          Data.X_t = Util.IO.load2('x_t.txt');
%          Nt = numel(Data.F_t);
%          Data.Y_t = zeros(Nt,1);
%          Data.Omega = 2.*pi.*ifftshift(TimeDelay.fft_freq(Nt))./Nt;
%          ErrF = [0.175, 0.1, 0.1]; % in flux/position units
%          LogL=TimeDelay.ast_delay_logl(Par,Data,Err)
% Reliable: 
%--------------------------------------------------------------------------

ColFreq  = 1;
ColPower = 2;

DefGamma = 2;

if nargin<6
    WinPS = [];
    if nargin<5
        Limits = [1.5 3;0 Inf;0 Inf; -2 2;-2 2; 5 100; 0 Inf; -2 2; -2 2];
        if nargin<4
            FitFlag = [NaN NaN NaN NaN NaN NaN NaN NaN NaN];
        end
    end
end

FlagNaN = isnan(FitFlag);
Par     = FitFlag;
Par(FlagNaN) = InPar;



Npar = numel(InPar);
if sum(isnan(FitFlag))~=Npar
    error('Numner of parameters to fit is not consistent with FitFlag input');
end



Skip = false;
for Ipar=1:1:Npar
    if Par(Ipar)<Limits(Ipar,1) || Par(Ipar)>Limits(Ipar,2)
        % parameter is out of bound
        
        Skip = true;
    end
end

%     if isnan(FitFlag)
%         % check if within limits
%         if Par(Ipar)<Limits(Ipar,1) || Par(Ipar)>Limits(Ipar,2)
%             % parameter is out of bound
%             LogL = Inf;
%         end
%     else
%         % do not fit parameter
%         % use value in FitFlag
%         Par(Ipar) = FitFlag(Ipar);
%     end
% end


if Skip
    LogL = Inf;
else
    %[gamma, A0, A1 x1, y1, Tau2, A2, x2, y2,...]
        
    gamma = Par(1);
    A0    = Par(2);
    A     = Par(3:4:end);
    Xi    = Par(4:4:end);
    Yi    = Par(5:4:end);
    Tau   = [0, Par(6:4:end)];
    
    
    sigma_F_t = ErrF(1);
    sigma_X_t = ErrF(2);
    sigma_Y_t = ErrF(3);
    

    F_t   = Data.F_t;
    F_w   = Data.F_w;
    X_t   = Data.X_t;
    Y_t   = Data.Y_t;
    Omega = Data.Omega;

    Nim    = numel(A);
    Nt     = numel(F_t);
    Nw     = numel(Omega);
    
    sigma_F_w = sigma_F_t.*sqrt(Nt);
    sigma_X_w = sigma_X_t.*sqrt(Nt);
    sigma_Y_w = sigma_Y_t.*sqrt(Nt);



    % input:
    % X_t - x(t)
    % Y_t - y(t)
    % F_t - F(t)
    % Nt  - number of observations (N)
    % sigma_X_t
    % sigma_Y_t
    % sigma_F_t


    % Alpha  = [1 0.66];
    % Alpha0 = 2  % (at position 0,0)
    % Xi     = [1 -1];
    % Yi     = [0 0];
    % Xi0    = 0;
    % Yi0    = 0;
    % Tau    = [0 14.7];
    % gamma  = 3;
    % 
    % Nim    = numel(Alpha);
    % 
    % %FT = @(t,f,Freq) sum(f(:).*exp(-2.*pi.*1i.*Freq(:).'.*t(:)),1).'; 
    %  
    % 
    % %TimeVec = (1:1:100).';
    % F_t     = rand(size(TimeVec));
    % X_t     = rand(size(TimeVec));
    % Y_t     = rand(size(TimeVec));




    %    sigma_X_t = 0.1;
    %    sigma_Y_t = 0.00;
    % 
    %    Y_t = zeros(size(X_t));
    %    
    % 
    % DeltaFreq = 1./range(TimeVec);
    % MaxFreq   = 0.5./median(diff(TimeVec));
    % FreqVec  = (0:DeltaFreq:MaxFreq).';
    % FreqVec  = [-flipud(FreqVec(2:end)); FreqVec];


    % Fouerir transforms:
    %Nt  = numel(F_t);
    %FreqVec = ifftshift(TimeDelay.fft_freq(Nt))./Nt;
    %Omega = 2.*pi.*FreqVec;
    %Nw  = numel(Omega);

    [Fjl] = TimeDelay.dft_matrix(Nt);

    %F_w = fft(F_t); % FT(TimeVec, F_t, FreqVec);
    %X_w = fft(X_t); % FT(TimeVec, X_t, FreqVec);
    %Y_w = fft(Y_t); % FT(TimeVec, Y_t, FreqVec);

    Gx_t = X_t.*F_t;
    Gy_t = Y_t.*F_t;
    Gx_w = fft(Gx_t);
    Gy_w = fft(Gy_t);

    %sigma_F_t = 0.175;   % sigma_F_prop.*fDC

    % Z = prod(abs(F_t));
    LogZ = sum(log(abs(F_t)));

    % for two images
    % Sigma_phi_w = (Alpha(1).^2 + Alpha(2).^2 + 2.*Alpha(1).*Alpha(2).*cos(Omega.*Tau(2))).*abs(Omega).^(-gamma);

    % for N images
    SumAT_w = zeros(Nw,1);
    for Iw=1:1:Nw
        SumAT_w(Iw) = 2.*sum(triu(A(:).'.*A(:).*cos(Omega(Iw).*(Tau(:).' - Tau(:))),1),'all');
    end
    Sigma_phi_w = (sum(A.^2) + SumAT_w).*abs(Omega).^(-gamma);
    Sigma_F_w   = diag(Sigma_phi_w) + sigma_F_w.^2; % keep this is a vector (representing a diagonal matrix)
    % convert into a diagonal matrix:
    Sigma_phi_w = diag(Sigma_phi_w);


    % for two images
    % Ax_w         = (Alpha(1).*Xi(1) + Alpha(2).*Xi(2).*exp(1i.*Omega.*Tau))./(Alpha(1) + Alpha(2).*exp(1i.*Omega.*Tau));

    % for N images
    AlphaTau_w = sum(A(:).*exp(1i.*Omega(:).'.*Tau(:)),1);

    % need to be a column vector!
    Ax_w = [sum(A(:).*Xi(:).*exp(1i.*Omega(:).'.*Tau(:)),1)./ AlphaTau_w].';
    Ay_w = [sum(A(:).*Yi(:).*exp(1i.*Omega(:).'.*Tau(:)),1)./ AlphaTau_w].';
    % convert into diagonal matrix
    %Ax_w = diag(Ax_w);
    %Ay_w = diag(Ay_w);

    % ???
    %Ax_w = realAx_w);
    %Ay_w = real(Ay_w);


    % 1. Ax_w*F_w,... should be a diagonal matrix
    % 2. number of elements in F_t and F_w must be identical
    % 3. The times in F_t should corresponds to the time in ifft(Ax_w.*F_w)...

    % F_wDCS - is DC subtrcated version of F_w (0 at frequency=0).
    F_wDCS = F_w;
    F_wDCS(Omega==0) = F_wDCS(Omega==0) - A0.*Nt;



    % original - doesn't work - need to subtracted DC
    Xx_t = ifft(Ax_w.*F_w)./F_t;
    Xy_t = ifft(Ay_w.*F_w)./F_t;

    Xx_t = ifft(Ax_w.*F_wDCS)./F_t;
    Xy_t = ifft(Ay_w.*F_wDCS)./F_t;


    % apply the dft transformation:
    Xx_w = Fjl*diag(Xx_t)*Fjl'  ./Nt./Nt;
    Xy_w = Fjl*diag(Xy_t)*Fjl'  ./Nt./Nt;

    % ???
    %Xx_w = real(Xx_w);
    %Xy_w = real(Xy_w);



    % \Gamma_{\hat{\epsilon}_{F}|\hat{F}}


    Gamma_eFgF = (sigma_F_w.^(-2) + 1./Sigma_phi_w).^(-1);    
    % \mu_{\hat{\epsilon}_{F}|\hat{F}}
    Mu_eFgF    = F_w.*sigma_F_w.^2./(sigma_F_w.^2 + diag(Sigma_phi_w));   

    % \Gamma_{\hat{G}|\hat{F}}
    % problem is in: sigma_X_w.^2.*F_w.*F_w'...
    % Ofer's calc:
    Gamma_x_GgF  = sigma_X_w.^2.* Fjl*diag(F_t.^2)*Fjl'./Nt + (Xx_w - diag(Ax_w))*Gamma_eFgF*(Xx_w - diag(Ax_w))';

    % my calc:
    Gamma_x_GgF  = sigma_X_w.^2.*F_w.*F_w' + (Xx_w - diag(Ax_w))*Gamma_eFgF*(Xx_w - diag(Ax_w))';
    Gamma_y_GgF  = sigma_Y_w.^2.*F_w.*F_w' + (Xy_w - diag(Ay_w))*Gamma_eFgF*(Xy_w - diag(Ay_w))';

    % \mu_{\hat{G}|\hat{F}}
    Mu_x_GgF     = Ax_w.*F_w + (Xx_w - Ax_w)*Mu_eFgF;   %<---- why isn't it a vector?
    Mu_y_GgF     = Ay_w.*F_w + (Xy_w - Ay_w)*Mu_eFgF;

    % calculate the log-likelihood
    % but remove zero frequency
    Flag = abs(Omega)>0;

    LogP_F = -0.5.*sum(log(2.*pi.*Sigma_F_w(Flag))) - sum(abs(F_w(Flag)).^2./(2.*Sigma_F_w(Flag)));


    % 1. do we need to remove zero frequency - no?!
    % 2. the answer is complex - why? beacuse ifft not symmetric?
    LogP_x = -TimeDelay.logdet(pi.*Gamma_x_GgF) - (Gx_w - Mu_x_GgF)' * inv(Gamma_x_GgF) * (Gx_w - Mu_x_GgF);
    % Consider using pinv to avoid division by zero:
    LogP_y = -TimeDelay.logdet(pi.*Gamma_y_GgF) - (Gy_w - Mu_y_GgF)' * inv(Gamma_y_GgF) * (Gy_w - Mu_y_GgF);
    if isnan(LogP_x) || isinf(LogP_x)
        LogP_x = 0;
    end
    if isnan(LogP_y) || isinf(LogP_y)
        LogP_y = 0;
    end


    LogL = LogP_F + real(LogP_x) + real(LogP_y);
    LogL = -LogL;

end





