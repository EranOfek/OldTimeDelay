function Fit=flux_delay_fit(varargin)
% Fit a time delay to a power spectrum of a process with time delay
% Package: +TimeDelay
% Description: Given a power spectrum of a time series, fit a time delay,
%              flux ratio, amplitude, and power-law index of the power
%              spectrum.
% Input  : * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            
% Output : - A structure containing the following fields
%            
% Dependency: Requires Util.fit package.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fit=TimeDelay.flux_delay_fit('PS',PS,'FreqVec',FreqVec)
% Reliable: 2
%--------------------------------------------------------------------------

InPar = inputParser;

addOptional(InPar,'PS',[]);
addOptional(InPar,'FreqVec',[]);

addOptional(InPar,'VecA1',logspace(-1,1,20).');
addOptional(InPar,'VecA2dA1',logspace(-1,1,20).');

addOptional(InPar,'TauRange',[10 50]);
addOptional(InPar,'TimeSpan',1000);
addOptional(InPar,'gamma',2);
addOptional(InPar,'RelaxGamma',false);
addOptional(InPar,'Norm',1);   % mean(Flux)
addOptional(InPar,'ErrF',0);   % mean(Flux)
addOptional(InPar,'ErrX',0);   % mean(Flux)
addOptional(InPar,'ErrY',0);   % mean(Flux)

addOptional(InPar,'WinPS',[]);   % mean(Flux)
addOptional(InPar,'Limits',[1.8 3;0 Inf;0 Inf; -2 2;-2 2; 5 100; 0 Inf; -2 2; -2 2]);
addOptional(InPar,'FitMethod','fit_perfreq');   % fit_perfreq


parse(InPar,varargin{:});
InPar = InPar.Results;


FreqVec   = InPar.FreqVec;
PS        = InPar.PS;

VecA1     = InPar.VecA1;
VecA2dA1  = InPar.VecA2dA1;
%VecTau    = (3:1:50).';
VecInvTau = (1./max(InPar.TauRange):1./(2.*InPar.TimeSpan):1./min(InPar.TauRange)); 
VecTau    = 1./VecInvTau;
Norm      = InPar.Norm;
ErrF      = [InPar.ErrF, InPar.ErrX, InPar.ErrY];
Limits    = InPar.Limits;

%VecA1     = VecA1.*Norm;

% remove zero frequency from fit + the lowest frequency...
Dfreq = min(abs(diff(FreqVec)));
Fn0 = abs(FreqVec)>(Dfreq+eps);
FreqVec = FreqVec(Fn0);
Omega     = 2.*pi.*FreqVec;
PS      = PS(Fn0);

%MatSigmaF = MatSigmaF(Fn0,:);

if InPar.RelaxGamma
   FitGamma = NaN;
else
    FitGamma = InPar.gamma;
end
   


Data.F_t   = F_t; % - Vector of total flux
Data.Time  = (1:1:numel(F_t)).';
Data.X_t   = x_t; %  - Vector of center of mass X position.
Data.Y_t   = zeros(size(x_t)); % - Vector of center of mass Y position.
%Data.Freq   %- Vector of frequencies
Data.Omega = freqs.*2.*pi; % Omega; %- Vector of 2.*pi.*Frequency
Data.F_w   = fft(F_t); % F_w; %- power spectrum - multiply by 1/sqrt(N)
Err=[sigma_F_t, sigma_X_t, sigma_Y_t];
Err = [5 0.1 0]


switch lower(InPar.FitMethod)
    case 'fit_perfreq'
        
        if InPar.RelaxGamma
            %            [gamma,      A0, A1  x1, y1, Tau2, A2,   x2, y2,...]
            GuessParH1 = [InPar.gamma 1   1   1   0   0     0.66  -1  0];
            GuessParH0 = GuessParH1;
        else
            GuessParH1 = [InPar.gamma 1   1   1   0   0     0.66  -1  0];
            GuessParH0 = GuessParH1;
        end
        
        Ntau = numel(VecTau);
        for Itau=1:1:Ntau
            Tau = VecTau(Itau);
            
            %           [gamma,       A0, A1  x1, y1, Tau2, A2, x2, y2,...]
            FitFlagH1 = [InPar.gamma  NaN NaN NaN NaN Tau   NaN NaN NaN];
            
            FlagNaN   = isnan(FitFlagH1);
            
            GuessParH1c = GuessParH1(FlagNaN);
            
            
            [BestParH1,Fit.LogLH1(Itau),Fit.ExitH1(Itau)] = Util.fit.fminsearch_my({@TimeDelay.ast_delay_logl,Data,Err,FitFlagH1,Limits,InPar.WinPS},GuessParH1c);
            
            
            FitFlagH0 = [InPar.gamma  0   NaN 0   0   0     0   0   0];
            
            [BestParH0,Fit.LogLH0(Itau),Fit.ExitH0(Itau)] = Util.fit.fminsearch_my({@TimeDelay.ast_delay_logl,Data,Err,FitFlagH0,Limits,InPar.WinPS},GuessParH0);
        end
        
        Fit.Tau       = VecTau;
        Fit.InvTau    = 1./VecTau;
        Fit.DeltaLogL = Fit.LogLH1 - Fit.LogLH0;
        [Fit.MinLogLH1, It] = min(Fit.LogLH1);
        Fit.Tau_MinLogLH1 = VecTau(It);
        [Fit.MinDeltaLogLH1, It] = min(Fit.DeltaLogL);
        Fit.Tau_DeltaLogL = VecTau(It);
        
        % rerun best fit
        Tau = Fit.Tau_DeltaLogL;
        FitFlagH1 = [NaN NaN Tau FitGamma];
        FitFlagH0 = [NaN 0   Tau FitGamma];
        [BestParH1,Fit.LogLH1(Itau),E_H1] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlagH1,Limits},GuessParH1);
        [BestParH0,Fit.LogLH0(Itau),E_H1] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlagH0,Limits},GuessParH0);
        Fit.BestParH1 = BestParH1;
        Fit.BestParH0 = BestParH0;
        Fit.gammaInit = InPar.gamma;
        
        % calc grid on best fit
        
        
        VecA1    = VecA1.*    Fit.BestParH1(1);
        
        [MatA1,MatA2dA1] = meshgrid(VecA1,VecA2dA1);
        MatA2 = MatA1.*MatA2dA1;
        
        if InPar.RelaxGamma
            gamma = Fit.BestParH1(end);
        else
            gamma = InPar.gamma;
        end
        
        Omega  = 2.*pi.*FreqVec(:);
        SigmaF = abs(Omega).^(-gamma).*( MatA1(:).'.^2 + MatA2(:).'.^2 + 2.*MatA1(:).'.*MatA2(:).'.*cos(Omega.*Fit.Tau_DeltaLogL) ) + ErrF.^2;

        % log determinant of a diagonal matrix
        LogL = -0.5.*sum(log(2.*pi.*SigmaF),1) - sum(PS(:)./(2.*SigmaF),1);

        % return the -LogL
        Fit.LogL     = reshape(-LogL,numel(VecA2dA1),numel(VecA1));
        Fit.VecA1    = VecA1;
        Fit.VecA2dA1 = VecA2dA1;
        
    
    case 'grid_perfreq'
        
        
        if InPar.RelaxGamma
            %            [gamma,      A0, A1  x1, y1, Tau2, A2,   x2, y2,...]
            GuessParH1 = [InPar.gamma 1   1   1   0   0     0.66  -1  0];
            GuessParH0 = GuessParH1;
        else
            GuessParH1 = [InPar.gamma 1   1   1   0   0     0.66  -1  0];
            GuessParH0 = GuessParH1;
        end
        
        VecX1 = (-2:0.1:2).';
        VecX2 = (-2:0.1:2).';
        VecA1 = (0.1:0.1:2).';
        VecA2dA1 = (0.3:0.1:1).';
        Nx1 = numel(VecX1);
        Nx2 = numel(VecX2);
        Na1 = numel(VecA1);
        Na2a1 = numel(VecA2dA1);
        
        WinPS = [];
        
        Ntau = numel(VecTau);
        for Itau=1:1:Ntau
            Tau = VecTau(Itau);
            
            %           [gamma,       A0, A1  x1, y1, Tau2, A2, x2, y2,...]
            FitFlagH1 = [InPar.gamma  NaN NaN NaN NaN Tau   NaN NaN NaN];
            
            FlagNaN   = isnan(FitFlagH1);
            
            GuessParH1c = GuessParH1(FlagNaN);
            
            LogL = zeros(Nx1,Nx2,Na1,Na2a1);
             for Ix1=1:1:Nx1
                 for Ix2=1:1:Nx2
                     for Ia1=1:1:Na1
                         for Ia2a1=1:1:Na2a1
                             %[gamma, A0, A1 x1, y1, Tau2, A2, x2, y2,...]
                             Par = [3 2 VecA1(Ia1) VecX1(Ix1) 0 Tau VecA1(Ia1).*VecA2dA1(Ia2a1) VecX2(Ix2) 0];
                             FitFlag = [NaN NaN NaN NaN NaN NaN NaN NaN NaN];
                                     Limits = [1.5 3.5;0 Inf;0 Inf; -3 3;-3 3; 5 100; 0 Inf; -3 3; -2 3];

                             LogL(Ix1,Ix2,Ia1,Ia2a1)=TimeDelay.ast_delay_logl(Par,Data,Err,FitFlag,Limits,WinPS);
                             
                         end
                     end
                 end
             end
             
             
            
            [BestParH1,Fit.LogLH1(Itau),Fit.ExitH1(Itau)] = Util.fit.fminsearch_my({@TimeDelay.ast_delay_logl,Data,Err,FitFlagH1,Limits,InPar.WinPS},GuessParH1c);
            
            
            FitFlagH0 = [InPar.gamma  0   NaN 0   0   0     0   0   0];
            
            [BestParH0,Fit.LogLH0(Itau),Fit.ExitH0(Itau)] = Util.fit.fminsearch_my({@TimeDelay.ast_delay_logl,Data,Err,FitFlagH0,Limits,InPar.WinPS},GuessParH0);
        end
        
        
        
        
    case 'fit_bestgrid'

        % greated grid of A1, A2/A1, Tau free parameters
        [MatA1,MatA2dA1,MatTau] = meshgrid(VecA1,VecA2dA1,VecTau);
        MatA2 = MatA1.*MatA2dA1;


        Ngamma = numel(InPar.gamma);
        BestLL = Inf;
        for Igamma=1:1:Ngamma
            gamma     = InPar.gamma(Igamma);

            % calculate SigmaF for all values of free parameters
            MatSigmaF = Norm.*abs(Omega).^(-gamma).*( MatA1(:).'.^2 + MatA2(:).'.^2 + 2.*MatA1(:).'.*MatA2(:).'.*cos(Omega.*MatTau(:).') ) + ErrF.^2;

            % convolve with window function (doesn't work)
            if ~isempty(InPar.WinPS)
                Npar = size(MatSigmaF,2);
                for Ipar=1:1:Npar
                    MatSigmaF(:,Ipar) = conv(MatSigmaF(:,Ipar),InPar.WinPS,'same');
                end
            end

      
            % calculate log likelihhod for all free parameters
            MatLogL = -(-0.5.*sum(log(2.*pi.*MatSigmaF)) - sum(PS./(2.*MatSigmaF)));
            
            [MinLL,MinInd] = min(MatLogL);
            if (MinLL<BestLL)

                [Sll,SIll] = sort(MatLogL);


                %clf;
                %semilogy(MatTau(SIll),Sll,'o')

                BestLL = MinLL;
                bestA1 = MatA1(MinInd);
                bestA2 = MatA2dA1(MinInd);
                bestTau= MatTau(MinInd);
                bestGamma = gamma;
            end
        end

      
        if InPar.RelaxGamma
            % fit gamma
            GuessPar = [sqrt(Norm).*bestA1   bestA1.*bestA2   bestTau, bestGamma];
            gamma    = bestGamma;
            FitFlag = [NaN NaN NaN NaN];

            [BestPar,LogLH1,E_H1,O] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlag,Limits},GuessPar);

            [H] = Util.fit.calc_hessian(@TimeDelay.flux_delay_logl,BestPar,BestPar.*0.3e-1,[FreqVec, PS],ErrF,FitFlag);

            % fit the null hypothesis (A2=0 / Tau=0)
            FitFlag = [NaN 0 0 NaN];
            [BestParH0,LogLH0,E_H0,O] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlag,Limits},GuessPar([1 4]));

        else
            % gamma is constant

            GuessPar = [sqrt(Norm).*bestA1   bestA1.*bestA2   bestTau];
            gamma    = bestGamma;
            FitFlag = [NaN NaN NaN gamma];
            
            [BestPar,LogLH1,E_H1,O] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlag,Limits},GuessPar);

            [H] = Util.fit.calc_hessian(@TimeDelay.flux_delay_logl,BestPar,BestPar.*0.3e-1,[FreqVec, PS],ErrF,FitFlag);

            % fit the null hypothesis (A2=0 / Tau=0)
            FitFlag = [NaN 0 0 gamma];
            [BestParH0,LogLH0,E_H0,O] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlag,Limits},GuessPar(1));

        end

        Fit.BestGrid  = [bestA1, bestA2, bestTau, bestGamma];
        Fit.BestGrifLogL = BestLL;
        Fit.BestPar   = BestPar;
        Fit.gammaInit = gamma;
        Fit.Cov       = inv(0.5.*H);
        Fit.Err       = sqrt(diag(inv(0.5.*H)));
        Fit.BestParH0 = BestParH0;
        Fit.LogLH0    = LogLH0;
        Fit.LogLH1    = LogLH1;
        Fit.DeltaLogL = LogLH0 - LogLH1;
        Fit.DeltaDof  = 2;
        Fit.ExitH1    = E_H1;
        Fit.ExitH0    = E_H0;
        
        
 
    case 'fit2_perfreq'
        
        Ntau = numel(VecTau);
        
        [MatA1,MatA2dA1] = meshgrid(VecA1,VecA2dA1);
        MatA2            = MatA1.*MatA2dA1;

        Nh = numel(MatA1);
        H = [ones(Nh), MatA1(:), MatA1(:).^2, MatA2dA1(:), MatA2dA1(:).^2, MatA1.*MatA2dA1];
        
        Ngamma = numel(InPar.Results.gamma);
        BestLL = Inf;
        for Igamma=1:1:Ngamma
            gamma     = InPar.Results.gamma(Igamma);

            for Itau=1:1:Ntau
                MatSigmaF = Norm.*abs(Omega).^(-gamma).*( MatA1(:).'.^2 + MatA2(:).'.^2 + 2.*MatA1(:).'.*MatA2(:).'.*cos(Omega.*VecTau(Itau).') ) + ErrF.^2;
                MatLogL = -(-0.5.*2.*pi.*sum(log(MatSigmaF)) - sum(PS./(2.*MatSigmaF)));
                
                % fit parabola
                ParPar = H\MatLogL;
            end
        end
        
        
       
        
    otherwise
        error('Unknown FitMethod option');
end
