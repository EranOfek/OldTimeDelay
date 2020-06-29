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
addOptional(InPar,'Err',0);   % mean(Flux)
addOptional(InPar,'WinPS',[]);   % mean(Flux)
addOptional(InPar,'Limits',[0 Inf; -eps 1; 5 100; 1.5 3.5]);   % mean(Flux)
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
ErrF      = InPar.Err;
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
   

switch lower(InPar.FitMethod)
    case 'fit_perfreq'
        
        if InPar.RelaxGamma
            GuessParH1 = [1 1 InPar.gamma];
            GuessParH0 = [1 InPar.gamma];
        else
            GuessParH1 = [1 1];
            GuessParH0 = [1];
        end
        
        Ntau = numel(VecTau);
        for Itau=1:1:Ntau
            Tau = VecTau(Itau);
            FitFlagH1 = [NaN NaN Tau FitGamma];
            FitFlagH0 = [NaN 0   Tau FitGamma];
            [BestParH1,Fit.LogLH1(Itau),Fit.ExitH1(Itau)] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlagH1,Limits,InPar.WinPS},GuessParH1);
            
            [BestParH0,Fit.LogLH0(Itau),Fit.ExitH0(Itau)] = Util.fit.fminsearch_my({@TimeDelay.flux_delay_logl,[FreqVec, PS],ErrF,FitFlagH0,Limits,InPar.WinPS},GuessParH0);
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
