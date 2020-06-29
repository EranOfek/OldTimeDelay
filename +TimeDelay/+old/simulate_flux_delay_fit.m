function [Fit,Res]=simulate_flux_delay_fit(varargin)
% 
% Package: +TimeDelay
% Description: 
% Input  : -
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Fit,Res]=TimeDelay.simulate_flux_delay_fit
%          [Fit,Res]=TimeDelay.simulate_flux_delay_fit('EqSpaced',false);
% Reliable: 
%--------------------------------------------------------------------------


InPar = inputParser;

addOptional(InPar,'TimeVec',(1:1024).');
addOptional(InPar,'Cyclic',false);
addOptional(InPar,'EqSpaced',true);
addOptional(InPar,'InterpMethod','pchip');
addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 0.4]); %2./3]);
addOptional(InPar,'Tau',[36.9]);
addOptional(InPar,'gamma',2.5);
addOptional(InPar,'RelaxGamma',false);
addOptional(InPar,'fDC',1);
addOptional(InPar,'Std2Flux',0.5);
addOptional(InPar,'Noise2Std',0.01);

addOptional(InPar,'sigmaFprop',0.01);
addOptional(InPar,'NoiseType','relative_to_relstd');
addOptional(InPar,'Plot',true);
addOptional(InPar,'Limits',[0 100; 0 1; 3 120; 1.5 2.5]);




parse(InPar,varargin{:});
InPar = InPar.Results;

VecA1   = logspace(-1,1,40).'; %(0:0.05:10).';
VecA2dA1= logspace(-2,2,30).'; %(0:0.05:5).';


if InPar.Cyclic
    [Res] = TimeDelay.rand_lc_tau('TimeVec',InPar.TimeVec,...
                                  'TimeResample',[],...
                                  'InterpMethod',InPar.InterpMethod,...
                                  'A0',InPar.A0,...
                                  'A',InPar.A,...
                                  'Tau',InPar.Tau,...
                                  'gamma',InPar.gamma,...
                                  'fDC',InPar.fDC,...
                                  'sigmaFprop',InPar.sigmaFprop,...
                                  'NoiseType',InPar.NoiseType);
    % relative_to_flux  | relative_to_std | relative_to_relstd
    PS = abs(Res.F_w).^2;
    FreqVec = Res.Omega./(2.*pi);
    
    WinPS = [];
else
    % non cyclic simulations
    
    if InPar.EqSpaced
        % non cyclic equally spaced
        [Res]=TimeDelay.rand_lc_tau_nes('TimeVec',InPar.TimeVec,...
                                           'A',InPar.A,...
                                           'Tau',InPar.Tau,...
                                           'gamma',InPar.gamma,...
                                           'Std2Flux',InPar.Std2Flux,...
                                           'Noise2Std',InPar.Noise2Std,...
                                           'CalcPS',true);
                                       
        PS = abs(Res.F_w).^2;
        FreqVec = Res.Omega./(2.*pi);    
        
        WinPS = [];
                                  
    else
        % SimN3
        %InPar.TimeVec = timeseries.random_time_sequence(5.*365,1,240,0.05,0.9);
        % SimN4
        InPar.TimeVec = timeseries.random_time_sequence(5.*365,2,240,0.05,0.8);
        
        [Res]=TimeDelay.rand_lc_tau_nes('TimeVec',InPar.TimeVec,...
                                           'A',InPar.A,...
                                           'Tau',InPar.Tau,...
                                           'gamma',InPar.gamma,...
                                           'Std2Flux',InPar.Std2Flux,...
                                           'Noise2Std',InPar.Noise2Std,...
                                           'CalcPS',true);      
                                       
           
        UseSection = true;
        if UseSection
            [PS,FT,Section]=TimeDelay.power_spectrum_sections([Res.TimeVec, Res.F_t]);   
            FreqVec = PS(:,1);
            PS      = PS(:,2);
            WinPS   = [];
        else
            % use window function
            PS = abs(Res.F_w).^2;
            FreqVec = Res.Omega./(2.*pi);    

            WinPS = abs(Res.W_w).^2;
        end
    end

end

TimeSpan  = range(InPar.TimeVec); %3.*365; %1731;
  
    
%ErrF      = InPar.sigmaFprop.*Res.MeanF; %.*sqrt(Npt);



Flag = abs(FreqVec)<0.25;

Fit=TimeDelay.flux_delay_fit('PS',PS(Flag),'FreqVec',FreqVec(Flag),'Err',Res.sigmaFhat,'WinPS',WinPS,...
                'VecA1',VecA1,'VecA2dA1',VecA2dA1,'gamma',InPar.gamma,...
                'TauRange',[3 120],'TimeSpan',TimeSpan,'RelaxGamma',InPar.RelaxGamma,...
                'Limits',InPar.Limits,'FitMethod','fit_perfreq');

% Fit1=TimeDelay.flux_delay_fit('PS',PS(Flag),'FreqVec',FreqVec(Flag),'Err',Res.sigmaFhat.*3,'winPS',[],...
%                 'VecA1',VecA1,'VecA2dA1',VecA2dA1,'gamma',InPar.gamma,...
%                 'TauRange',[3 120],'TimeSpan',TimeSpan,'RelaxGamma',InPar.RelaxGamma,...
%                 'Limits',InPar.Limits,'FitMethod','fit_perfreq');
%             
% [Fit.MinLogLH1, Fit1.MinLogLH1]
% [Fit.Tau_DeltaLogL, Fit1.Tau_DeltaLogL]

            %'fit_bestgrid'); %'fit_perfreq');
            


if InPar.Plot
    Omega = FreqVec.*2.*pi;
    SigmaF = abs(Omega).^(-InPar.gamma).*( InPar.A(1).^2 + InPar.A(2).^2 + 2.*InPar.A(1).*InPar.A(2).*cos(Omega.*InPar.Tau(1)) ) + Res.sigmaFhat.^2;
    
    BestA1 = Fit.BestParH1(1);
    BestA2 = Fit.BestParH1(2);
    BestTau = Fit.Tau_DeltaLogL;
    if (numel(Fit.BestParH1)>2)
        BestGamma = Fit.BestParH1(3);
    else
        BestGamma = Fit.gammaInit;
    end
    SigmaF_calc = abs(Omega).^(-BestGamma).*( BestA1.^2 + BestA2.^2 + 2.*BestA1.*BestA2.*cos(Omega.*BestTau) ) + Res.sigmaFhat.^2;
    
    % add NaN between negative and positive frequencies
    I0 = find(FreqVec<0,1,'first');
    FreqVec     = [FreqVec(1:I0-1); NaN; FreqVec(I0:end)];
    SigmaF      = [SigmaF(1:I0-1); NaN; SigmaF(I0:end)];
    SigmaF_calc = [SigmaF_calc(1:I0-1); NaN; SigmaF_calc(I0:end)];
    PS          = [PS(1:I0-1); NaN; PS(I0:end)];

    semilogy(FreqVec,PS)
    hold on;
    semilogy(FreqVec,SigmaF,'-','LineWidth',2)
    
    semilogy(FreqVec,SigmaF_calc,'-','LineWidth',2)
    %semilogy(FreqVec,SigmaF_calc.*0.5,'-','LineWidth',2)
    %PS = sum([Section.CommonPS],2);

    H=legend('PS','$\Sigma_{F}$ actual','$\Sigma_{F}$ fit','Location','NorthEast');
    H.Interpreter = 'latex';
end
   
