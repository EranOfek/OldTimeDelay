function [Res,NES]=rand_lc_tau(varargin)
% Generate a light curve generated from a combination of time-delayed LC.
% Package: +TimeDelay
% Description: Generate a light curve generated from a combination of
%              time-delayed light curves with a specific power spectrum.
% Input  : * Arbibtrary number of ...,key,val,...
%            The following keywords are available:
%            'TimeVec' - Default is (1:240).
%            'TimeResample' - Default is []. If provided, then the equally
%                             spaced time series will be resample to this 
%                             not equally spaced (NES) time series.
%            'InterpMethod' - Default is 'pchip'.
%            'A0'      - Default is 0.
%            'A'       - Default is [1 2./3]
%            'Tau'     - Default is [23.1]
%            'gamma'   - Default is 2.
%            'fDC'     - Default is 4.
%            'sigmaFprop' - Default is 0.01.
%            'NoiseType'  - Default is 'relative_to_relstd'
% Output : - A structure with the following fields:
%            'phi_t' - Noisless light curve as a function of time.
%            'F_t'   - Light curve including noise.
%            'phi_w' - FT of phit_t.
%            'F_w'   - FT of F_t.
%            'Omega' - 2.*pi.*Freq of FT.
%            'Time'  - 
%            'MeanF'
%            'RelStdF'
%            'sigmaFhat'
%            'sigmaF'
%          - A structure array with the non equally spaced time series.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,NES]=TimeDelay.rand_lc_tau
% Reliable: 
%--------------------------------------------------------------------------


InPar = inputParser;

addOptional(InPar,'TimeVec',(1:240).');
addOptional(InPar,'OverFactor',false);
addOptional(InPar,'TimeResample',[]);
addOptional(InPar,'InterpMethod','pchip');

addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 2./3]);
addOptional(InPar,'Tau',[23.1]);
addOptional(InPar,'gamma',2);
addOptional(InPar,'fDC',4);
addOptional(InPar,'sigmaFprop',0.01);
addOptional(InPar,'NoiseType','relative_to_relstd');  % relative_to_flux  | relative_to_std | relative_to_relstd

parse(InPar,varargin{:});
InPar = InPar.Results;

DT = diff(InPar.TimeVec);
if (isempty(InPar.TimeResample) && range(DT)>(10.*eps))
    error('Non equally spaced time series');
end

DT = DT(1);

Nt = numel(InPar.TimeVec);

if InPar.OverFactor
    OverFactor = 10;
    
    Nover = OverFactor.*Nt;
else
    Nover = Nt;
end


Omega   = ifftshift(TimeDelay.fft_freq(Nover)).*2.*pi./(Nover.*DT);
[f_red] = TimeDelay.rand_psd(Omega.^-InPar.gamma);

f_red(1) = InPar.fDC.*Nover;   %???

f_multi = f_red.*(InPar.A(1) + InPar.A(2).*exp(1i.*Omega.*InPar.Tau(1)));
% add galaxy (constant) at frequency zero
f_multi(1) = f_multi(1) + InPar.A0.*Nover;   %???   % phi

Res.phi_t = ifft(f_multi);   % noisless LC of all images


% calculations with OverFactor light curve are done
% return to original light curve
Res.phi_t = Res.phi_t(1:Nt);

%sigmaF    = InPar.sigmaFprop.*sum(InPar.A).*InPar.fDC;
%sigmaFhat = sigmaF.*sqrt(Nover);  % ???  note we are using N as this is the overfactor sampled series
% recalculate omega without OverFactor
Omega   = ifftshift(TimeDelay.fft_freq(Nt)).*2.*pi./(Nt.*DT);


switch lower(InPar.NoiseType)
    case 'relative_to_flux'
        % noise is already relative to Flux
    case 'relative_to_std'
        StdPhi = std(Res.phi_t);
        
        sigmaF = InPar.sigmaFprop.*StdPhi;
        sigmaFhat = sigmaF.*sqrt(Nt);  %???
        
    case 'relative_to_relstd'
        RelStdPhi = std(Res.phi_t)./mean(Res.phi_t);
        sigmaF = InPar.sigmaFprop.*RelStdPhi;
        sigmaFhat = sigmaF.*sqrt(Nt);   %???
    otherwise
        error('Unknown NoiseType option');
end


Res.F_t   = Res.phi_t + randn(Nt,1).*sigmaF; %.*mean(Res.phi_t); %   sum(InPar.A).*InPar.fDC;
%Res.F_t   = Res.phi_t + randn(N,1).*InPar.sigmaF.*sum(InPar.A).*InPar.fDC;
Res.phi_w = fft(Res.phi_t);
Res.F_w   = fft(Res.F_t);
Res.Omega = Omega;
Res.Time  = InPar.TimeVec;
Res.MeanF = mean(Res.F_t);
Res.StdF  = std(Res.F_t);
Res.RelStdF = Res.StdF./Res.MeanF;
Res.NoiseToStd = sigmaF./Res.StdF;
Res.sigmaFhat = sigmaFhat;
Res.sigmaF    = sigmaF;


% deal with non-evenly spaced time series
NES = [];
if nargout>1 && ~isempty(InPar.TimeResample)
    NES.F_t   = interp1(InPar.TimeVec,Res.F_t,InPar.TimeResample,InPar.InterpMethod);
    NES.Time  = InPar.TimeResample;
    
    % calculate the Fourier Transform
    NES.F_w   = sum(NES.F_t(:).*exp(-1i.*Omega(:).'.*NES.Time(:)),1).';
    NES.Omega = Omega;
    
end

