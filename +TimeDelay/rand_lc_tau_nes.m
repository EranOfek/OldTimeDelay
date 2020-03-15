function [Res,PS]=rand_lc_tau_nes(varargin)
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
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=TimeDelay.rand_lc_tau_nes
% Reliable: 
%--------------------------------------------------------------------------


InPar = inputParser;

addOptional(InPar,'TimeVec',(1:240).');
addOptional(InPar,'InterpMethod','pchip');
addOptional(InPar,'MaxFreq',[]);

addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 0.4]);
addOptional(InPar,'Tau',[36.9]);
addOptional(InPar,'gamma',2);
addOptional(InPar,'Std2Flux',0.3);
addOptional(InPar,'Noise2Std',0.03);
addOptional(InPar,'CalcPS',true);


addOptional(InPar,'NoiseType','relative_to_relstd');  % relative_to_flux  | relative_to_std | relative_to_relstd

parse(InPar,varargin{:});
InPar = InPar.Results;

MinDT = min(diff(InPar.TimeVec));

ES_TimeVec = (min(InPar.TimeVec)-2.*InPar.Tau:0.5.*MinDT:max(InPar.TimeVec)+2.*InPar.Tau).';
N_ES       = numel(ES_TimeVec);    % number of ES points
N          = numel(InPar.TimeVec); % number of points

TS = Util.stat.rand_ps(ES_TimeVec,[InPar.gamma 1],zeros(N_ES,1));
% interpolate LC
LC1 = interp1(TS(:,1),TS(:,2),TS(:,1)          , InPar.InterpMethod);
LC2 = interp1(TS(:,1),TS(:,2),TS(:,1)+InPar.Tau, InPar.InterpMethod);

% plot(TS(:,1),LC1);
% hold on;
% plot(TS(:,1),LC2);

LC  = InPar.A(1).*LC1 + InPar.A(2).*LC2;
% interpolate into NES grid
LC_NES = interp1(TS(:,1),LC,InPar.TimeVec,InPar.InterpMethod);


StdLC = std(LC_NES);
% add constant flux
LC_NES = LC_NES + StdLC./InPar.Std2Flux;
% add noise
LC_NES = LC_NES + StdLC.*InPar.Noise2Std.*randn(N,1);

% plot(InPar.TimeVec,LC_NES);

Res.TimeVec = InPar.TimeVec;
Res.F_t     = LC_NES;
Res.SigmaF  = StdLC.*InPar.Noise2Std;
Res.MeanF   = mean(Res.F_t);
Res.StdF    = std(Res.F_t);
Res.sigmaF  = StdLC.*InPar.Noise2Std;
Res.N       = N;
Res.N_ES    = N_ES;

% calculate the power spectrum
DeltaFreq = 0.5./range(Res.TimeVec);
if isempty(InPar.MaxFreq)
    % choose max freq from data
    InPar.MaxFreq = 0.5./median(diff(InPar.TimeVec));
end

if (InPar.CalcPS)
    FreqVec  = (0:DeltaFreq:InPar.MaxFreq).';
    FreqVec  = [-flipud(FreqVec(2:end)); FreqVec];
    Res.Nfreq = numel(FreqVec);
    Res.F_w   = sum(Res.F_t(:).*exp(-2.*pi.*1i.*FreqVec(:).'.*Res.TimeVec(:)),1).'./2;
    Res.Freq  = FreqVec;
    Res.Omega = 2.*pi.*Res.Freq;
    Res.sigmaFhat = Res.SigmaF.*sqrt(Res.N);
end





