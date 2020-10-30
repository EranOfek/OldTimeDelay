function LC=random_lc(varargin)
% Return a red-power-spectrum random LCs with possible time delays
% Package: +TimeDelay
% Description: Construct a random light curve which contains several
%              shifted components with light curve generated from a red
%              power-law power spectrum process (omega^-gamma).
%              For example, this can be used to generate the light curve of
%              a lensed quasar.
% Input  : * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Time' - If scalar, then this is the number of points in an
%                     equally light curve.
%                     If vector this is a non-equally spaced light curve
%                     interpolated from an equally spaced light curve with
%                     sampling of 'SampT'.
%            'SampT' - Time sampling. Default is 1.
%                     This is used when the equaly spaced light curve is
%                     generated.
%            'gamma' - Power-law power spectrum: Amp ~ Frequency^-gamma.
%                     Default is 2.
%            'A0'    - The galaxy (constant light) flux.
%                     Default is 0.
%            'A'     - Vector of the mean fluxes of the various components.
%                      Default is [1 2./3]
%            'TimeDelay' - Vector of time delays for each image.
%                      Default is [0 14.7].
%            'fDC'   - DC flux component. Default is 7.
%            'Err'   - Relative error to add to (only) the combined light
%                      curve.
%            'InterpMethod' - Default is 'linear'.
% Output : - A structure containing the following fields
%            'RelErr' - Relative error.
%            'ErrF'   - Flux error.
%            'LC1_Flux'- Flux of the primary component without noise.
%            'LC_Flux' - Combined flux of all components with noise.
%            'LC_FT  ' - shifted FFT of the FFT of the equaly spaced LC.
%            'LC_winFT'- window function of FFT (only for non equaly
%                        spaced).
%            'FreqVec' - Vector of frequencies corresponding to LC_sFFT.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: LC=TimeDelay.random_lc
%          TR=timeseries.random_time_sequence(700,1,250,0,0.8);
%          LC=TimeDelay.random_lc('Time',TR);
% Reliable: 2
%--------------------------------------------------------------------------

InPar = inputParser;

addOptional(InPar,'Time',731);
addOptional(InPar,'gamma',2);
addOptional(InPar,'A0',0);
addOptional(InPar,'A',[1 2./3]);
addOptional(InPar,'TimeDelay',[0 14.7]);
addOptional(InPar,'fDC',7);
addOptional(InPar,'Err',0.03);
addOptional(InPar,'SampT',1);
addOptional(InPar,'InterpMethod','linear');

parse(InPar,varargin{:});


% light curves times
if numel(InPar.Results.Time)==1
    % scalar option
    % number of equay spaced data points
    
    TimeES = (1:InPar.Results.SampT:InPar.Results.Time).';
    Time   = TimeES;
    IsEqSpaced = true;
else
    % vector option - assume non equaly spaced
    Time = InPar.Results.Time;
    MinT = min(Time);
    MaxT = max(Time);
    
    TimeES = (MinT:InPar.Results.SampT:MaxT+InPar.Results.SampT).';
    IsEqSpaced = false;
end

Nt      = numel(TimeES);  % number of equaly spaced times
RangeT  = range(TimeES);
FreqVec = TimeDelay.fft_freq(Nt)./Nt; %RangeT;  % used to be RangeT!!!
FreqVec = FreqVec(:);
Nf      = numel(FreqVec);

% angular frequency
Omega   = 2.*pi.*FreqVec;
% zero frequency
IF0     = find(FreqVec==0);

% shifted versions of the frequency
% zero frequency is in position=1
shOmega = ifftshift(Omega);     
shFreqVec = ifftshift(FreqVec);
    

% genertae random power low light curve
RandPhase = rand(Nf,1).*2.*pi;
%RandAmp   = randn(Nf,1).*FreqVec.^(-0.5.*InPar.Results.gamma);
RandAmp   = randn(Nf,1).*abs(2.*pi.*FreqVec).^(-0.5.*InPar.Results.gamma);

% FFT of single random LC 
LC1_FFT      = RandAmp.*(cos(RandPhase) + 1i.*sin(RandPhase));
LC1_FFT(IF0) = Nt.*InPar.Results.fDC;
LC1_sFFT     = ifftshift(LC1_FFT);
LC1_Flux     = ifft(LC1_sFFT,'symmetric'); % complex conjugate symmetric ifft

% combined sources
% shOmega and LC1_sFFT is a column vector by construction
LC_sFFT    = [InPar.Results.A(:).'.*exp(1i.*shOmega.*InPar.Results.TimeDelay(:).')].*LC1_sFFT;
% sum all images
LC_sFFT    = sum(LC_sFFT,2);
LC_sFFT(1) = LC_sFFT(1) + Nt.*InPar.Results.A0;   % add galaxy light (A0)
LC_Flux    = ifft(LC_sFFT,'symmetric'); % complex conjugate symmetric ifft
    


if ~IsEqSpaced
    % interpolate back to non equlay spaced
    LC.Time     = Time;
    LC.LC1_Flux = interp1(TimeES,LC1_Flux,Time,InPar.Results.InterpMethod);
    LC.LC_Flux  = interp1(TimeES,LC_Flux, Time,InPar.Results.InterpMethod);
    
    % calculate the non-equaly spaced Fourier Transform
    LC.LC_FT    = sum(LC.LC_Flux(:).*exp(-1i.*shOmega(:).'.*LC.Time(:)),1).';
    % window function of fft
    LC.LC_winFT = sum(exp(-1i.*shOmega(:).'.*LC.Time(:)),1).';
else
    % equaly spaced
    LC.Time     = Time;
    LC.LC1_Flux = LC1_Flux;
    LC.LC_Flux  = LC_Flux;
    LC.LC_winFT = [];
    %LC.LC_FT    = fft(LC.LC_Flux);
end

LC.FreqVec  = shFreqVec;

% if ~IsEqSpaced
%     LC.LC_Flux = interp1(LC.Time, LC.LC_Flux, TimeES, InPar.Results.InterpMethod);
%     LC.Time = TimeES;
%     NN = ~isnan(LC.LC_Flux);
%     LC.LC_Flux = LC.LC_Flux(NN);
%     LC.Time    = LC.Time(NN);
%     Nt         = numel(LC.Time);
%     LC.FreqVec = ifftshift(TimeDelay.fft_freq(Nt)./RangeT);
%     
%     LC.LC_FT = fft(LC.LC_Flux);
%     LC.LC_winFT = [];
% else
%     LC.FreqVec  = shFreqVec;
% end

LC.RelErr   = InPar.Results.Err;
LC.ErrF     = mean(LC.RelErr .* LC.LC_Flux);
LC.LC_Flux  = LC.LC_Flux + LC.ErrF.*randn(numel(LC.LC_Flux),1);
LC.LC_FT    = fft(LC.LC_Flux);