function [FT,FreqVec,FT2]=fft_withphase(Time,Data)
% Calculate the log-likelihood for two images flux time delay fit
% Package: TimeDelay
% Description: Calculate the log-likelihood for two images flux time delay
%              fit.
% Input  : - Vector of free paramaeters [A1, A2, Tau [gamma]]
%            gamma is optional.
%          - A two column vector of [Frequency, PowerSpectrum].
%          - A flux error scalar.
%          - gamma. Default is 2.
% Output : - Minus Log likelihhod.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=rand(100,1); [FT,FreqVec,FT2] = TimeDelay.fft_withphase([],Data);
% Reliable: 
%--------------------------------------------------------------------------

Nt = numel(Data);

if isempty(Time)
    Time = (0:1:Nt-1).'; %(-Nt./2:1:(Nt-1)./2).';
end

TimeRange = range(Time);

FT = fft(Data);
FreqVec = ifftshift(TimeDelay.fft_freq(Nt)./Nt);
Nf      = numel(FreqVec);
%FV      = (0:1:Nf-1);
%TV      = (0:1:Nt-1).';

if (nargout>2)
    FT2 = sum(Data.*exp(-1i.*2.*pi.*FreqVec.'.*Time),1).';
end
