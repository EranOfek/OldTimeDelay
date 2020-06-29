function [F_w,LC]=rand_psd(VarOmega)
% Generate a complex random frequency realization with real ifft.
% Package: +TimeDelay
% Description: 
% Input  : -
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N=240;
%          Omega=ifftshift(TimeDelay.fft_freq(N))./N .*2.*pi;
%          F_w=TimeDelay.rand_psd(Omega.^-2)
% Reliable: 
%--------------------------------------------------------------------------

N   = numel(VarOmega);
if mod(N,2)==0
    FN2 = floor(N./2)-1;
    F_w = zeros(N,1);
    Freq = (2:1:FN2+1);
    F_w(Freq) = randn(FN2,1).*abs(VarOmega(Freq)).^0.5./sqrt(2);
    F_w(Freq) = F_w(Freq) + 1i.*randn(FN2,1).*abs(VarOmega(Freq)).^0.5./sqrt(2);
    F_w(ceil(N./2 +1)+1:end) = flipud(conj(F_w(Freq)));
    F_w(ceil(N./2 +1))=0;
else
    FN2 = floor(N./2);
    F_w = zeros(N,1);
    Freq = (2:1:FN2+1);
    F_w(Freq) = randn(FN2,1).*abs(VarOmega(Freq)).^0.5./sqrt(2);
    F_w(Freq) = F_w(Freq) + 1i.*randn(FN2,1).*abs(VarOmega(Freq)).^0.5./sqrt(2);
    F_w(ceil(N./2 +1):end) = flipud(conj(F_w(Freq)));
end

if (nargout>1)
    LC = ifft(F_w);
end