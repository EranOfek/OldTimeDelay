function [F_w,LC]=grad_fun(Fun)
% 
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

Fun = @TimeDelay.flux_delay_logl;
AddPar = {PS, ErrF, FitFlag, Limits};

LogL = Fun(Par,AddPar{:});