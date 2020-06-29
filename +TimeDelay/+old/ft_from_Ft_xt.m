function f_t=ft_from_Ft_xt(F_t,x_t,errF_t,errx_t,xVec,AlphaVec)
% Reconstruct f(t) in time domain, from F(t) and x(t)
% Package: +TimeDelay
% Description: Given the total flux of a two image system (F(t)), and its
%              center of light position (x(t)), reconstruct f(t) in the
%              time domain, from F(t) and x(t), the errors and the flux
%              ratios and images position.
% Input  : - Total flux of a two image system and a host, as a function of
%            time - F(t).
%          - Center of light position of the system, as a function of
%            time - x(t).
%          - Vector of total flux error as a function of time.
%          - Vector of total centeroid errots as a function of time.
%          - A vector of [x0, x1, x2]. In principle each element is a two D
%            vector, but here they are treated as scalars.
%          - A vector of [Alpha0, Alpha1, Alpha2].
% Output : -  The light curve of the individual source f(t).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------


Alpha0 = AlphaVec(1);
Alpha1 = AlphaVec(2);
Alpha2 = AlphaVec(3);
x0     = xVec(1);
x1     = xVec(2);
x2     = xVec(3);

f_t = ((F_t - errF_t).*(x_t - x2 - errx_t) + Alpha0.*(x2-x0))./(Alpha1.*(x1 - x2));
