function FreqVec=ast_delay_fit_simple(T,F,x,y,errF_t,errx_t,erry_t)
% 
% Package: +TimeDelay
% Description: 
% Input  : - Number of points.
% Output : - Vector of frequencies.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FreqVec=TimeDelay.fft_freq(4)
% Reliable: 2
%--------------------------------------------------------------------------

if nargin==0
    % simulation mode
    gamma = 2.5;
    T = (1:1:240).';
    TS = Util.stat.rand_ps(T,[gamma 1],zeros(size(T)));
    f1 = TS(:,2);
    errF_t = 0.01;
    errx_t = 0.01;
    erry_t = 0.01;
    xVec = [1 -1];
    yVec = [0  0.0];
    A0   = 0.1;
    AlphaVec = [1 0.5];
    Tau      = 14.7;
    
    ResFs = TimeDelay.photometric_astrometric_model(T,f1,'A0',A0,...
                                                        'A',AlphaVec,...
                                                        'Tau',Tau,...
                                                        'x0',0,...
                                                        'y0',0,...
                                                        'x',xVec,...
                                                        'y',yVec,...
                                                        'eps_x',errx_t,...
                                                        'eps_F',errF_t,...
                                                        'gamma',gamma);
    
    %
    T   = ResFs.T;
    F_t = ResFs.F + randn(size(T)).*errF_t;
    x_t = ResFs.X + randn(size(T)).*errx_t;
    y_t = ResFs.Y + randn(size(T)).*erry_t;
    
end

x0 = 0;
y0 = 0;
x1vec = (-2:0.1:2).';
x2vec = (-2:0.1:2).';
a0vec = (0:0.1:1).';
a1vec = (0.1:0.1:2).';
a2vec = (0.1:0.1:2).';
tauvec= (10:1:50).';
Nx1   = numel(x1vec);
Nx2   = numel(x2vec);
Na0   = numel(a0vec);
Na1   = numel(a1vec);
Na2   = numel(a2vec);
Ntau  = numel(tauvec);
Nx1.*Nx2.*Na0.*Na1.*Na2.*Ntau


        tic;
f1_rec = TimeDelay.ft_from_Ft_xt(F_t,x_t,errF_t,errx_t,[0 xVec],[A0 AlphaVec]);

ResF = TimeDelay.photometric_astrometric_model(T,f1_rec,'A0',A0,...
                                                        'A',AlphaVec,...
                                                        'Tau',Tau,...
                                                        'x0',x0,...
                                                        'y0',y0,...
                                                        'x',xVec,...
                                                        'y',yVec,...
                                                        'eps_x',errx_t,...
                                                        'eps_F',errF_t,...
                                                        'gamma',gamma);
                                                    

   [X,Fval]=Util.fit.fminsearch_my({Fun,2},0.7);                                                  
%
Nt = numel(ResF.T);
No = numel(F_t);
Resid = ResF.F - F_t(No-Nt+1:end);

std(Resid)

toc

