function plots_for_paper
%%

T = (1:1:300).';
Nt = numel(T);
gamma = 3;

f1 = TimeDelay.simulate_lc_from_ps(T,[gamma 1],zeros(Nt,1));
f1 = f1(:,2);

InPar.A0 = 0;
InPar.A  = [1 0.5];
InPar.Tau = 14;
InPar.x0  = 0;
InPar.y0  = 0;
InPar.x   = [1 -1];
InPar.y   = [0 0];
InPar.StdF = 0.2;
InPar.eps_x_abs = 0.005;
InPar.eps_F_rel = 0.005;


ResLC=TimeDelay.timedelayed_lc(T,f1,'A0',InPar.A0,'A',InPar.A,...
                          'Tau',InPar.Tau,...
                          'x0',InPar.x0,'y0',InPar.y0,...
                          'x',InPar.x,'y',InPar.y,...
                          'StdF',InPar.StdF',...
                          'eps_x_abs',InPar.eps_x_abs,'eps_F_rel',InPar.eps_F_rel);
                     
FitPar = [InPar.A0   InPar.A(1)  InPar.A(2)  InPar.x0   InPar.x(1)   InPar.x(2)   InPar.y0   InPar.y(1)   InPar.y(2)    gamma];  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
VecA1 = (0.5:0.05:1.5);
VecA2dA1 = (0.4:0.05:1);

                      
Res=TimeDelay.fit_scan_alpha_astrometric_flux(ResLC.T, ResLC.F_t,ResLC.x_t,ResLC.y_t,ResLC.eps_F_abs,ResLC.eps_x_abs,...
                    'Tau',InPar.Tau,'FitPar',FitPar,'VecA1',VecA1,'VecA2dA1',VecA2dA1);
                

contour(Res.A1,Res.A2dA1,Res.LL_H1);

