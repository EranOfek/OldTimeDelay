function test_pl

%%

got to DeltaT=0.1...

% using simulation 11
InPar=select_parameters(11);
        
InPar.Gamma = 2.5;
InPar.A(2) = 0;
InPar.Tau = 0;

Nsim = 1000;

FitPar = [NaN 0 NaN];
DefPar = [1   0 3];
Limits = [1e-5 1000;0 1000;       1.0 4.5];
VecInvTau = Inf;
Min_w = [2.*pi./200, 2.*pi./20];

AllFit = zeros(Nsim,1);
AllLL  = zeros(Nsim,1);


for Isim=1:1:Nsim
    [ResLC,ResG]=generateLC(InPar);
    Flag = ResLC.w>Min_w(1) & ResLC.w<Min_w(2);
    [Par,ParErr,Chi2,Dof] = Util.fit.fitpow(ResLC.w(Flag),abs(ResLC.F_w(Flag)).^2,1);
    AllFit(Isim) = Par(2);
    
    % fit using like
    Res=TimeDelay.fit_flux(ResLC.T,ResLC.F_t,0,'FitPar',FitPar,'DefPar',DefPar,'VecInvTau',VecInvTau,'Limits',Limits,'Min_w',Min_w);
    AllLL(Isim) = Res.BestPar_H1(2);
end

mean(AllFit)
std(AllFit)

mean(AllLL)
std(AllLL)



end % end main fun


%%  Simulation parameters
function InPar=select_parameters(SimName)
% generate parameters for specific simulation number
% Input  : - 11, 15

    InPar.Cyclic = false;
    InPar.x0  = 0;
    InPar.y0  = 0;
    InPar.y   = [0.0  0.0];
    InPar.f_dc = 50;
    InPar.DeltaT  = 1;
    InPar.StdMeanRange = [0.1 0.15];

    switch SimName
        case 11
            % used
            InPar.Tau = 25;
            InPar.A0  = 0;
            InPar.A   = [1 0.5];
            InPar.x   = [0.1 -0.4];
            InPar.Gamma = 2.0;
            InPar.TotTime = 1000;
            InPar.sigma_x = 0.02;
            InPar.sigma_F_rel = 0.02.*sum(InPar.A);
        case 15

            % simulation for non evenly spaced case
            InPar.DeltaT  = 0.1;

            InPar.Tau = 25;
            InPar.A0  = 0;
            InPar.A   = [1 0.5];
            InPar.x   = [0.1 -0.4];
            InPar.Gamma = 2.0;
            InPar.TotTime = 1000;
            InPar.sigma_x = 0.02;
            InPar.sigma_F_rel = 0.02.*sum(InPar.A);

    end
end




function [ResLC,ResG]=generateLC(InPar)


    [ResLC,ResG]=TimeDelay.rand_lensed('A0',InPar.A0,'A',InPar.A,'Tau',InPar.Tau,...
                                    'x0',InPar.x0,'y0',InPar.y0,'x',InPar.x,'y',InPar.y,...
                                    'f_dc',InPar.f_dc,'Gamma',InPar.Gamma,...
                                    'TotTime',InPar.TotTime,...
                                    'DeltaT',InPar.DeltaT,...
                                    'sigma_x',InPar.sigma_x,...
                                    'sigma_F_rel',InPar.sigma_F_rel,...
                                    'Cyclic',InPar.Cyclic,...
                                    'Validate',true,...
                                    'StdMeanRange',InPar.StdMeanRange);
end





