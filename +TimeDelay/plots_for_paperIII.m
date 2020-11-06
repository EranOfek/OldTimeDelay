function plots_for_paperIII
%% single light curve example
SimName = 111;

Fun = @(Par,X) Par(1).*X.^(-Par(2));
FunChi2 = @(Par,X,Yobs,Err) sum(((Yobs - Par(1).*X.^(-Par(2)))./Err).^2);
FunChi2 = @(Par,X,Yobs,Err) sum(((Yobs - (Par(1).*X.^(-Par(2)) + Err.^2)  )./(Err.^2)).^2);
FunChi2 = @(Par,X,Yobs,Err) sum(((Yobs - (Par(1).*X.^(-Par(2)) )  )./(Err.^2)).^2);

%FunChi2 = @(Par,X,Yobs,Err) sum(((Yobs - (Par(1).*X.^(-Par(2)) + Err.^2)  )./(1)).^2);
%FunChi2 = @(Par,X,Yobs,Err) sum(((Yobs - (Par(1).*X.^(-Par(2)) + Err.^2)  )./(1)).^2);

%FunChi2 = @(Par,X,Yobs,RelErr) sum(( (Yobs - (Par(1).*X.^(-Par(2))+Yobs.*RelErr)  )./(Yobs.*RelErr)).^2);
FunChi2wE = @(Par,X,Yobs) sum(((Yobs - Par(1).*X.^(-Par(2)))./(Yobs.*Par(3))).^2);


InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = false;
InPar.AliasFactor = 1;

[ResLC,ResG]=TimeDelay.rand_lensed(InPar);

ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));

Flag = ResLC.w>0;

loglog(ResLC.w(Flag),abs(ResLC.F_w(Flag)).^2)
hold on;

% max like fit
ResLL = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
loglog(ResLC.w(Flag), ResLL.BestPar_H0(1).*ResLC.w(Flag).^(-ResLL.BestPar_H0(2)))


% simple power-law fit to the PS
X = log10(ResLC.w(Flag));
Y = log10( abs(ResLC.F_w(Flag)).^2);
Par = polyfit(X,Y,1);
Yfit = polyval(Par,X);
loglog(10.^X, 10.^Yfit)


% non linear fit with errors
%[Beta,Chi2,E,O]=Util.fit.fminsearch_chi2(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, ResLC.sigma_F_hat,{Fun},[1 2]);
%loglog(ResLC.w(Flag), Beta(1).*ResLC.w(Flag).^(-Beta(2)))

RelErr = ResLC.sigma_F_rel;
Err    = ResLC.sigma_F_hat;
FlagNL = ResLC.w>0 & ResLC.w<3e-1;
[Beta,Chi2,E,O]=Util.fit.fminsearch_my({FunChi2, ResLC.w(FlagNL), abs(ResLC.F_w(FlagNL)).^2, Err}, [1 2]);

% RelErrVec = [0.001 0.003 0.01 0.02 0.03 0.05 0.1 0.3 1 3 10];
% Nre = numel(RelErrVec)
% for Ire=1:Nre
%     RelErr = RelErrVec(Ire);
% end

%[Beta,Chi2,E,O]=Util.fit.fminsearch_my({FunChi2wE, ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2}, [1.5 2 0.02]);


loglog(ResLC.w(Flag), Beta(1).*ResLC.w(Flag).^(-Beta(2)) )




%[Par,ParErr,Chi2,Dof] = Util.fit.fitpow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1);
%loglog(ResLC.w(Flag), Par(1).*ResLC.w(Flag).^(-Par(2)))

% fit using non-lin chi^2
%[ResNL]=Util.fit.fit_pow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1,'Par0',[1 -2]);
%loglog(ResLC.w(Flag), ResNL.Par(1).*ResLC.w(Flag).^(ResNL.Par(2)))

%%

SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = false;
InPar.AliasFactor = 1;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

%     % linear fit to all data points
%     X = log10(ResLC.w(Flag));
%     Y = log10( abs(ResLC.F_w(Flag)).^2);
%     Par_Lin(Isim,:) = polyfit(X,Y,1);
%     
%     
%     % linear fot to binned data points
%     BinVec = logspace(log10(min(ResLC.w(Flag))), log10(max(ResLC.w(Flag))),30);
%     B      = timeseries.binning([ResLC.w(Flag), ResLC.F_w(Flag)], BinVec,[NaN NaN],{'MeanBin',@mean,@std,@numel});
%     
%     % non linear fit with errors
%     [Beta,Chi2,E,O]=Util.fit.fminsearch_my({FunChi2, ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, ResLC.sigma_F_hat}, [1 2]);
%     Par_NLe(Isim,:) = Beta;
    
    
end

save Res_PL_NEM_Alias1.mat Res ParLL InPar

%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 1;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

%     % linear fit to all data points
%     X = log10(ResLC.w(Flag));
%     Y = log10( abs(ResLC.F_w(Flag)).^2);
%     Par_Lin(Isim,:) = polyfit(X,Y,1);
%     
%     
%     % linear fot to binned data points
%     BinVec = logspace(log10(min(ResLC.w(Flag))), log10(max(ResLC.w(Flag))),30);
%     B      = timeseries.binning([ResLC.w(Flag), ResLC.F_w(Flag)], BinVec,[NaN NaN],{'MeanBin',@mean,@std,@numel});
%     
%     % non linear fit with errors
%     [Beta,Chi2,E,O]=Util.fit.fminsearch_my({FunChi2, ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, ResLC.sigma_F_hat}, [1 2]);
%     Par_NLe(Isim,:) = Beta;
    
    
end

save Res_PL_EM_Alias1.mat Res ParLL InPar


%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 10;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end

    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;
end

save Res_PL_EM_Alias10.mat Res ParLL InPar
 
%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 5;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    
end

save Res_PL_EM_Alias5.mat Res ParLL InPar

%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 10;
InPar.Gamma = 2.5;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    
end

save Res_PL_EM_Alias10_gamma25.mat Res ParLL InPar

%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 10;
InPar.Gamma = 3.0;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    
end

save Res_PL_EM_Alias10_gamma30.mat Res ParLL InPar

%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 10;
InPar.Gamma = 3.5;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    
end

save Res_PL_EM_Alias10_gamma35.mat Res ParLL InPar


%%
SimName = 111;

Nsim = 1000;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0.02;
InPar.EndMatching = true;
InPar.AliasFactor = 10;
InPar.Gamma = 1.8;

Fun = @(Par,X) Par(1).*X.^(-Par(2));


ParLL      = zeros(Nsim,2);
%Par_Lin    = zeros(Nsim,2);
%Par_NLe    = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    try
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    catch
        [ResLC,ResG]=TimeDelay.rand_lensed(InPar);
    end
    
    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    
end

save Res_PL_EM_Alias10_gamma15.mat Res ParLL InPar


%%
load Res_PL_NEM_Alias1.mat
Res_PL_N1 = Res(1:Nsim);
ParPL_N1 = reshape([Res_PL_N1.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias1.mat
Res_PL_E1 = Res(1:Nsim);
ParPL_E1 = reshape([Res_PL_E1.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias10.mat
Res_PL_E10 = Res(1:Nsim);
ParPL_E10 = reshape([Res_PL_E10.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias5.mat
Res_PL_E5 = Res(1:Nsim);
ParPL_E5 = reshape([Res_PL_E5.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias10_gamma25.mat
Res_PL_E10_g25 = Res(1:Nsim);
ParPL_E10_g25 = reshape([Res_PL_E10_g25.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias10_gamma30.mat
Res_PL_E10_g30 = Res(1:Nsim);
ParPL_E10_g30 = reshape([Res_PL_E10_g30.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias10_gamma35.mat
Res_PL_E10_g35 = Res(1:Nsim);
ParPL_E10_g35 = reshape([Res_PL_E10_g35.BestPar_H1],[2 1000])';

load Res_PL_EM_Alias10_gamma15.mat
Res_PL_E10_g18 = Res(1:Nsim);
ParPL_E10_g18 = reshape([Res_PL_E10_g18.BestPar_H1],[2 1000])';

X = (1:0.1:3);
[N1] = histcounts(ParPL_N1(:,2),X);
[N2] = histcounts(ParPL_E1(:,2),X);
[N3] = histcounts(ParPL_E10(:,2),X);
[N4] = histcounts(ParPL_E5(:,2),X);
[N5] = histcounts(ParPL_E10_g25(:,2),X);
[N6] = histcounts(ParPL_E10_g30(:,2),X);
[N7] = histcounts(ParPL_E10_g35(:,2),X);
[N8] = histcounts(ParPL_E10_g18(:,2),X);

Xc = (X(1:end-1) + X(2:end)).*0.5;
plot(Xc,N1./sum(N1));
hold on;
plot(Xc,N2./sum(N2));
plot(Xc,N3./sum(N3));
plot(Xc,N4./sum(N4));
%plot(Xc,N5./sum(N5));
%plot(Xc,N6./sum(N6));
%plot(Xc,N7./sum(N7));
%plot(Xc,N8./sum(N8));

legend('No/Alias1','EM/Alias1','EM/Alias10','EM/Alias5')

%%
clear Mg;
Mg(1) = median(ParPL_E10(:,2));
Mg(2) = median(ParPL_E10_g25(:,2));
Mg(3) = median(ParPL_E10_g30(:,2));
Mg(4) = median(ParPL_E10_g35(:,2));
Mg(5) = median(ParPL_E10_g18(:,2));
Vg = [2 2.5 3 3.5 1.8];

plot(Vg,Mg,'o')





%%



%hist([ParLL(:,2),-ParRegular(:,2), -ParBin(:,2),-ParNL(:,2), -ParNLB(:,2)],(-1:0.1:3.5))
%hist([ParLL(:,2), -Par_Lin(:,1), Par_NLe(:,2)],(-1:0.1:3.5))
%hist([ParLL(:,2), Par_NLe(:,2)],(-1:0.1:4.5))
[N1,X]=histcounts(ParLL1(:,2),(1:0.02:3));
[N2,X]=histcounts(ParLL2(:,2),(1:0.02:3));
[N3,X]=histcounts(ParLL3(:,2),(1:0.02:3));
Xc = (X(1:end-1) + X(2:end)).*0.5;
plot(Xc,N1./sum(N1));
hold on;
plot(Xc,N2./sum(N2));
plot(Xc,N3./sum(N3));


axis([1 3 0 500])
H=legend('max$\mathcal{L}$','Lin. Fit','Binned Lin. Fit','Fit','Fit Binned');
H.Interpreter='latex';

H = xlabel('$\gamma$');
H.FontSize = 18;
H.Interpreter='latex';

H = ylabel('Number');
H.FontSize = 18;
H.Interpreter='latex';

%print PowerLawFit_EM_Alias10_DifferentMethods.eps -depsc2


%% error estimation

[ResLC,ResG]=TimeDelay.rand_lensed(InPar);

gammaVec = [1:0.01:4].';
Ngv = numel(gammaVec);

FitPar = [NaN     2        1./1000        1]
DefPar = [1       2        1./1000        1]

clear Res
for Igv=1:1:Ngv
    FitPar(2) = gammaVec(Igv);
    DefPar(2) = gammaVec(Igv);
    Res(Igv) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat,'DefPar',DefPar,'FitPar',FitPar);
end
plot(gammaVec,[Res.LL_H1]-[Res.LL_H0])
hold on;

GaussProb = 1-2.*normcdf([1:1:5],0,1,'upper');
Npar = 2;  % Tau, Alpha2
Level = 0.5.*chi2inv(GaussProb,Npar);

plot([1 3],Level(1).*ones(1,2),'k--');
plot([1 3],Level(2).*ones(1,2),'k--');
plot([1 3],Level(3).*ones(1,2),'k--');
axis([1 3 0 20]);

%%

end


%%  Simulation parameters
function InPar=select_parameters(SimName)
% generate parameters for specific simulation number
% Input  : - 11, 15

    InPar.Cyclic = false;
    InPar.x0  = 0;
    InPar.y0  = 0;
    InPar.y   = [0.0  0.0];    InPar.f_dc = 50;
    InPar.DeltaT  = 1;
    InPar.StdMeanRange = [0.1 0.15];
    InPar.AliasFactor  = 10;
    InPar.EndMatching  = true;
    
    switch SimName
        case 111
            % used
            InPar.Tau = 0;
            InPar.A0  = 0;
            InPar.A   = [1 0];
            InPar.x   = [0 0];
            InPar.Gamma = 2.0;
            InPar.TotTime = 1000;
            InPar.sigma_x = 0.02;
            InPar.sigma_F_rel = 0.02;
        

    end
end


%% Plot examples for simulated light curves
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
                                    'StdMeanRange',InPar.StdMeanRange,...
                                    'AliasFactor',InPar.AliasFactor,...
                                    'EndMatching',InPar.EndMatching);
end
