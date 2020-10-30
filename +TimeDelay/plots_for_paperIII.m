function plots_for_paperIII
%% single light curve example
SimName = 111;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0;
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






[Par,ParErr,Chi2,Dof] = Util.fit.fitpow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1);
loglog(ResLC.w(Flag), Par(1).*ResLC.w(Flag).^(-Par(2)))

% fit using non-lin chi^2
[ResNL]=Util.fit.fit_pow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1,'Par0',[1 -2]);
loglog(ResLC.w(Flag), ResNL.Par(1).*ResLC.w(Flag).^(ResNL.Par(2)))

%%
SimName = 111;

Nsim = 100;

InPar=select_parameters(SimName);
InPar.sigma_F_rel = 0;
InPar.EndMatching = false;
InPar.AliasFactor = 1;


ParLL      = zeros(Nsim,2);
Par_Lin    = zeros(Nsim,2);

ParRegular = zeros(Nsim,2);
ParBin     = zeros(Nsim,2);

ParNL      = zeros(Nsim,2);
ParNLB     = zeros(Nsim,2);

for Isim=1:1:Nsim
    Isim
    
    [ResLC,ResG]=TimeDelay.rand_lensed(InPar);

    ParPoly = polyfit(ResLC.T,ResLC.F_t,1);
    ResLC.F_t = ResLC.F_t - polyval(ParPoly,ResLC.T);
    ResLC.F_w = fft(ResLC.F_t) ./ sqrt(numel(ResLC.F_t));
    
    
    % max like fit
    Res(Isim) = TimeDelay.fit_fluxBPL(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat);
    ParLL(Isim,:) = [Res(Isim).BestPar_H0];
    
    Flag = ResLC.w>0;

    % linear fit to all data points
    X = log10(ResLC.w(Flag));
    Y = log10( abs(ResLC.F_w(Flag)).^2);
    Par_Lin(Isim,:) = polyfit(X,Y,1);
    
    
    % linear fot to binned data points
    BinVec = logspace(log10(min(ResLC.w(Flag))), log10(max(ResLC.w(Flag))),30);
    B      = timeseries.binning([ResLC.w(Flag), ResLC.F_w(Flag)], BinVec,[NaN NaN],{'MeanBin',@mean,@std,@numel});
    
    
    
%     % linearized fit to all points
%     Flag = (ResLC.w)>0;
%     [Par,ParErr,Chi2,Dof] = Util.fit.fitpow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1);
%     ParRegular(Isim,:) = Par(:).';
%     
%     % linaerized fit to binned points
%     Bin_w = logspace(log10(1./1000),log10(pi),10);
%     Nbin = numel(Bin_w);
%     Mean = zeros(Nbin-1,2);
%     for Ibin=1:1:Nbin-1
%         Flag = ResLC.w>Bin_w(Ibin) & ResLC.w<=Bin_w(Ibin+1);
%         Mean(Ibin,1)   = mean(abs(ResLC.F_w(Flag)).^2);
%         Mean(Ibin,2) = mean(ResLC.w(Flag));
%     end
%     FF = ~isnan(sum(Mean,2));
%     [Par,ParErr,Chi2,Dof] = Util.fit.fitpow(Mean(FF,1),Mean(FF,2),1);
%     ParBin(Isim,:) = Par(:).';
%     
%     % fit using non-lin chi^2
%     [ResNL]=Util.fit.fit_pow(ResLC.w(Flag), abs(ResLC.F_w(Flag)).^2, 1,'Par0',[1 -2]);
%     ParNL(Isim,:) = ResNL.Par;
%     
%     % fit using non-lin chi^2 on binned data
%     [ResNLB]=Util.fit.fit_pow(Mean(FF,1), Mean(FF,2), 1,'Par0',[1 -2]);
%     ParNLB(Isim,:) = ResNLB.Par;
    
end



%%

%hist([ParLL(:,2),-ParRegular(:,2), -ParBin(:,2),-ParNL(:,2), -ParNLB(:,2)],(-1:0.1:3.5))
hist([ParLL(:,2), -Par_Lin(:,1)],(-1:0.1:3.5))
axis([-0.5 2.5 0 1000])
H=legend('max$\mathcal{L}$','Lin. Fit','Binned Lin. Fit','Fit','Fit Binned');
H.Interpreter='latex';

H = xlabel('$\gamma$');
H.FontSize = 18;
H.Interpreter='latex';

H = ylabel('Number');
H.FontSize = 18;
H.Interpreter='latex';

print PowerLawFit_EM_Alias10_DifferentMethods.eps -depsc2


'a'


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
