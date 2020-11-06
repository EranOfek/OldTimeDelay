function plots_for_paperII(Plot,Nsim)
% PLot options:
%       'single'
%     
% Example: TimeDelay.plots_for_paperII('single');

%%


switch lower(Plot)
    case 'single'
        
        %-------------------------------------------
        %% Plot examples for simulated light curves
        %-------------------------------------------
        InPar=select_parameters(1);
        InPar.EndMatching = false;
        % simulate LC
        [ResLC,ResG]=generateLC(InPar);

        %[ResLC.F_t,Slope]=TimeDelay.end_matching(ResLC.T,ResLC.F_t);
        %[ResLC.x_t,Slope]=TimeDelay.end_matching(ResLC.T,ResLC.x_t);
        %%

        figure(1)
        plot(ResLC.T,ResLC.F_t,'LineWidth',3)
        hold on;

        f_t=TimeDelay.reconstruct_ft(ResLC.F_t,ResLC.x_t,InPar.x0,InPar.x,InPar.A0,InPar.A(1),ResLC.sigma_F_hat,ResLC.sigma_x_hat);
        plot(f_t,'Color',[0.8 0.8 0.8],'LineWidth',3);

        plot(ResLC.f1_t,'LineWidth',2)
        plot(ResLC.f2_t,'LineWidth',2)



        legend('Combined','Reconstructed','Image 1','Image 2');
        H = xlabel('Time [days]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('Flux [arbitrary]');
        H.FontSize = 18;
        H.Interpreter = 'latex';

        %print Sim_3_LC.eps -depsc2
        %%



        figure(2)
        plot(ResLC.x_t,'LineWidth',2,'Color',[0.9 0.7 0.05]);
        hold on;
        plot(ResLC.chi_x_t,'k-','LineWidth',2);

        legend('x(t)','\chi(t)','combined');
        H = xlabel('Time [days]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('Position [arcsec]');
        H.FontSize = 18;
        H.Interpreter = 'latex';

        %print Sim_3_xt.eps -depsc2

        %%
        if 1==0
        figure(3);
        plot(ResLC.x_t,ResLC.y_t,'.','MarkerSize',8)
        H = xlabel('x(t) [arcsec]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('y(t) [arcsec]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        %print Sim5_xt_yt.eps -depsc2

        %Gx = ResLC.F_t.*ResLC.x_t;
        %Gy = ResLC.F_t.*ResLC.y_t;
        %figure(2);
        %plot(Gx,Gy,'.')

        figure(4);
        mGx = (ResLC.F_t-mean(ResLC.F_t)).*(ResLC.x_t-mean(ResLC.x_t));
        mGy = (ResLC.F_t-mean(ResLC.F_t)).*(ResLC.y_t-mean(ResLC.y_t));
        plot(mGx,mGy,'.')
        H = xlabel('G$_x$(t) [Flux arcsec]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('G$_y$(t) [Flux arcsec]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        %print Sim5_Gx_Gy.eps -depsc2

        %plot(ResLC.T,ResLC.x_t)
        %semilogy(ResLC.w,abs(ResLC.F_w))

        end
        
        %% Ghat/Fhat example

        figure(5);
        G = ResLC.F_t.*ResLC.x_t;
        Ghat = fft(G);
        Fhat = fft(ResLC.F_t);
        w    = TimeDelay.fft_freq(numel(ResLC.F_t),1);

        Ahat = (InPar.A(1).*InPar.x(1) + InPar.A(2).*InPar.x(2).*exp(1i.*ResLC.w.*InPar.Tau))./ ...
                (InPar.A(1) + InPar.A(2).*exp(1i.*ResLC.w.*InPar.Tau));

        [~,SI] = sort(w);
        
        semilogy(w(SI),abs(Ghat(SI)./Fhat(SI)));
        hold on
        semilogy(w(SI),abs(Ahat(SI)))
        %loglog(w(SI),abs(Ahat(SI)).*[w(SI).^(-InPar.Gamma.*0.5)]')
        axis([0 0.25 3e-2 3e1])

    case 'alpha'
        %-------------------------
        %% Plot examples alpha vs alpha ratio
        %-------------------------------------------
        InPar=select_parameters(1);
        InPar.EndMatching = false;
        % simulate LC
        [ResLC,ResG]=generateLC(InPar);

        
        FitPar = [InPar.A0   InPar.A(1)  InPar.A(2)  InPar.x0   InPar.x(1)   InPar.x(2)    InPar.Gamma];  % [A0, A1, A2, x0, x1, x2, gamma]
        VecA1 = logspace(log10(0.3),log10(10),10);
        VecA2dA1 = logspace(log10(0.1),log10(0.96),10);

        Res=TimeDelay.fit_scan_alpha_astrometric_flux(ResLC.T, ResLC.F_t, ResLC.x_t, ResLC.sigma_F_hat, ResLC.sigma_x_hat,...
                            'Tau',InPar.Tau,'FitPar',FitPar,'VecA1',VecA1,'VecA2dA1',VecA2dA1,'Min_w',[2.*pi./200]);


        %%
        figure(1);
        Data = Res.LL_xF;
        Min=min(Data(:));
        GaussProb = 1-2.*normcdf([1:1:5],0,1,'upper');
        Npar = 2;
        Level = 0.5.*chi2inv(GaussProb,Npar);
        contour(Res.A1,Res.A2dA1,Data'-Min,Level);
        H=colorbar;
        H.Label.String = '$\Delta$ ln[$\mathcal{L}$(x,F)]';
        H.Label.Interpreter = 'latex';
        H = xlabel('$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$\alpha_{2}$/$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        set(gca,'XS','log','YS','log')
        set(gca,'XTick',[0.1 0.3 1 3 10]);
        set(gca,'YTick',[0.05 0.1 0.3 0.5 0.9]);

        print Ast_Sim1_LLxF_A1vsA2A1.eps -depsc2
        %print Sim1_LLxF_A1vsA2A1_p_20_0_1_05_0_p01_m04_3.eps -depsc2

        %
        figure(2);
        Data = Res.LL_GF;
        Min=min(Data(:));
        GaussProb = 1-2.*normcdf([1:1:5],0,1,'upper');
        Npar = 2;
        Level = 0.5.*chi2inv(GaussProb,Npar);
        contour(Res.A1,Res.A2dA1,Data'-Min,Level);
        H=colorbar;
        H.Label.String = '$\Delta$ ln[$\mathcal{L}$(G$|$F)]';
        H.Label.Interpreter = 'latex';
        H = xlabel('$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$\alpha_{2}$/$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        set(gca,'XS','log','YS','log')
        set(gca,'XTick',[0.1 0.3 1 3 10]);
        set(gca,'YTick',[0.05 0.1 0.3 0.5 0.9]);

        %print Sim1_LLGF_A1vsA2A1_p_20_0_1_05_0_p01_m04_3.eps -depsc2

        %
        figure(3);
        Data = Res.LL_F;
        Min=min(Data(:));
        GaussProb = 1-2.*normcdf([1:1:5],0,1,'upper');
        Npar = 2;
        Level = 0.5.*chi2inv(GaussProb,Npar);
        contour(Res.A1,Res.A2dA1,Data'-Min,Level);
        H=colorbar;
        H.Label.String = '$\Delta$ ln[$\mathcal{L}$(F)]';
        H.Label.Interpreter = 'latex';
        H = xlabel('$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$\alpha_{2}$/$\alpha_{1}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        set(gca,'XS','log','YS','log')
        set(gca,'XTick',[0.1 0.3 1 3 10]);
        set(gca,'YTick',[0.05 0.1 0.3 0.5 0.9]);

        %print Sim1_LLF_A1vsA2A1_p_20_0_1_05_0_p01_m04_3.eps -depsc2


    case 'a0a1a2'

        %-----------------------------------
        %% Fit TimeDelay with unknown A0/A1/A2
        %-----------------------------------    
        InPar=select_parameters(1);
        InPar.EndMatching = false;
        InPar.AliasFactor = 1;

        VecInvTau = [(1./50:1./200:1./5)];
        VecInvTau = sort([VecInvTau, 1./20]);
        Nsim = 1

        FitPar = [NaN       NaN         NaN           InPar.x0   InPar.x(1)    InPar.x(2)    InPar.Gamma];
        DefPar = [InPar.A0  InPar.A(1)  InPar.A(2)    InPar.x0   InPar.x(1)    InPar.x(2)    InPar.Gamma];



        clear Res
        clear ResLC
        clear AllLC
        for Isim=1:1:Nsim
            [Isim, Nsim]
            
            % simulate LC
            [ResLC,ResG]=generateLC(InPar);

            tic;
            Res(Isim) = TimeDelay.fit_astrometric_flux(ResLC.T,ResLC.F_t,ResLC.x_t,ResLC.y_t,ResLC.sigma_F_hat,ResLC.sigma_x_hat,...
                                                 'TwoD',false,...
                                                 'DefPar',DefPar,...
                                                 'FitPar',FitPar,...
                                                 'Min_w',[2.*pi./200],...
                                                 'VecInvTau',VecInvTau);
            toc
            AllLC(Isim) = ResLC;                                 
        end

        %save -v7.3 Res.mat Res ResLCa

        % the fitted parameters distribution
        OutBestPar = zeros(Nsim,sum(isnan(FitPar)));
        DeltaL     = zeros(Nsim,1);
        BestTau    = zeros(Nsim,1);

        for Isim=1:1:Nsim
            [~,MinI] = min(abs(Res(Isim).Tau - InPar.Tau));
            [DeltaL(Isim),MinItau] = min(Res(Isim).LL_H1-Res(Isim).LL_H0);
            BestTau(Isim)  = Res(Isim).Tau(MinItau);
            OutBestPar(Isim,:) = Res(Isim).BestPar_H1(MinI,:);
        end

        %save -v7.3 Res_A012.mat Res OutBestPar AllLC

        for I=1:1:numel(Res)
            plot(1./Res(I).Tau,Res(I).LL_H1-Res(I).LL_H0)
            hold on;
        end
        %axis([0.02 0.2 -800 200]);

        H = xlabel('$1/\tau$ [day$^{-1}$]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$-\ln{\mathcal{L}(x,F,|\tau,\alpha_0,\alpha_1,\alpha_2)}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';

        %print Sim2_LogL_A0A1A2.eps -depsc2

        %%
        Flag=OutBestPar(:,1)>1;
        Flag=DeltaL(:,1)>-chi2inv(0.9973,2)./2;
        plot(OutBestPar(~Flag,2),OutBestPar(~Flag,3)./OutBestPar(~Flag,2),'.','MarkerSize',10)
        hold on
        plot(OutBestPar(Flag,2),OutBestPar(Flag,3)./OutBestPar(Flag,2),'.','MarkerSize',10)   

        H = xlabel('$\alpha_2$');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$\alpha_2/\alpha_1$');
        H.FontSize = 18;
        H.Interpreter = 'latex';

        print Sim2_LogL_A0A1A2_ParsA1A2.eps -depsc2

        %%

        H = plot(1./BestTau,DeltaL,'.');
        axis([0.02 0.21 -800 600]);
        H = xlabel('$1/\tau$ [day$^{-1}$]');
        H.FontSize = 18;
        H.Interpreter = 'latex';
        H = ylabel('$-\Delta{\ln{\mathcal{L}}}$');
        H.FontSize = 18;
        H.Interpreter = 'latex';


        [Hh]=plot.hist_ofplot('NbinX',20,'NbinY',20)

        print Sim2_LogL_vsBestTau.eps -depsc2


        
        

    otherwise
        error('Unknwon Plot option');
end






end % end main function


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
    InPar.AliasFactor  = 10;
    InPar.EndMatching  = true;
    
    switch SimName
        case 1
            % used
            InPar.Tau = 35;
            InPar.A0  = 0;
            InPar.A   = [1 0.5];
            InPar.x   = [0.2 -0.8];
            InPar.Gamma = 2.0;
            InPar.TotTime = 1000;
            InPar.sigma_x = 0.02;
            InPar.sigma_F_rel = 0.02.*sum(InPar.A);
            
        case 2
            % used
            InPar.Tau = 15;
            InPar.A0  = 0;
            InPar.A   = [1 0.5];
            InPar.x   = [0.1 -0.4];
            InPar.Gamma = 2.0;
            InPar.TotTime = 1000;
            InPar.sigma_x = 0.02;
            InPar.sigma_F_rel = 0.02.*sum(InPar.A);    
            
        case 5

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


function [ResLC,Section,PS]=generateLCuneq(InPar)

    %T=timeseries.random_time_sequence; 
    T=timeseries.random_time_sequence(6.*365,1,270,0.05,0.8);
    ResLC=TimeDelay.rand_lensed_uneq(T,InPar);
    
    N = numel(T);
    LC = [ResLC.T, ResLC.F_t, ResLC.Base.sigma_F_hat.*ones(N,1)];
    [PS,Section]=TimeDelay.power_spectrum_sections(LC);
    
end