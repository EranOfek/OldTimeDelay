function [MaxC,MaxAlpha,Alpha,CC,ResSim]=corr_rotinv(X,Y,Nsim)
% Rotationaly invariant Pearson correlation coef. including false alarm pr.
% Package: +TimeDelay
% Description: Given X and Y vectors, the Pearson correlation coef. is
%              calculated for each rotation of X,Y (default for (0:pi./100:pi).'
%              radians. Than the maximum correlation is selected and we
%              define this number as the Rotationaly invariant Pearson
%              correlation coef.
%              If more than 4 output arguments are requested than also
%              calculate the false alarm probablity using a permutation
%              test on the Y vector.
% Input  : - X vector.
%          - Y vector.
%          - Number of permutation simulations. Default is 100.
% Output : - The Rotationaly invariant Pearson correlation coef. defined
%            as the maximum correlation over all rotation angles in range 0
%            to pi.
%          - The rotation angle [radians] in which the max. correlation was
%            obtained.
%          - Vector of rotation angles tested.
%          - Vector of correlation coef. per rotation angle.
%          - A structure array with the permutation simulations results.
%            The following fields are available:
%            .FAProb - False alarm probability.
%            .VecMaxC - Vector of max corr. for each permutation.
%            .VecMaxA - Vector of max corr rot. angle for each permutation.
%      By : Eran O. Ofek                      Oct 2020
% Exanple: R=mvnrnd([0;0],[2 1; 1 1],100);
% [MaxC,MaxAlpha,Alpha,CC]=TimeDelay.corr_rotinv(R(:,1),R(:,2))
% [MaxC,MaxAlpha,~,~,ResSim]=TimeDelay.corr_rotinv(R(:,1),R(:,2))

if nargin<3
    Nsim = 1000;
end

Px  = sum(X.^2);
Py  = sum(Y.^2);
Cxy = sum(X.*Y);

Alpha = (0:pi./100:pi).';
Nalpha= numel(Alpha);
CC    = zeros(Nalpha,1);
CC = (0.5.*sin(2.*Alpha).*(Py-Px) + Cxy.*cos(2.*Alpha))./sqrt( (Px.*cos(Alpha).^2 + Py.*sin(Alpha).^2 + sin(2.*Alpha).*Cxy) .* (Px.*sin(Alpha).^2 + Py.*cos(Alpha).^2 - sin(2.*Alpha).*Cxy) );

[MaxC,MaxInd] = max(CC);
MaxAlpha = Alpha(MaxInd);

if nargout>4
    % simulations
    Npt = numel(X);
    VecMaxC = zeros(Nsim,1);
    VecMaxA = zeros(Nsim,1);
    for Isim=1:1:Nsim
        Ind = randperm(Npt);
        Ypn = Y(Ind);
        [VecMaxC(Isim),VecMaxA(Isim),~,~]=TimeDelay.corr_rotinv(X,Ypn);
    end
    ResSim.FAProb  = sum(VecMaxC>MaxC)./Nsim;
    ResSim.VecMaxC = VecMaxC;
    ResSim.VecMaxA = VecMaxA;
end
