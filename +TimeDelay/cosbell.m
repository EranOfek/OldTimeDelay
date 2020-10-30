function W=cosbell(T,FracRise)
% Cosine-bell weight function
% Package: +TimeDelay
% Input  : - Times.
%          - Fraction of rise.
% Output : - Cosine bell weight function
% Example: W=TimeDelay.cosbell((1:100),0.1);

MinT   = min(T);
MaxT   = max(T);
RangeT = range(T);
FracT  = FracRise.*RangeT;

% flag start and end points
FS = T<(MinT + FracT);
FE = T>(MaxT - FracT);

W  = ones(size(T));
W(FS) = sin((T(FS)-MinT).*0.5.*pi./FracT);
W(FE) = sin((MaxT-T(FE)).*0.5.*pi./FracT);









