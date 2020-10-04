function Freq=fft_freq(N, D)
% Return the frequencies off fft with N points and Delta time step D.
% Package: 
% Input  : - Number of data points
%          - Delta time.
% Output : - Vector of frequencies.
% Example: Freq=fft_freq(5, 1);


if mod(N, 2) == 0
    Freq = [0:(N/2-1), -N/2:-1] / (D.*N);
else
    Freq = [0:(N-1)/2, -(N-1)/2:-1] / (D.*N);
end