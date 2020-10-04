function FreqVec=fftsh_freq(N)
% Return the frequencies corresponding to fftshift(fft(vec_of_size_N))
% Package: +TimeDelay
% Description: Return the frequencies corresponding to
%              fftshift(fft(vec_of_size_N)), without dviding by the total
%              time span.
%              A version of this function also exist in +Util/+fft/
% Input  : - Number of points.
% Output : - Vector of frequencies.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FreqVec=TimeDelay.fftsh_freq(4)
% Reliable: 2
%--------------------------------------------------------------------------


if (N.*0.5)==floor(0.5.*N)
    % even
    FreqVec = (-N.*0.5:1:N.*0.5-1).';
else
    FreqVec = (-N.*0.5+0.5:1:N.*0.5-0.5).';
end

