% funciton to calculate the continuous fourier transform using the FFT

% must assume that the function is centered around 0

% define this in a different way
function [ft xm] = ctfft(f,t)

m = length(t);

a = t(end)*2;
beta = a/m;
xm = (t/beta)*2*pi/(m*beta);

ft = ((-1).^[1:m]).*beta.*fft((-1).^[1:m].*f,[],2);
return