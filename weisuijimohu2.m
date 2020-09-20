%% 伪随机码（PRN）/最大长度序列（MLS）码
%Compute the ambiguity function utilizing the FFT
%through combining multiple range cuts.
clc
clear
close all
clear u ambig
N = 31;	%15, 31

switch N
	case 15
		PRN = [1 -1 -1 -1 1 1 1 1 -1 1 -1 1 1 -1 -1];
	case 31
		PRN = [1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 -1 1 1 ...
				-1 -1 -1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 -1];
	otherwise
		disp 'There doesn''t exist such PRN code.';
	return;
end
tau = N;
samp_num = N*10;
nfft = 2^ceil(log(samp_num)/log(2));
u(1:nfft) = 0;
j = 0;

for index = 1:10:samp_num
	j = j+1;
	u(index:index+10-1) = PRN(j);
end

if 1
	delay = linspace(0, 10*tau, nfft);
else
	delay = linspace(-tau, tau, nfft);
end
	freq_del = 8/tau/100;
	j = 0;
	vfft = fft(u, nfft);
for freq = -4/tau:freq_del:4/tau
	j = j+1;
	exf = exp(1i*2*pi*freq*delay);
	u_times_exf = u.*exf;
	ufft = fft(u_times_exf, nfft);
	prod = ufft.*conj(vfft);
	ambig(:, j) = fftshift(abs(ifft(prod))');
end
freq = -4/tau:freq_del:4/tau;
delay = linspace(-N, N, nfft);

figure;
mesh(freq, delay, ambig/max(max(ambig)));
colormap jet;
grid on
colorbar;
axis tight;
title('31-bit M sequence code');
xlabel('Doppler Shift');
ylabel('Delay');
zlabel('Amplitude');
% 		xlabel('frequency');
% 		ylabel('delay');
% 		zlabel('ambiguity function a PRN code');

figure;
contour(freq, delay, ambig/max(max(ambig)));
colormap jet;
colorbar;
grid on
title('31-bit M sequence code');
xlabel('Doppler Shift');
ylabel('Delay');
% 		xlabel('frequency');
% 		ylabel('delay');

[row, col] = size(ambig);
range_ambiguity = ambig(:, ceil(col/2));
velocity_ambiguity = ambig(ceil(row/2), :);

figure;
plot(delay, range_ambiguity/max(max(ambig)), 'b');
axis tight;
grid on;
title('Range Ambiguity');
xlabel('Delay');
ylabel('|\chi(\tau,0)|');
% 		xlabel('delay');
% 		ylabel('normalized amibiguity cut for f=0');

figure;
plot(freq, velocity_ambiguity/max(max(ambig)), 'b');
grid on;
axis tight;
title('Velocity Ambiguity');
xlabel('Doppler Shift');
ylabel('|\chi(0,f_d)|');
