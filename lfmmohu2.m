%% 线性调频信号
clc
clear 
close all

tp = 10e-6;
B = 1e6;
fc = 10e6;
k = B/tp;
N = 90;
[t, fd] = meshgrid(linspace(-1, 1,N), linspace(-20, 20, 2*N));
A = (1-abs(t*tp)/tp).*abs(sin(pi*(fd/tp*tp+B*t*tp).*(1-abs(t*tp)/tp))/pi./(fd/tp*tp+B*t*tp)./(1-abs(t*tp)/tp));
A_norm = A/max(max(A));

figure;
mesh(t, fd, A_norm);
colormap jet;
grid on
colorbar;
title('LFM');
xlabel({'Delay'});
ylabel({'Doppler Shift'});
zlabel('Amplitude');

figure;
imagesc(A_norm);
colormap jet;
colorbar;
grid on
set(gca, 'XTick', [1 N/2 N]);
set(gca, 'XTickLabel', {-1 0 1});
set(gca, 'YTick', [1 N N*2]);
set(gca, 'YTickLabel', {-20 0 20});
title('LFM');
xlabel('Delay');
ylabel('Doppler Shift');
	
figure;
contour(t, fd, A_norm);
colormap jet;
title('LFM');
xlabel('Delay');
ylabel('Doppler Shift');
colorbar;
grid on
	
[row, col] = size(A_norm);
range_ambiguity = A_norm(row/2, :);
velocity_ambiguity = A_norm(:, col/2);
figure;plot(t, range_ambiguity, 'b');grid on;
title('Range Ambiguity');
xlabel('Delay');
ylabel('|\chi(\tau,0)|');

figure;
plot(fd, velocity_ambiguity, 'b');
grid on;
title('Velocity Ambiguity');
xlabel('Doppler Shift');
ylabel('|\chi(0,f_d)|');
