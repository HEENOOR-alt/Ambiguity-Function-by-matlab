clc
clear
close all
%cw
tp = 100e-6;
fc = 10e9;
N = 90;
[t, fd] = meshgrid(linspace(-1, 1, N), linspace(-10, 10, 2*N));
A = abs(sin(pi*fd/tp*tp.*(1-abs(t*tp)/tp))/pi./fd/tp/tp./(1-abs(t*tp)/tp)).*(1-abs(t*tp)/tp);
A_norm = A/max(max(A));
% figure;surf(t, fd, A_norm, 'EdgeColor', 'none');
figure;
mesh(t, fd, A_norm);
colormap jet;
colorbar;
title('CW signal');
xlabel({'Delay(Second)'});
ylabel({'Doppler Shift(Hertz)'});
zlabel('Amplitude');

figure;
imagesc(A_norm);
colormap jet;
colorbar;
set(gca, 'XTick', [1 N/2 N]);
set(gca, 'XTickLabel', {-1 0 1});
set(gca, 'YTick', [1 N N*2]);
set(gca, 'YTickLabel', {-10 0 10});
title('CW signal');
xlabel('Delay(Second)');
ylabel('Doppler Shift(Hertz)');

figure;
contour(t, fd, A_norm);
grid on
colormap jet;
colorbar;
title('CW signal');
xlabel('Delay(Second)');
ylabel('Doppler Shift(Hertz)');

[row, col] = size(A_norm);
range_ambiguity = A_norm(row/2, :);
velocity_ambiguity = A_norm(:, col/2);

figure;
plot(t, range_ambiguity, 'b');
grid on;
title('Range Ambiguity');
xlabel('Delay');
ylabel('|\chi(\tau,0)|');

figure;
plot(fd, velocity_ambiguity, 'b');
grid on;
title('Velocity Ambiguity');
xlabel('Doppler Shift');
ylabel('|\chi(0,f_d)|');
