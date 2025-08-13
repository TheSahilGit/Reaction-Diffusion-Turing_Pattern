clear; clc; close all;

% Parameters (match those used in Fortran)
L = 2*pi;
Nx = 128; Ny = 128;
dx = L / Nx; dy = L / Ny;

x = linspace(0, L - dx, Nx);
y = linspace(0, L - dy, Ny);
[X, Y] = meshgrid(x, y);

% Define time steps to load (adjust file names if needed)
steps = [1000, 100000, 200000, 300000, 400000, 500000];

% Plot settings
figure('Position',[100 100 1400 800]);

for k = 1:length(steps)
    step = steps(k);
    filename = sprintf('data/u_%d.dat', step);
    u = load(filename);

    subplot(2,3,k)
    surf(X, Y, u, 'EdgeColor', 'none');
    title(num2str(step));
    xlabel('x'); ylabel('y'); zlabel('u');
    axis square; view(2); colorbar;
    xticks([0 pi 2*pi])
    xticklabels([0 "\pi" "2\pi"])
    set(gca, "FontSize", 26)
end

sgtitle('Snapshots of u(x,y) at Different Time Steps', 'FontSize', 22);
