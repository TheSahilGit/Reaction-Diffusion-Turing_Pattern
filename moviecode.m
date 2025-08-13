clear; clc; close all;

% Parameters (should match simulation)
L = 32;
Nx = 128; Ny = 128;
dx = L / Nx; dy = L / Ny;
x = linspace(0, L - dx, Nx);
y = linspace(0, L - dy, Ny);
[X, Y] = meshgrid(x, y);

% Frame saving
plot_every = 5000;        % Must match Fortran code
start_step = 1000;
end_step   = 5000000;

% Output video
v = VideoWriter('u_evolution.avi');
v.FrameRate = 15;
open(v);

% Fixed-size figure
fig = figure('Units','pixels','Position',[100 100 800 800]);
set(fig, 'Resize', 'off');


for step = start_step:plot_every:end_step
    filename = sprintf('data/u_%d.dat', step);
    
    if exist(filename, 'file') ~= 2
        warning('Missing file: %s. Skipping.', filename);
        continue;
    end
    
    u = load(filename);
    
    surf(X, Y, u, 'EdgeColor', 'none');
    shading interp;
    axis([0 L 0 L min(u(:)) max(u(:))]);
    axis square;
    view(2);
    colorbar;
    
    title(sprintf('u(x,y), t = %.3f', step * 1e-3), 'FontSize', 18);
    xlabel('x'); ylabel('y');
    
    drawnow;
    frame = getframe(fig);  % consistent size due to fixed window
    writeVideo(v, frame);
end

close(v);
fprintf('Movie saved.\n');
