%% inputs and declarations
clc;
X=512;
Y=512;
c1=1500;
c2=2500;
rho1=1000;
rho2=600;
source_freq=2.5e6;
angle=15;
%% resizing and cropping
sound_Big=[c1*ones(X,Y/2),c2*ones(X,Y/2)];
sound_Big=imrotate(sound_Big,angle);
Nx=X/2;
Ny=Y/2;
% sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,Y/2-Ny/2+1:Y/2+Ny/2);
dx=200e-6;
dy=200e-6;
kgrid=kWaveGrid(Nx,dx,Ny,dy);
density_Big=[rho1*ones(X,Y/2),rho2*ones(X,Y/2)];
density_Big=imrotate(density_Big,angle);
% density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,Y/2-Ny/2+1:Y/2+Ny/2);
medium.sound_speed=sound_crop;
medium.density=density_crop;
kgrid.makeTime(medium.sound_speed);
%% define a time varying sinusoidal source
source_mag = 5;                             % [Pa]
source.p0 = zeros(Nx,Ny);
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
source.p_mask = makeLine(Nx, Ny,[Nx/2-Nx/32+1,1],[Nx/2+Nx/32,1]);
%% wavelength,aperture width and fraunhoff zone
WL=c1/source_freq;
D=length(source.p_mask(:,1)==1);
r=D^2*200e-3/4*WL;
%% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);
%% create a display mask to display the transducer
display_mask = source.p_mask;
sensor.mask = [1, 1, Nx, Ny].';
%% set the record mode capture the final wave-field and the statistics at
% each sensor point 
sensor.record = {'p_final', 'p_max', 'p_rms'};
% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};
%% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor,input_args{:});
% =========================================================================
% VISUALISATION
% =========================================================================
%% add the source mask onto the recorded wave-field
sensor_data.p_final(source.p_mask ~= 0) = 1;
sensor_data.p_max(source.p_mask ~= 0) = 1;
sensor_data.p_rms(source.p_mask ~= 0) = 1;

%% plot the final wave-field
figure;
subplot(1, 3, 1);
imagesc(kgrid.y_vec *10 , kgrid.x_vec*10, sensor_data.p_final, [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Final Wave Field');

%% plot the maximum recorded pressure
subplot(1, 3, 2);
imagesc(kgrid.y_vec*10 , kgrid.x_vec*10 , sensor_data.p_max,[-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Maximum Pressure');

%% plot the rms recorded pressure
subplot(1, 3, 3);
imagesc(kgrid.y_vec*10, kgrid.x_vec*10 , sensor_data.p_rms, [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('RMS Pressure');
scaleFig(2, 2);