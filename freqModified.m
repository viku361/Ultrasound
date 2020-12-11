clc;
c1=input('Enter velocity of sound in first medium');
c2=input('Enter velocity of sound in second medium');
rho1=input('Enter density of sound in first medium');
rho2=input('Enter density of sound in second medium');
source_freq=input('Enter center frequency in MHz');
angle=input('Enter angle of orientation');
resolution=input('Enter resolution in powers of e-6');
D=input('Enter transducer diameter in mms');
interface=input('Type 1 for straight or 2 for spherical');
% option=input('Type 1 for near field or 2 for far field');
%% resizing and cropping
grid_source=D*1000/resolution;
dx=resolution*10^-6;
dy=resolution*10^-6;
WL=c1/(source_freq*10^6);
N=(D^2)*10^-6/(4*WL);
NearFieldgrid=N/dx;
formatSp = 'near field is %4.2f \n';
fprintf(formatSp,NearFieldgrid);
if interface==1
    if option ==1
        grid_near=round(2*N*10^6/resolution,-2);
        sound_Big=[c1*ones(grid_near,grid_near/2),c2*ones(grid_near,grid_near/2)];
        sound_Big=imrotate(sound_Big,angle);
        Nx=grid_near/2;
        Ny=grid_near/2;
        % sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        sound_crop= sound_Big(grid_near/2-Nx/2+1:grid_near/2+Nx/2,grid_near/2-Ny/2+1:grid_near/2+Ny/2);
        kgrid=kWaveGrid(Nx,dx,Ny,dy);
        density_Big=[rho1*ones(grid_near,grid_near/2),rho2*ones(grid_near,grid_near/2)];
        density_Big=imrotate(density_Big,angle);
        % density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        density_crop= density_Big(grid_near/2-Nx/2+1:grid_near/2+Nx/2,grid_near/2-Ny/2+1:grid_near/2+Ny/2);
        medium.sound_speed=sound_crop;
        medium.density=density_crop;
        kgrid.makeTime(medium.sound_speed);
    elseif option==2
        grid_far=round(5*N*10^6/resolution,-2);
        sound_Big=[c1*ones(grid_far,grid_far/2),c2*ones(grid_far,grid_far/2)];
        sound_Big=imrotate(sound_Big,angle);
        Nx=grid_far/2;
        Ny=grid_far/2;
        kgrid=kWaveGrid(Nx,dx,Ny,dy);
        % sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        sound_crop= sound_Big(grid_far/2-Nx/2+1:grid_far/2+Nx/2,grid_far/2-Ny/2+1:grid_far/2+Ny/2);
        density_Big=[rho1*ones(grid_far,grid_far/2),rho2*ones(grid_far,grid_far/2)];
        density_Big=imrotate(density_Big,angle);
        % density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        density_crop= density_Big(grid_far/2-Nx/2+1:grid_far/2+Nx/2,grid_far/2-Ny/2+1:grid_far/2+Ny/2);
        medium.sound_speed=sound_crop;
        medium.density=density_crop;
        kgrid.makeTime(medium.sound_speed);
    else
        print('You have typed wrong button');
    end
else 
    if option==1
%       grid_near=round(2*N*10^6/resolution,-2);
        scan_area=input('What will be the depth you want to scan in terms of mm');
        grid_size=scan_area*1000/resolution;
        formatSpec = 'Forming grid %4.2f X %8.3f \n';
        X=grid_size;
        Y=2*grid_size;
        fprintf(formatSpec,X,Y);
        fx=input('Enter focal point(x coordinate) in terms of grid points ');
        fy=input('Enter focal point(y coordinate) in terms of grid_points');
        rad_in=input('Enter radius in terms of X grid not greater than 1');
        radius=rad_in*X;
        sound_Big = c2*makeDisc(X,Y,fx,fy,radius);
        sound_Big=imrotate(sound_Big,angle);
        sound_Big(sound_Big==0)=c1;
        Nx=X/2;
        Ny=Y/2;
        % sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        sound_crop= sound_Big(X/2-Nx/2+1:X/2+Nx/2,Y/2-Ny/2+1:Y/2+Ny/2);
        kgrid=kWaveGrid(Nx,dx,Ny,dy);
        density_Big = rho2*makeDisc(X,Y,fx,fy,radius);
        density_Big=imrotate(density_Big,angle);
        density_Big(density_Big==0)=rho1;
        % density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,5*Y/8-Ny/2+1:5*Y/8+Ny/2);
        density_crop= density_Big(X/2-Nx/2+1:X/2+Nx/2,Y/2-Ny/2+1:Y/2+Ny/2);
        medium.sound_speed=sound_crop;
        medium.density=density_crop;
        kgrid.makeTime(medium.sound_speed);
    else
        print('You have typed wrong button');
    end
end  
%% define a time varying sinusoidal source
source_mag = 1;                             % [Pa]
source.p0 = zeros(Nx,Ny);
source.p = source_mag * sin(2 * pi * source_freq * 10^6 * kgrid.t_array);
source.p_mask = makeLine(Nx, Ny,[Nx/2-grid_source/2+1,1],[Nx/2+grid_source/2,1]);
% %% wavelength,aperture width and fraunhoff zone
% WL=c1/source_freq;
% D=length(source.p_mask(:,1)==1)*dx;
% r=D^2/4*WL;
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



