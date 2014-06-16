
function [rflags] = flags_for_run
% create global flags to run choosen type of simulations
rflags=struct('NAME_OF_FLAGS',1);



% rflags.NAME_OF_PARAMS = ['default']; %Number of turns
YES=1; yes=1; NO=0; no=0;
rflags.name_of_job=[''];
rflags.PLOTS =yes;


% LASER parameteres
rflags.angle_deg=7.5;% in degree
rflags.angle=rflags.angle_deg*(pi/180); %initial scattered angle [rad]
rflags.pulseE=0.2; %laser puse energy [J]
rflags.sigLr=28/2; % given in [mu m] micro meter like 2 weist w0=28;
rflags.laserwl=515; % laser wavelenth [nm] nano meters
rflags.sigt=1.5; %pulse length [ps]
rflags.shifting_laser_x = 0;  %
rflags.shifting_laser_y = 0;  %
rflags.shifting_laser_s = 0;  %
rflags.shifting_laser_t = 0;%rflags.shifting_laser_t;  %