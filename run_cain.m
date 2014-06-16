
clear all; close all; clc;
make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];
mkdir(home_dir)

for aa=[0]%0:10:100
rflags.shifting_laser_t=0;%1e-12*3e8*aa;

BASE_DIRECTORY = [pwd '/new_beam/'];

[laser_parameters] = cain_parameteres;
beam_parameters=laser_parameters;

%  beam_load=load([pwd '/initial_beam_data.dat'],'-mat');
%  beam_phasespace_1=beam_load.beam_phasespace;
beam_phasespace=dlmread(['eli_360MeV_double_A_d.out.asci'],'',10,0);


mkdir(BASE_DIRECTORY);
mkdir([BASE_DIRECTORY 'initial_beam_plot/']);

% DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
mkdir(DIRECTORY_FOR_CAIN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for turn_number=1:1:2
    turn_number
    
    
    [beam_phasespace, X_RAY_photons_property]=start_cain_shell(beam_phasespace,turn_number);%;,factor_def);
    
%     X_RAY_phot_energy=X_RAY_photons_property(:,8);
%     
%     number_of_phot1=size(X_RAY_phot_energy);
%     number_of_scatered_photons(turn_number)=number_of_phot1(1);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end% end for turn_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% photons_plots(BASE_DIRECTORY,turn_number)
end
%%%%%%%%%%%%%%%%%%%%Save all sum photons data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(COMPTON_BACK_SCATERING==1)
%     % if(SAVE_all_sum_X_RAY_photons_property_DATA==1)%
%     for turn_number=(1+frozen_turn_number):1:beam_parameters.NUMBER_OF_TURNS
%         x_ray_prop=load([DIRECTORY_FOR_CAIN 'cain_output_photons_' num2str(turn_number) '.dat']);
%         %photon_spectrum=x_ray_prop(:,ENERGY_OF_PARTICLE);
%         if(SAVE_photons_turn_by_turn==0)
%             system(['rm ' DIRECTORY_FOR_CAIN 'cain_output_photons_' num2str(turn_number) '.dat']); % delete cain_output_photons_' num2str(turn_number) '.dat
%         end
%         if(SAVE_all_sum_X_RAY_photons_property_DATA==1)%
%             all_sum_X_RAY_photons_property=[all_sum_X_RAY_photons_property;x_ray_prop];
%         end
%     end
%
%     if(SAVE_all_sum_X_RAY_photons_property_DATA==1)%
%         %         save([BASE_DIRECTORY 'ALL_SUM_X_RAY_data_' num2str(beam_parameters.NUMBER_OF_TURNS,3) '_' num2str(beam_parameters.NUMBER_OF_MACROPARTICLES) '_' num2str(beam_parameters.pulseE) '.dat'],'all_sum_X_RAY_photons_property');
%         save([BASE_DIRECTORY 'ALL_SUM_X_RAY_data_angle_' num2str(beam_parameters.angle) '.dat'],'all_sum_X_RAY_photons_property');
%     end
%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
