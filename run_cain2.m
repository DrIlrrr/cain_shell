% one beam several events
clear all; close all; clc;
make_path
start_date=datestr(now);
%for new global
global home_dir DIRECTORY_FOR_CAIN BASE_DIRECTORY;
global rflags
[rflags] = flags_for_run;
home_dir=[pwd '/CAIN/'];

rflags.PLOTS =0;

for var_for_scan=[0]%0:10:100
    
BASE_DIRECTORY = [pwd '/test_1_MeV/'];

beam_phasespace=dlmread(['eli_lowen_oned_WP_newlayout_track_up_new_check_newsol_check.w5.asci'],'',35,0);

mkdir(BASE_DIRECTORY);
mkdir([BASE_DIRECTORY 'initial_beam_plot/']);
DIRECTORY_FOR_CAIN = [BASE_DIRECTORY 'cain_tmp/'];
mkdir(DIRECTORY_FOR_CAIN);

 beam_phasespace(:,6)=beam_phasespace(:,6)-114.7/0.511;
 
 [beam_phasespace] = defocusing_beam(beam_phasespace,1);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[beam_property]=formating_beam_for_cain(beam_phasespace,1);

for turn_number=1:1:10
    turn_number
    
    if turn_number==1
         [beam_property]=formating_beam_for_cain(beam_phasespace,1);
    end
   
    
    [nothing] = start_cain(beam_property,turn_number);
    
    
%     [beam_phasespace, X_RAY_photons_property]=start_cain_shell(beam_phasespace,turn_number);%;,factor_def);
    
%     X_RAY_phot_energy=X_RAY_photons_property(:,8);
%     
%     number_of_phot1=size(X_RAY_phot_energy);
%     number_of_scatered_photons(turn_number)=number_of_phot1(1);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end% end for turn_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% photons_plots(BASE_DIRECTORY,turn_number)
end

% check_flux
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
