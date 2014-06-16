clear all; close all; clc;

out_folder = [pwd '/new_beam_for_candarelli_1_MeV/'];
for turn_number=[20000]
mkdir([ out_folder 'photon_plots_' num2str(turn_number)])

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

all_data_phot=[];
full_photons=[];

    full_photons_1=[];
    full_photons_2=[];
    full_photons_3=[];
    full_photons_4=[];
    full_photons_5=[];
    full_photons_6=[];
    full_photons_7=[];
    full_photons_8=[];
    full_photons_9=[];
    full_photons_10=[];
    full_photons_11=[];
    full_photons_12=[];
    full_photons_13=[];
    full_photons_14=[];
for ni=1:1:turn_number
    ni
    photons_data=dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
    
    
  
    full_photons_1=[full_photons_1;photons_data(:,1)];
    full_photons_2=[full_photons_2;photons_data(:,2)];
    full_photons_3=[full_photons_3;photons_data(:,3)];
    full_photons_4=[full_photons_4;photons_data(:,4)];
    full_photons_5=[full_photons_5;photons_data(:,5)];
    full_photons_6=[full_photons_6;photons_data(:,6)];
    full_photons_7=[full_photons_7;photons_data(:,7)];
    full_photons_8=[full_photons_8;photons_data(:,8)];
    full_photons_9=[full_photons_9;photons_data(:,9)];
    full_photons_10=[full_photons_10;photons_data(:,10)];
    full_photons_11=[full_photons_11;photons_data(:,11)];
    full_photons_12=[full_photons_12;photons_data(:,12)];
    full_photons_13=[full_photons_13;photons_data(:,13)];
    full_photons_14=[full_photons_14;photons_data(:,14)];
    
end
full_photons=[full_photons_1';full_photons_2';full_photons_3'./turn_number;full_photons_4';full_photons_5';full_photons_6';...
    full_photons_7';full_photons_8';full_photons_9';full_photons_10';full_photons_11';full_photons_12';full_photons_13';full_photons_14';];

fid = fopen([out_folder 'photons' num2str(turn_number) '.dat'],'w');%save beam for cain standart
fprintf(fid,' %i    %i       %1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e % 1.12e \n',full_photons);
fclose(fid);
% 
% 
 xx=load([out_folder 'photons' num2str(turn_number) '.dat']);
% 
figure(1)
 hist(xx(:,8))

end