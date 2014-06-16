function photons_plots(out_folder,turn_number)

close all;

%out_folder=['new_5000_const_em_1'];

mkdir([ out_folder 'photon_plots'])

%  1  2         3     4    5    6    7     8      9        10       11    12 13 14
%  K GEN NAME Weight T(m) X(m) Y(m) S(m) E(eV) Px(eV/c) Py(eV/c) Ps(eV/c) Sx Sy Ss

all_data_phot=[];
full_spectrum=[];
phot_angle=[];
weigth=0;
histmat=[];
new_m=[];
x_phot=[];
y_phot=[];
xp_phot=[];
yp_phot=[];

for ni=1:1:turn_number
    photons_data=dlmread([out_folder 'cain_tmp/cain_output_photons_' num2str(ni) '.dat'],'',1,0);%read electrons from cain
    
    x_phot=[x_phot;photons_data(:,5)];
    y_phot=[y_phot;photons_data(:,6)];
    xp_phot=[xp_phot;photons_data(:,9)];
    yp_phot=[yp_phot;photons_data(:,10)];
    
    weigth=photons_data(1,3);
%     phot_number=size(photons_data,1)*weigth;
%     zx=size(photons_data(:,8));
%     number_of_photons=zx(1)*weigth;
    
    
    full_spectrum=[full_spectrum;photons_data(:,8)./1e3];
    % phot_angle=[phot_angle;sign(photons_data(:,9)).*atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
    phot_angle=[phot_angle;atan(sqrt(photons_data(:,9).^2+photons_data(:,10).^2)./photons_data(:,11))];
    
end

number_of_photons=length(full_spectrum)*weigth*weigth;
%SPEED_OF_LIGHT=3e8;
el_angel=4e-5
aa=find(abs(phot_angle)<el_angel);

ifig=1;

bandwith_cm=[];
num_phot_th=[];
for ni=1:1:100
el_angel=ni*0.5e-6;
aa_1=find(abs(phot_angle)<el_angel);

num_phot_th(ni)=length(full_spectrum(aa_1))*weigth/5000;
bandwith_cm(ni)=std(full_spectrum(aa_1))/mean(full_spectrum(aa_1));

end

figure(ifig)
ifig=ifig+1;
subplot(2,1,1)
plot((1:1:100).*0.5e-6,num_phot_th./turn_number,'-xb')
grid on
ylabel('number photons (macro)')
subplot(2,1,2)
plot((1:1:100).*0.5e-6,bandwith_cm,'--sb')
ylabel('bandwith')
grid on
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);







figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta [rad]')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum,50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
 xlabel('photons energy (KeV)')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]})
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



figure(ifig)
ifig=ifig+1;
subplot(1,2,1)
hold on
hist(phot_angle(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('theta')
% xlabel('photons energy (KeV)')
subplot(1,2,2)
hold on
hist(full_spectrum(aa),50)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
%xlabel('theta')
 xlabel('photons energy (KeV)')
 suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



nbin_plot=30;
figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(0,max(full_spectrum),nbin_plot),smooth(hist(full_spectrum,nbin_plot)*weigth),'-r','LineWidth',2)
hold off
grid on
ylim([0 max(smooth(hist(full_spectrum,nbin_plot)*weigth))])
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')

title({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]})
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);




figure(ifig)
ifig=ifig+1;
hold on
plot(phot_angle,full_spectrum,'.b')
hold off
grid on
set(gca,'FontSize',16)
ylabel('energy of scattered photons')
xlabel('scattered angle')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]});
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);






figure(ifig)
ifig=ifig+1;
hold on
plot(phot_angle(aa),full_spectrum(aa),'.b')
hold off
grid on
set(gca,'FontSize',16)
ylabel('energy of scattered photons')
xlabel('scattered angle')
suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);

figure(ifig)
ifig=ifig+1;
hold on
plot(linspace(0,max(full_spectrum)),histc(full_spectrum(aa),linspace(0,max(full_spectrum)))*weigth,'-r','LineWidth',2)
% plot([30:0.5:47],histc(full_spectrum(aa),[30:0.5:47])*weigth,'-r','LineWidth',2)
hold off
grid on
set(gca,'FontSize',16)
ylabel('number of scattered photons')
xlabel('photons energy (KeV)')
% suptitle(['For theta<' num2str(el_angel) ' [rad]']);
suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})
text(93.1658986175114, 110886076.949367,{['bandwidth=' num2str(std(full_spectrum(aa))/mean(full_spectrum(aa)))];...
    ['FWHM=' num2str(0)]})
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);









ifig=ifig+1;
figure(ifig)
subplot 221
set(gca,'FontSize',16)
plot(x_phot,y_phot,'.b')
xlabel('x [m]')
ylabel('y [m]')
subplot 222
set(gca,'FontSize',16)
plot(x_phot,xp_phot,'.b')
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
set(gca,'FontSize',16)
plot(y_phot,yp_phot,'.b')
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
plot(xp_phot,yp_phot,'.b')
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')
suptitle({['For all theta'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ]});
filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);



ifig=ifig+1;
figure(ifig)

subplot 221
set(gca,'FontSize',16)
plot(x_phot(aa),y_phot(aa),'.b')
xlabel('x [m]')
ylabel('y [m]')
subplot 222
set(gca,'FontSize',16)
plot(x_phot(aa),xp_phot(aa),'.b')
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
set(gca,'FontSize',16)
plot(y_phot(aa),yp_phot(aa),'.b')
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
plot(xp_phot(aa),yp_phot(aa),'.b')
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')

suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})


filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);








ifig=ifig+1;
figure(ifig)
subplot 221

set(gca,'FontSize',16)
hold on
plot(x_phot,y_phot,'.b')
plot(x_phot(aa),y_phot(aa),'.r')
hold off
xlabel('x [m]')
ylabel('y [m]')
subplot 222
hold on
set(gca,'FontSize',16)
plot(x_phot,xp_phot,'.b')
plot(x_phot(aa),xp_phot(aa),'.r')
hold off
xlabel('x [m]')
ylabel('Px [eV/c]')
subplot 223
hold on
set(gca,'FontSize',16)
plot(y_phot,yp_phot,'.b')
plot(y_phot(aa),yp_phot(aa),'.r')
hold off
xlabel('y [m]')
ylabel('Px [eV/c]')
subplot 224
set(gca,'FontSize',16)
hold on
plot(xp_phot,yp_phot,'.b')
plot(xp_phot(aa),yp_phot(aa),'.r')

hold off
xlabel('Px [eV/c]')
ylabel('Py [eV/c]')

suptitle({['For theta<' num2str(el_angel) ' [rad]'];['Total number of photons=' num2str(number_of_photons/turn_number,'%10.2e') ];['number scattered photons in theta<' num2str(length(full_spectrum(aa))*weigth/turn_number,'%10.2e') ]})

filename = [ out_folder 'photon_plots/photons_plot_' num2str(ifig) ];
fname = [ filename '.png'];
print('-dpng', fname);








%
%
% for nbin=[101];
%
%
%
%     % xedges = linspace(-2,2,nbin); yedges = linspace(0,50,nbin);
%     xedges = linspace(-50,50,nbin); yedges = linspace(5,50,nbin);
%     histmat = hist2(phot_angle, full_spectrum, xedges, yedges);
%
%     new_m=zeros(nbin,nbin);
%     for nni=1:1:nbin-1
%         new_m(nni,:)=histmat(nni,:)./((pi/2).*(abs((xedges(nni+1)).^2-(xedges(nni)).^2)));
%
%     end
%
%     figure(ifig)
%     ifig=ifig+1;
% %     set(pcolor(xedges,yedges,new_m'),'EdgeColor','none')
%  surf(xedges,yedges,new_m')
%     % colormap([1 1 1;0.857142865657806 0.857142865657806 1;0.714285731315613 0.714285731315613 1;0.571428596973419 0.571428596973419 1;0.428571432828903 0.428571432828903 1;0.28571429848671 0.28571429848671 1;0.142857149243355 0.142857149243355 1;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0]);
%     shading interp
%     colorbar
%     set(gca,'FontSize',16)
%     title({['with IBS']; })
%     set(gca,'FontSize',16)
%     xlabel('angle (mrad)');
%     ylabel('photons energy (KeV)');
%     filename = [pwd '/flux_plot/Poster'];
%     fname = [ filename '.png'];
%     print('-dpng', fname);
%
%
%
%     leng=zeros(nbin,1);
%
%     stda=zeros(nbin,1);
%
%     for lni=1:1:nbin
%         lni
%
%
%         if length(find(new_m(lni,:)))>0
%             leng(lni)=max(find(new_m(lni,:)))-min(find(new_m(lni,:)));
%         else
%             leng(lni)=0;
%         end
%
%
%         stda(lni)=std(new_m(lni,:));
%
%     end
%
%
%     figure(ifig)
%     ifig=ifig+1;
%     hold on
%     %  plot(new_m_1(37,:),'.r')
%     plot(xedges,leng*(yedges(3)-yedges(2)),'-or')
%
%     hold off
%     set(gca,'FontSize',16)
%
%     xlabel('angle (mrad)');
%     ylabel('photons energy (KeV)');
%     title(['nbin' num2str(nbin) ])
%     filename = [pwd '/flux_plot/broad_en_slice'];
%     fname = [ filename '.png'];
%     print('-dpng', fname);
%
%     figure(ifig)
%     ifig=ifig+1;
%     hold on
%     %  plot(new_m_1(37,:),'.r')
%     plot(xedges,smooth(leng)*(yedges(2)-yedges(1)),'--or','LineWidth',2)
%
%     hold off
%     grid on
%     set(gca,'FontSize',16)
%
%     xlabel('angle (mrad)');
%     ylabel('width of slice (KeV)');
%     %title(['smoothed nbin' num2str(nbin) ])
%     filename = [pwd '/flux_plot/broad_en_slice_smoothed'];
%     fname = [ filename '.eps'];
%     print('-depsc', fname);
%
%
%     figure(ifig)
%     ifig=ifig+1;
%     hold on
%     %  plot(new_m_1(37,:),'.r')
%     plot(xedges,stda*(yedges(3)-yedges(2)),'.r')
%
%     hold off
%     set(gca,'FontSize',16)
%     legend('CBS IBS','CBS',0)
%     xlabel('angle (mrad)');
%     ylabel('photons energy (KeV)');
%     title(['std nbin' num2str(nbin) ])
%     filename = [pwd '/flux_plot/std_of_slice'];
%     fname = [ filename '.png'];
%     print('-dpng', fname);
%
%
%
% end
%
% % % figure(ifig)
% % % ifig=ifig+1;
% % % hold on
% % % plot(smooth(hist(all_phot_1(en_1,10),nbin_plot)/length(all_phot_1)),'-r','LineWidth',2)
% % % plot(smooth(hist(all_phot_2(en_2,10),nbin_plot)/length(all_phot_2)),'-b','LineWidth',2)
% % % hold off
% % % grid on
% % % set(gca,'FontSize',16)
% % % ylabel('number of scattered photons')
% % % xlabel('Py')
% % % title(['for ph en=' num2str(ph_en) ])
% % % legend('CBS IBS','CBS',0)
% % % filename = [pwd '/flux_plot/Py_photons_' num2str(ph_en) ];
% % % fname = [ filename '.eps'];
% % % print('-depsc', fname);
% % %
% % % figure(ifig)
% % % ifig=ifig+1;
% % % hold on
% % % plot(smooth(hist(all_phot_1(en_1,10),nbin_plot)/length(all_phot_1)),'-r','LineWidth',2)
% % % plot(smooth(hist(all_phot_2(en_2,10),nbin_plot)/length(all_phot_2)),'-b','LineWidth',2)
% % % hold off
% % % grid on
% % % set(gca,'FontSize',16)
% % % ylabel('number of scattered photons')
% % % xlabel('Ps')
% % % title(['for ph en=' num2str(ph_en) ])
% % % legend('CBS IBS','CBS',0)
% % % filename = [pwd '/flux_plot/Ps_photons_' num2str(ph_en) ];
% % % fname = [ filename '.eps'];
% % % print('-depsc', fname);
% % end
% %
% %
% %
% %
% %
% for aaa=[0.01 0.1:0.1:1]
%     lml=zeros(nbin,1);
%
%     for lni=1:1:nbin
%         lni
%
%         %    val_1=find(new_m_1(lni,:));
%         %    val_2=find(new_m_2(lni,:));
%         %     leng_1(lni)=length(val_1)
%         %    leng_2(lni)=length(val_2)
%         if length(find(new_m(lni,:)))>0
%
%             m1=max(new_m(lni,:))*aaa
%             lml_1(lni)=length(find(new_m(lni,:)>m1))
%         else
%             leng_1(lni)=0;
%         end
%
%
%     end
%
%
%
%     figure(ifig)
%     ifig=ifig+1;
%     hold on
%     %  plot(lml_1,'.r')
%     %  plot(lml_2,'.b')
%     plot(xedges,smooth(lml)*(yedges(3)-yedges(2)),'--or','LineWidth',2)
%
%     hold off
%     grid on
%     set(gca,'FontSize',16)
%
%     xlabel('angle (mrad)');
%     ylabel('Spectrum FWHM (KeV)');
%     title(['estm=' num2str(aaa) ])
%     filename = [pwd '/flux_plot/fwms_' num2str(aaa*10)];
%     fname = [ filename '.eps'];
%     print('-depsc', fname);
% end
