clear all
close all

%define region switch
RegionName = 'Total';
%RegionName = 'GulfOfFinland';
%RegionName = 'BayOfBothnia';
%RegionName = 'BothniaSea';
%RegionName = 'box';
nbox = '04';% 02 03 04 05 07 08 13 15
folder = 'test';
std_folder = 'test';
switch RegionName
    case('Total')
        z1=importdata( strcat(std_folder,'/stdalone_total.txt') );
        z2=importdata( strcat(folder,'/coupled_total.txt') );
        myLab1=['Ice area [km^2]    x10^3'];
        myLab2=['Ice volume [km^3]'];
        fileOut2=['ice_vol_total.jpg'];
        fileOut1=['ice_area_total.jpg'];
        IceAreaNorm = 1000.;
    case('GulfOfFinland');
        z1=importdata( strcat(std_folder,'/stdalone_Gulf_of_Finland.txt') );
        z2=importdata( strcat(folder,'/coupled_Gulf_of_Finland.txt') );
        myLab1=['Ice area [km^2]    x10^3'];
        myLab2=['Ice volume [km^3]'];
        fileOut2=['ice_vol_GoF_new.jpg'];
        fileOut1=['ice_area_GoF_new.jpg'];
        IceAreaNorm = 1000.;
    case('BayOfBothnia');
        z1=importdata( strcat(std_folder,'/stdalone_Bothnian_Bay.txt') );
        z2=importdata( strcat(folder,'/coupled_Bothnian_Bay.txt') );
        myLab1=['Ice area [km^2]    x10^3'];
        myLab2=['Ice volume [km^3]'];
        fileOut2=['ice_vol_BoB.jpg'];
        fileOut1=['ice_area_BoB.jpg'];
        IceAreaNorm = 1000.;
    case('BothniaSea');
        z1=importdata( strcat(std_folder,'/stdalone_Bothnian_Sea.txt') );
        z2=importdata( strcat(folder,'/coupled_Bothnian_Sea.txt') );
        myLab1=['Ice area [km^2]    x10^3'];
        myLab2=['Ice volume [km^3]'];
        fileOut2=['ice_vol_BS.jpg'];
        fileOut1=['ice_area_BS.jpg'];
        IceAreaNorm = 1000.;
    case('box');
        z1=importdata([std_folder '/stdalone_box' nbox '.txt']);
        z2=importdata([folder '/coupled_box' nbox '.txt']);
        myLab1=['Ice area [km^2]         '];
        myLab2=['Ice volume [km^3]'];
        fileOut2=['ice_vol_box' nbox '.jpg'];
        fileOut1=['ice_area_box' nbox '.jpg'];
        IceAreaNorm = 1.;      
end;
N1=90;
N2=max(size(z1));
%time1=z1(N1:N2,1);
volumetot(:,1)=z1(N1:N2,5);
volumetot(:,2)=z2(N1:N2,5);
areatot(:,1)=z1(N1:N2,4);
areatot(:,2)=z2(N1:N2,4);
time=datenum(z1(N1:N2,1),z1(N1:N2,2),z1(N1:N2,3));
%time=1:max(size(areatot));
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.6, 0.5],'color','w');
plot(time,[areatot(:,1) areatot(:,2)]/IceAreaNorm,'Linewidth',3),legend('stand alone','coupled'),...
    ylabel(myLab1);
datetick('x','mmm','keepticks')
set(gca,'fontsize',18);
saveas(gcf,fileOut1);


figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.6, 0.5],'color','w');
plot(time,[volumetot(:,1) volumetot(:,2)],'Linewidth',3),legend('stand alone','coupled'),...
    ylabel(myLab2)
datetick('x','mmm','keepticks');
set(gca,'fontsize',18);
saveas(gcf,fileOut2);



