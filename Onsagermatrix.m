%this version only consider the actin & velocity, ev \codt gradA
clear;clc;
tic
kgapc='k06_';kgap=0.6;cas='1';
%tv=3:203;
dt=5;

%import data
R=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\R_',kgapc,cas,'.txt']);
F=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\F_',kgapc,cas,'.txt']);
Fx=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\vgradF\Fx_',kgapc,cas,'.txt']);
Fy=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\vgradF\Fy_',kgapc,cas,'.txt']);
ux=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\vgradF\ux_',kgapc,cas,'.txt']);
uy=importdata(['Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\vgradF\uy_',kgapc,cas,'.txt']);


%parameter
kect=0.04827;A1=0.05;k2=0.5361;G1=0.01;Km=10;ka=0.08/10*0.3+kgap*0.15/10;%ka=0.0024;
SF=3e-5;kd=0.02484;S=-0.0045;
kb=1.380649*1e-23;Tem=293.15;
%%%%%%%%%%%


it1=0;
for ix1=1:5:10201
for tv=4:3:193
    
y1=R(ix1,tv-1);y2=F(ix1,tv-1);y31=1e6*sqrt(ux(ix1,tv-1)^2+uy(ix1,tv-1)^2);if y1<0||y2<0||y31<1e-3,continue;end
jRp(1)=kect*y1/(A1+y1);jRm(1)=k2*y2*y1/(G1+y1)+kgap*y1/(Km+y1)-S;
jFp(1)=ka*y1+SF;jFm(1)=kd*y2;
jMp(1)=Fx(ix1,tv-1)*ux(ix1,tv-1)/y31;jMm(1)=Fy(ix1,tv-1)*uy(ix1,tv-1)/y31;

y1=R(ix1,tv);y2=F(ix1,tv);y32=1e6*sqrt(ux(ix1,tv)^2+uy(ix1,tv)^2);if y1<0||y2<0||y32<1e-3,continue;end
jRp(2)=kect*y1/(A1+y1);jRm(2)=k2*y2*y1/(G1+y1)+kgap*y1/(Km+y1)-S;
jFp(2)=ka*y1+SF;jFm(2)=kd*y2;
jMp(2)=Fx(ix1,tv)*ux(ix1,tv)/y32;jMm(2)=Fy(ix1,tv)*uy(ix1,tv)/y32;

y1=R(ix1,tv+1);y2=F(ix1,tv+1);y33=1e6*sqrt(ux(ix1,tv+1)^2+uy(ix1,tv+1)^2);if y1<0||y2<0||y33<1e-3,continue;end
jRp(3)=kect*y1/(A1+y1);jRm(3)=k2*y2*y1/(G1+y1)+kgap*y1/(Km+y1)-S;
jFp(3)=ka*y1+SF;jFm(3)=kd*y2;
jMp(3)=Fx(ix1,tv+1)*ux(ix1,tv+1)/y33;jMm(3)=Fy(ix1,tv+1)*uy(ix1,tv+1)/y33;

Aff=[log(jRp(1)/jRm(1)),log(jFp(1)/jFm(1)),jMp(1)+jMm(1);
     log(jRp(2)/jRm(2)),log(jFp(2)/jFm(2)),jMp(2)+jMm(2)
     log(jRp(3)/jRm(3)),log(jFp(3)/jFm(3)),jMp(3)+jMm(3)];

p1=max(max(isnan(Aff)))*max(max(isinf(Aff)))*(rank(Aff)<3);
if p1~=1
   it1=it1+1;
   L1=Aff\[jRp(1)-jRm(1);jRp(2)-jRm(2);jRp(3)-jRm(3)]; %L11 L12 L13
   L2=Aff\[jFp(1)-jFm(1);jFp(2)-jFm(2);jFp(3)-jFm(3)]; %L21 L22 L23
   L3=Aff\[y31*F(ix1,tv-1);y32*F(ix1,tv);y33*F(ix1,tv+1)]; %L31 L32 L33
   Lij(it1,1:9)=[L1',L2',L3']; %
   %scalue
   Lnij(it1,1)=Lij(it1,2)/(Lij(it1,1)+Lij(it1,5)); %L12/(L11+L22)
   Lnij(it1,2)=Lij(it1,4)/(Lij(it1,1)+Lij(it1,5)); %L21/(L11+L22)
   Lnij(it1,3)=Lij(it1,6)/(Lij(it1,5)+Lij(it1,9)); %L23/(L22+L33)
   Lnij(it1,4)=Lij(it1,8)/(Lij(it1,5)+Lij(it1,9)); %L32/(L22+L33)
   Lnij(it1,5)=Lij(it1,2)/(Lij(it1,1)+Lij(it1,5)+Lij(it1,9)); %L12/(L11+L22+L33)
   Lnij(it1,6)=Lij(it1,4)/(Lij(it1,1)+Lij(it1,5)+Lij(it1,9)); %L21/(L11+L22+L33)
   Lnij(it1,7)=Lij(it1,6)/(Lij(it1,1)+Lij(it1,5)+Lij(it1,9)); %L23/(L11+L22+L33)
   Lnij(it1,8)=Lij(it1,8)/(Lij(it1,1)+Lij(it1,5)+Lij(it1,9)); %L32/(L11+L22+L33)
   Lnij(it1,9)=abs(Lij(it1,2)-Lij(it1,4))/(Lij(it1,1)+Lij(it1,5)); %abs(L12-L21)/(L11+L22)
   Lnij(it1,10)=abs(Lij(it1,6)-Lij(it1,8))/(Lij(it1,5)+Lij(it1,9)); %abs(L23-L32)/(L22+L33)
   Lnij(it1,11)=abs(Lij(it1,3)-Lij(it1,7))/(Lij(it1,1)+Lij(it1,9)); %abs(L13-L31)/(L11+L33)
end

%add filter
if L1(1)>1||L1(1)<0||L2(2)>0.1||L2(2)<0||L1(2)>0.01||L1(2)<-0.01||L2(1)>0.01||L2(1)<-0.01||L3(3)>1||L3(3)<0||L3(2)>1e-2||L3(2)<-1e-2
it1=it1-1;
end


end
end
%%%

toc


store1=[mean(Lij(:,1)),mean(Lij(:,2)),mean(Lij(:,3)),mean(Lij(:,4)),mean(Lij(:,5)),mean(Lij(:,6)),mean(Lij(:,7)),mean(Lij(:,8)),mean(Lij(:,9));
        mean(Lnij(:,1)),mean(Lnij(:,2)),mean(Lnij(:,3)),mean(Lnij(:,4)),mean(Lnij(:,5)),mean(Lnij(:,6)),mean(Lnij(:,7)),mean(Lnij(:,8)),0];

store1(3,1:3)=[mean(Lnij(:,9)),mean(Lnij(:,10)),mean(Lnij(:,11))];

load('Z:\Sheng\simulation\comsol\Sep21_noise_test\result_file\RF_noise_R8_corrected_reactional_flow\data_reciprocal_test1\vgradF\mfile\data1.mat')
%pointcloud plot
figure(1)
x1=zeros(1001,1)+kgap;
y1=Lnij(4000:5000,1);y2=Lnij(4000:5000,9);
hold on
scatter(x1,y2,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha','0.04')
plot(s(:,1),s(:,21),'Color',[0.64 0.08 0.18],'LineStyle','-','Marker','o','Markersize',12,'LineWidth',3)

axis square
xlim([-0.02 0.62])
ylim([-0.05 0.45])
set(gca,'fontsize',20)
set(gca,'LineWidth',3)
plot([-0.05 0.65],[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',3)
box on;

figure(2)
x1=zeros(1001,1)+kgap;
y1=Lnij(4000:5000,3);y2=Lnij(4000:5000,10);
hold on
scatter(x1,y2,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha','0.04')
plot(s(:,1),s(:,22),'Color',[0 0 1],'LineStyle','-','Marker','o','Markersize',12,'LineWidth',3)

axis square
xlim([-0.02 0.62])
ylim([-0.5 5.5])
set(gca,'fontsize',20)
set(gca,'LineWidth',3)
plot([-0.05 0.65],[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',3)
box on;

figure(3)
x1=zeros(1001,1)+kgap;
y1=Lnij(4000:5000,3);y2=Lnij(4000:5000,11);
%scatter(x1,y1,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha','0.01')
hold on
scatter(x1,y2,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha','0.04')
plot(s(:,1),s(:,23),'Color',[1 0.41 0.16],'LineStyle','-','Marker','o','Markersize',12,'LineWidth',3)

axis square
xlim([-0.02 0.62])
ylim([-1 11])
set(gca,'fontsize',20)
set(gca,'LineWidth',3)
plot([-0.05 0.65],[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',3)
box on;