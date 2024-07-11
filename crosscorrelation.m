clear;clc;%this version consider the dx in the experiment

step=10;cr1=0.1;
addpath(genpath('Z:\Sheng\MechnoChem_data_analysis\strain\July_6_22\gap_1000'))
%read pigment data
linkg='Z:\Sheng\MechnoChem_data_analysis\oocyte_raw_data\granule_pos2\';
sindex=importdata('Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\230109_corr\matlab\data_index.txt');
%sindex=importdata('Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\230109_corr\matlab\data_index1.txt');
sda=sindex.data;
ste=sindex.textdata;

%read actin data
linke='Z:\Sheng\MechnoChem_data_analysis\oocyte_raw_data\';
No2='_C2';

%for iii=1:size(sda,1)
for iii=1:43 %iii=43, 255-017
clearvars -except iii store1 phit wt tau1 tau2 ng1 step cr1 linkg sindex sda ste linke No2

ntmax=1;
%vtrc=0.08;
vtrc=1;
vtrc2=0.0;

filename=ste(iii,2);filename=char(filename);fileno='_C1';
xL=sda(iii,2);yL=sda(iii,3);
ddt=sda(iii,5);dx=sda(iii,4);%new
maxdis=2.4/6*0.266*ddt/dx;
pos_lst=load([linkg,filename,'_pos.dat']);
t1=max(pos_lst(:,3))-1;
meanIr=zeros(t1,1);meanIa=zeros(t1,1);

tr = track(pos_lst, maxdis); %result = track( positionlist, maxdisp, param )
A=sortrows(tr,4);


%filter remove isolated particles
np=max(A(:,4));
i1=1;
for i=1:np
    A1=A;
    A2=(A1(:,4)==i);
    idx = find(A2~=0, 1, 'first');
    nt=sum(A2);
    if nt>ntmax
    A3(i1:(i1+nt-1),1:4)=A1(idx:(idx+nt-1),1:4);
    i1=i1+nt;
    end
end
A=A3;
numelement=max(size(A));
%end of filter


for i=1:(numelement-1)

if (A(i+1,3)-A(i,3)<0)
A(i,6)=A(i-1,6);A(i,7)=A(i-1,7);
elseif (A(i+1,3)-A(i,3)==0)   %isolate single granule, no need when adding filter
A(i,6)=0;A(i,7)=0;    
else
A(i,6)=dx*(A(i+1,1)-A(i,1))/(A(i+1,3)-A(i,3))/ddt; %%%%%%modified consider dx here
A(i,7)=dx*(A(i+1,2)-A(i,2))/(A(i+1,3)-A(i,3))/ddt; %%%%%%modified consider dx here
end

vmag=sqrt(A(i,6)^2+A(i,7)^2);
if (vmag>vtrc||vmag<vtrc2)
    A(i,6)=0;A(i,7)=0;
end
end
C=sortrows(A,3);

%%record the row index for each time step
t0=1;%t1=75;
dt=1;
st=zeros(t1,3);
for i=1:numelement
j=C(i,3)+1;
st(j:t1,3)=st(j:t1,3)+1;
end

st(1,1)=1;st(1,2)=1;

for k=(t0+1):t1
    st(k,1)=k;
    st(k,2)=st(k-1,3)+1;
end

for k=1:t1
C(st(k,2):st(k,3),8)=C(st(k,2):st(k,3),6)-mean(C(st(k,2):st(k,3),6));
C(st(k,2):st(k,3),9)=C(st(k,2):st(k,3),7)-mean(C(st(k,2):st(k,3),7));
C(st(k,2):st(k,3),10)=sqrt(C(st(k,2):st(k,3),6).^2+C(st(k,2):st(k,3),7).^2); %before truncation of shifting
C(st(k,2):st(k,3),11)=sqrt(C(st(k,2):st(k,3),8).^2+C(st(k,2):st(k,3),9).^2); %after truncation of shifting
end


%calculate the shear

sx=max(size(0:step:floor(xL)-1));
sy=max(size(0:step:floor(yL)-1));
strain=zeros(sy,sx);
srate=zeros(t1,sy,sx);
shearr=zeros(t1,sy,sx);
Mediss=zeros(t1,sy,sx);

for k=1:t1

scale_factor = 1;
xx1=C(st(k,2):st(k,3),1);
yy1=C(st(k,2):st(k,3),2);
uu1=C(st(k,2):st(k,3),6)*scale_factor;
vv1=C(st(k,2):st(k,3),7)*scale_factor;

uu1(isnan(uu1))=0;
vv1(isnan(vv1))=0;
vm=sqrt(uu1.^2+vv1.^2);

[xq,yq]=meshgrid(0:step:floor(xL),0:step:floor(yL)); %resolution of interpolation
vqx=griddata(xx1,yy1,uu1,xq,yq,'linear');
vqy=griddata(xx1,yy1,vv1,xq,yq,'linear');
vq = griddata(xx1,yy1,vm,xq,yq,'linear');

srate(k,2:sy-1,2:sx-1) = ...
    (vqy(3:sy,2:sx-1) - vqy(1:sy-2,2:sx-1))/2/dx/step ...
    +(vqx(2:sy-1,3:sx) - vqx(2:sy-1,1:sx-2))/2/dx/step; 

srate(isnan(srate))=0;

shearr(k,2:sy-1,2:sx-1) = ...
    (vqy(2:sy-1,3:sx) - vqy(2:sy-1,1:sx-2))/2/dx/step ...
    +(vqx(3:sy,2:sx-1) - vqx(1:sy-2,2:sx-1))/2/dx/step; 

shearr(isnan(shearr))=0;
shearr(abs(shearr)>1)=0; %trucation
srate(abs(srate)>1)=0; %trucation
%strain=strain+srate;

Mediss(k,:,:)=shearr(k,:,:).*shearr(k,:,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% calculate the cross correlation
rec1=0;
for i=2:sy-1
%    for j=2:sy-1
%for i=18:23
    for j=2:sx-1
        x0=(i-1)*step-floor(step/2);x1=(i-1)*step+floor(step/2);
        y0=(j-1)*step-floor(step/2);y1=(j-1)*step+floor(step/2);        
        for tt=1:t1
            img_k2 = imread([linke,filename,No2,'.tif'], tt) ; % read the image
            meanIr(tt)=srate(tt,i,j);
            meanIa(tt)=mean(mean(img_k2(x0:x1,y0:y1)));
        end
        
        %f1 = fit((1:t1)'*ddt,meanIr(1:t1),'Fourier1');
        f2 = fit((1:t1)'*ddt,meanIa(1:t1),'Fourier1');
        %cof1=coeffvalues(f1);
        cof2=coeffvalues(f2);
        
        clear normr norma
        %normr(1:t1)=(meanIr(1:t1)-cof1(1))*1/sqrt(cof1(2)^2+cof1(3)^2);
        normr(1:t1)=meanIr(1:t1);
        norma(1:t1)=(meanIa(1:t1)-cof2(1))*1/sqrt(cof2(2)^2+cof2(3)^2);
        
        [c,lag]=xcorr(norma,normr);
        pm=round((size(lag,2)+1)/2);
        ct0(i,j,iii)=c(pm);ct(i,j,:)=c(pm:size(c,2))*(ct0(i,j,iii)>cr1);
        if ct0(i,j,iii)>0.1
            rec1=rec1+1;
        end

    end
end

crrx1(iii,:)=lag(pm:size(c,2))*ddt;
crrx2(iii,:)=sum(sum(ct))/rec1;
%%%%% measure correlation variables
xf1=crrx1(iii,:);yf1=crrx2(iii,:);
f3 = fit(xf1(1:11)',yf1(1:11)','Fourier1');
cof3=coeffvalues(f3);
%phit(iii,1)=abs(acos(cof3(2)/sqrt(cof3(2)^2+cof3(3)^2))/cof3(4));
%f4 = fit(xf1(1:21)',yf1(1:21)','Fourier1');
%cof4=coeffvalues(f4);
%wt(iii,1)=cof3(4);
        
%measure tau1
f5=fit(xf1(1:15)',yf1(1:15)',['exp(-a*x)*(',num2str(cof3(1)),'+',num2str(cof3(2)),'*cos(',num2str(cof3(4)),'*x)+',num2str(cof3(3)),'*sin(',num2str(cof3(4)),'*x))'],'Start', 0);
cof5=coeffvalues(f5);
tau1(iii,1)=1/cof5;
%measure tau2
ii=1;xf2(ii,1)=xf1(1,1);yf2(ii,1)=yf1(1,1);
for ki=2:size(xf1,2)-1
if yf1(ki)>=yf1(ki-1)&&yf1(ki)>=yf1(ki+1)
    ii=ii+1;
    xf2(ii,1)=xf1(1,ki);yf2(ii,1)=yf1(1,ki);
end
end

f8 = fit(xf2,yf2,'Exp1','Startpoint',[0,0]);cof8=coeffvalues(f8);
tau2(iii,1)=-1/cof8(2);
%measure w and phi
f6 = fit(xf1',yf1','Exp1','Startpoint',[0,0]);cof6=coeffvalues(f6);
%plot(xf1,exp(xf1*cof6(2))*max(yf1));

xp2=crrx1(iii,:);yp2=yf1-(cof6(1)*exp(cof6(2)*xf1));
f7 = fit(xp2',yp2','Fourier1');cof7=coeffvalues(f7);

phit(iii,1)=asin(cof7(3)/sqrt(cof7(2)^2+cof7(3)^2))/cof7(4);
wt(iii,1)=cof7(4);
%%%%%%%%% end of measuring corr variables

%%%%%%% end of calculate cross correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(xf1,yf1);
hold on
plot(xf1,exp(-xf1/tau1(iii,1))*max(yf1));
plot(xf1,cof6(1)*exp(cof6(2)*xf1)+cof7(1)+cof7(2)*cos(cof7(4)*xf1)+cof7(3)*sin(cof7(4)*xf1));
plot(xf1,exp(-xf1/tau2(iii,1))*max(yf1));
ng1(iii,1)=rec1;
store1(iii,1:5)=[phit(iii,1),wt(iii,1),tau1(iii,1),tau2(iii,1),ng1(iii,1)];

%saveas(gcf,['Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\230109_corr\matlab\plot\',filename,'.png'])
close all
end

