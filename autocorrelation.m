clear;clc;%this version consider the dx in the experiment
tic
step=10;cr1=0.0;
addpath(genpath('Z:\Sheng\MechnoChem_data_analysis\strain\July_6_22\gap_1000'))
%read pigment data
linkg='Z:\Sheng\MechnoChem_data_analysis\oocyte_raw_data\granule_pos2\';
sindex=importdata('Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\230109_corr\matlab\data_index.txt');
%sindex=importdata('Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\230109_corr\matlab\data_index1.txt');
sda=sindex.data;
ste=sindex.textdata;

%read actin data
linke='Z:\Sheng\MechnoChem_data_analysis\oocyte_raw_data\';
No1='_C1';No2='_C1';

%for iii=1:size(sda,1)
for iii=1:43 %iii=43, 255-017
clearvars -except iii store1 phit wt tau1 tau2 ng1 step cr1 linkg sindex sda ste linke No2 No1


filename=ste(iii,2);filename=char(filename);fileno='_C1';
xL=sda(iii,2);yL=sda(iii,3);
ddt=sda(iii,5);dx=sda(iii,4);%new
t1=sda(iii,6);
meanIr=zeros(t1,1);meanIa=zeros(t1,1);
sx=max(size(0:step:floor(xL)-1));
sy=max(size(0:step:floor(yL)-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% calculate the cross correlation
rec1=0;
for i=2:sy-1
    for j=2:sx-1
        x0=(i-1)*step-floor(step/2);x1=(i-1)*step+floor(step/2);
        y0=(j-1)*step-floor(step/2);y1=(j-1)*step+floor(step/2);        
        for tt=1:t1
            img_k2 = imread([linke,filename,No1,'.tif'], tt) ; % read the image
            meanIr(tt)=mean(mean(img_k2(x0:x1,y0:y1)))-mean(mean(img_k2(:,:)));
            meanIa(tt)=mean(mean(img_k2(x0:x1,y0:y1)))-mean(mean(img_k2(:,:)));
        end
        
        %f1 = fit((1:t1)'*ddt,meanIr(1:t1),'Fourier1');
        f2 = fit((1:t1)'*ddt,meanIa(1:t1),'Fourier1');
        %cof1=coeffvalues(f1);
        cof2=coeffvalues(f2);
        
        clear normr norma
        normr(1:t1)=(meanIr(1:t1)-cof2(1));%*1/sqrt(cof2(2)^2+cof2(3)^2);
        norma(1:t1)=(meanIa(1:t1)-cof2(1));%*1/sqrt(cof2(2)^2+cof2(3)^2);
        
        [c,lag]=xcorr(norma,normr,'coef');
        %[c,lag]=xcorr(norma,normr);
        pm=round((size(lag,2)+1)/2);
        ct0(i,j,iii)=c(pm);ct(i,j,:)=c(pm:size(c,2))*(ct0(i,j,iii)>cr1);
        if ct0(i,j,iii)>cr1
                rec1=rec1+1;
        end

%%%%%% remove some bad cases
jj=1;
for ki=2:size(ct(i,j,:),3)-1
if ct(i,j,ki)>=ct(i,j,ki-1)&&ct(i,j,ki)>=ct(i,j,ki+1)
    jj=jj+1;
end
end

if jj<3
    rec1=rec1-1;ct(i,j,:)=0;
end
%%%%%% end of remove bad cases


    end
end

crrx1(iii,:)=lag(pm:size(c,2))*ddt;
crrx2(iii,:)=sum(sum(ct))/rec1;
%%%%% measure correlation variables
xf1=crrx1(iii,:);yf1=crrx2(iii,:);
f3 = fit(xf1(1:11)',yf1(1:11)','Fourier1');
cof3=coeffvalues(f3);
        
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

if ii<3
    f9 = fit(xf1(1:30)',yf1(1:30)','Exp1','Startpoint',[0,0]);cof9=coeffvalues(f9);
    store1(iii,3)=-1/cof9(2);
    store1(iii,5)=rec1;
    continue;
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

%%%%% calculating std
ki1=1;
for i=1:sy-1
    for j=1:sx-1
    if ct(i,j,:)~=0
        %for ki=1:size(ct(i,j,:),3)
        rec2(ki1,:)=ct(i,j,:);
        %end
        ki1=ki1+1;
    end
    end
end

ki2=size(ct(i,j,:),3);
for ki=1:ki2
std1(ki)=std(rec2(:,ki));    
end


%%%%% end of calculating std

%%%% output result
plotshaded(xf1,[yf1-std1;yf1+std1],'k');
hold on
plot(xf1,yf1,'k','linewidth',3);
%plot(xf1,exp(-xf1/tau1(iii,1))*max(yf1));
%plot(xf1,cof6(1)*exp(cof6(2)*xf1)+cof7(1)+cof7(2)*cos(cof7(4)*xf1)+cof7(3)*sin(cof7(4)*xf1));
plot(xf1,exp(-xf1/tau2(iii,1))*max(yf1),'linewidth',3);
ng1(iii,1)=rec1;
store1(iii,1:5)=[phit(iii,1),wt(iii,1),tau1(iii,1),tau2(iii,1),ng1(iii,1)];
axis square
set(gca,'linewidth',4)


saveas(gcf,['Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\Oct_11_2022\cross_corr\autocorrelation_Jan27_2023\plot5_shade\',filename,No1,'.fig'])
saveas(gcf,['Z:\Sheng\MechnoChem_data_analysis\ori_corr_analysis\Oct_11_2022\cross_corr\autocorrelation_Jan27_2023\plot5_shade\',filename,No1,'.png'])
close all


end

toc