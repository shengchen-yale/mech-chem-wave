clear;clc;% calculate strain rates and viscous dissipation
tic
%step=50
%average over all time

addpath(genpath('Z:\Sheng\MechnoChem_data_analysis\strain\July_6_22\gap_1000'));
sindex=importdata('Z:\Sheng\MechnoChem_data_analysis\oocyte_raw_data\granule_pos2\data_index.txt');
sda=sindex.data;
ste=sindex.textdata;

ks=zeros(size(sda,1),10);
for iii=1:size(sda,1)
%for iii=9:13
ntmax=1;
vtrc=0.08;
vtrc2=0.0;

filename=ste(iii,2);filename=char(filename);fileno='_C1';
xL=sda(iii,2);yL=sda(iii,3);
ddt=sda(iii,5);dx=sda(iii,4);%new
maxdis=2.4/6*0.266*ddt/dx;
pos_lst=load([filename,'_pos.dat']);
t1=max(pos_lst(:,3))-1;

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


Atr=sortrows(C,4);vst=zeros(np,4);
vmax1=0.25;Nvr=25;probv=zeros(Nvr,4);

for i=1:np
    A2=(Atr(:,4)==i);
    idx = find(A2~=0, 1, 'first');
    nt=sum(A2);
    if isempty(idx)
        continue
    end
    
    vst(i,1)=mean(Atr(idx:(idx+nt-1),10));
    vst(i,2)=max(Atr(idx:(idx+nt-1),10));
    vst(i,3)=mean(Atr(idx:(idx+nt-1),11)); %after correction of shifting
    vst(i,4)=max(Atr(idx:(idx+nt-1),11));  %after correction of shifting
    
    ii=floor(vst(i,:)*Nvr/vmax1)+1;
    probv(ii(1),1)=probv(ii(1),1)+1;
    probv(ii(2),2)=probv(ii(2),2)+1;
    probv(ii(3),3)=probv(ii(3),3)+1;
    probv(ii(4),4)=probv(ii(4),4)+1;
    
end

%plot(1/Nvr*vmax1:1/Nvr*vmax1:vmax1,probv(:,1)/np)
%hold on
%plot(1/Nvr*vmax1:1/Nvr*vmax1:vmax1,probv(:,2)/np)
%plot(1/Nvr*vmax1:1/Nvr*vmax1:vmax1,probv(:,3)/np)
%plot(1/Nvr*vmax1:1/Nvr*vmax1:vmax1,probv(:,4)/np)

probv=probv/np;
probv(:,5)=(1/Nvr*vmax1:1/Nvr*vmax1:vmax1)';

%print out the velocity field
%fid = fopen([filename,fileno,'_vel.txt'],'wt');
%for ii = 1:size(probv,1)
%    fprintf(fid,'%g\t',probv(ii,:));
%    fprintf(fid,'\n');
%end
%fclose(fid);

%fitting to guassian distribution

f = fit(probv(:,5),probv(:,1),'gauss1');
cof=coeffvalues(f);
ks(iii,1)=cof(2);

f = fit(probv(:,5),probv(:,2),'gauss1');
cof=coeffvalues(f);
ks(iii,2)=cof(2);

f = fit(probv(:,5),probv(:,3),'gauss1');
cof=coeffvalues(f);
ks(iii,3)=cof(2);

f = fit(probv(:,5),probv(:,4),'gauss1');
cof=coeffvalues(f);
ks(iii,4)=cof(2);


%calculate the shear
step=10;
sx=max(size(0:step:floor(xL)));
sy=max(size(0:step:floor(yL)));
strain=zeros(sy,sx);
srate=zeros(t1,sy,sx);
shearr=zeros(t1,sy,sx);
srate2=zeros(t1,sy,sx);
rotr=zeros(t1,sy,sx);

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

%dilation
srate(k,2:sy-1,2:sx-1) = ...
    (vqy(3:sy,2:sx-1) - vqy(1:sy-2,2:sx-1))/2/dx/step ...
    +(vqx(2:sy-1,3:sx) - vqx(2:sy-1,1:sx-2))/2/dx/step; 

srate(isnan(srate))=0;

%shear 1
shearr(k,2:sy-1,2:sx-1) = ...
    (vqy(2:sy-1,3:sx) - vqy(2:sy-1,1:sx-2))/2/dx/step ...
    +(vqx(3:sy,2:sx-1) - vqx(1:sy-2,2:sx-1))/2/dx/step;

shearr(isnan(shearr))=0;

%shear 2
srate2(k,2:sy-1,2:sx-1) = ...
    -(vqy(3:sy,2:sx-1) - vqy(1:sy-2,2:sx-1))/2/dx/step ...
    +(vqx(2:sy-1,3:sx) - vqx(2:sy-1,1:sx-2))/2/dx/step; 

srate2(isnan(srate2))=0;

%rotation
rotr(k,2:sy-1,2:sx-1) = ...
    -(vqy(2:sy-1,3:sx) - vqy(2:sy-1,1:sx-2))/2/dx/step ...
    +(vqx(3:sy,2:sx-1) - vqx(1:sy-2,2:sx-1))/2/dx/step;

rotr(isnan(rotr))=0;

shearr(abs(shearr)>1)=0; %trucation
srate(abs(srate)>1)=0; %trucation
srate2(abs(srate2)>1)=0; %trucation
rotr(abs(rotr)>1)=0; %trucation
%strain=strain+srate;

Mediss(k,:,:)=shearr(k,:,:).*shearr(k,:,:);
end

%%%%%%%%%%%%%%%%% average over nearest frame
srate_n=srate;
shearr_n=shearr;
srate2_n=srate2;
rotr_n=rotr;

for ii=1:sx
    for jj=1:sy
        srate(:,jj,ii)=mean(srate_n(:,jj,ii));
        shearr(:,jj,ii)=mean(shearr_n(:,jj,ii));
        srate2(:,jj,ii)=mean(srate2_n(:,jj,ii));
        rotr(:,jj,ii)=mean(rotr_n(:,jj,ii));
    end
end

%%%%%%%%%%%%%%%%% end of average

%Mediss=shearr.^2;
Asrate2=abs(shearr);
Asrate1=abs(srate);
Asrate3=abs(srate2);
Asrate4=abs(rotr);
%calculate mean strain rate
sr1=zeros(t1,10);

for ii=1:t1
    sr1(ii,1)=mean(mean(Asrate2(ii,:,:)));
    pp=isinf(1./Asrate2(ii,:,:));
    [ix,iy,iz]=size(Asrate2);
    sump=iy*iz-sum(sum(pp));
    if sump == 0
        sr1(ii,2)=0;
    else
        sr1(ii,2)=sum(sum(Asrate2(ii,:,:)))/sump;
    end
    
    sr1(ii,3)=mean(mean(Mediss(ii,:,:)));
    pp=isinf(1./Mediss(ii,:,:));
    [ix,iy,iz]=size(Mediss);
    sump=iy*iz-sum(sum(pp));
    if sump == 0
        sr1(ii,4)=0;
    else
        sr1(ii,4)=sum(sum(Mediss(ii,:,:)))/sump;
    end
    
    sr1(ii,5)=mean(mean(Asrate1(ii,:,:)));
    pp=isinf(1./Asrate1(ii,:,:));
    [ix,iy,iz]=size(Asrate1);
    sump=iy*iz-sum(sum(pp));
    if sump == 0
        sr1(ii,6)=0;
    else
        sr1(ii,6)=sum(sum(Asrate1(ii,:,:)))/sump;
    end
    
    %shear2
    sr1(ii,7)=mean(mean(Asrate3(ii,:,:)));
    pp=isinf(1./Asrate3(ii,:,:));
    [ix,iy,iz]=size(Asrate3);
    sump=iy*iz-sum(sum(pp));
    if sump == 0
        sr1(ii,8)=0;
    else
        sr1(ii,8)=sum(sum(Asrate3(ii,:,:)))/sump;
    end
    
    %rotate
    sr1(ii,9)=mean(mean(Asrate4(ii,:,:)));
    pp=isinf(1./Asrate4(ii,:,:));
    [ix,iy,iz]=size(Asrate4);
    sump=iy*iz-sum(sum(pp));
    if sump == 0
        sr1(ii,10)=0;
    else
        sr1(ii,10)=sum(sum(Asrate4(ii,:,:)))/sump;
    end
end
    
store(1,1)=mean(sr1(:,5));
store(1,2)=mean(sr1(:,6));
store(1,3)=mean(sr1(:,1)); %avarage over all
store(1,4)=mean(sr1(:,2)); %remove the zeros 
store(1,5)=mean(sr1(:,3));
store(1,6)=mean(sr1(:,4));
%shear2
store(1,7)=mean(sr1(:,7));
store(1,8)=mean(sr1(:,8));
%rotate
store(1,9)=mean(sr1(:,9));
store(1,10)=mean(sr1(:,10));

ks(iii,5)=store(1,1);
ks(iii,6)=store(1,2);
ks(iii,7)=store(1,3);
ks(iii,8)=store(1,4);
ks(iii,9)=store(1,5);
ks(iii,10)=store(1,6);

%shear2,rotate
ks(iii,11)=store(1,7);
ks(iii,12)=store(1,8);
ks(iii,13)=store(1,9);
ks(iii,14)=store(1,10);
end
toc