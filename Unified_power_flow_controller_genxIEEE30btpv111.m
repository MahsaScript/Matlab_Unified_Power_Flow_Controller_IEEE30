
clc
clear
t0=clock;
load E30BTKMX.m
AT=E30BTKMX;
nb=AT(1,1);
itr=AT(1,2);
rtr=AT(1,3);
tp=rtr;
nrtr=AT(1,4);
trl=AT(1,5);
sz=AT(1,6);
sv=AT(1,7);
pq=AT(2,1);
g=AT(2,2);
rm=nb-g-sv;
tol=AT(2,3);
l=rtr+nrtr+trl;
ld=nb-g;
AL=AT(3:l+sz+2,:);
AD=AT(l+sz+3:l+sz+2+pq,:);
Ag=AT(l+sz+3+pq:l+sz+2+pq+g,:);
conv=3.14159/180.0;
se=AL(:,1);
re=AL(:,2);

for pl=1:l
    for ql=1:l
        if (AL(ql,1)==se(pl))&(AL(ql,2)==re(pl))
            y(ql)=(1/(AL(ql,3)+AL(ql,4)*i));   
        end
    end 
end
y;

for pl=1:l
    for ql=1:l
        if (se(ql)==se(pl,:))&(re(ql)==re(pl,:))
            Y(se(pl,:),re(pl,:))=-y(ql);
            Y(re(pl,:),se(pl,:))=-y(ql);
        end
    end
end
Y;

for p=1:nb
    for q=1:nb
        Y(p,p)=0.0;
    end
end

for p=1:nb
    for q=1:nb
        if (p~=q)
            Y(p,p)=(Y(p,p)-Y(p,q));
        end
    end
end
YTTNT=Y;

for pt=1:rtr
    for qt=1:rtr
        if (se(qt)==se(pt,:))&(re(qt)==re(pt,:))
            wt(qt)=AL(qt,7);
        end
    end
end

for pt=1:rtr
    for qt=1:rtr
        if (se(qt)==se(pt,:))&(re(qt)==re(pt,:))
            Y(se(pt,:),re(pt,:))=-y(qt)/(wt(qt)*wt(qt));
        end
    end
end

for p=1:nb
    for q=1:nb
        Y(p,p)=0.0;
    end
end

for p=1:nb
    for q=1:nb
        if (p~=q)
            Y(p,p)=(Y(p,p)-Y(p,q));
        end
    end
end

for pt=1:rtr
    for qt=1:rtr
        if (se(qt)==se(pt,:))&(re(qt)==re(pt,:))
            Y(re(pt,:),se(pt,:))=-y(qt)/wt(qt);
            Y(se(pt,:),re(pt,:))=-y(qt)/wt(qt);
        end
    end
end
Y;


if sz~=0
for pz=l+1:l+sz
   for qz=l+1:l+sz
          if (AL(qz,1)==se(pz))
            yz(qz-l+1)=1/(AL(qz,3)+AL(qz,4)*i);  
         end    
      end 
   end
   
yz;

YZ=zeros(nb,nb);

for pz=l+1:l+sz
    if(YZ(AL(pz,1),AL(pz,1))==(0+0*j))
        for qz=l+1:l+sz
            if (AL(pz,1)==AL(qz,1)) 
                YZ(AL(pz,1),AL(pz,1))=YZ(AL(pz,1),AL(pz,1))+yz(qz-l+1);
            end    
        end 
    end
end
YZ;                           
else
end

YB=zeros(nb,nb);
for p=1:nb
    for ql=1:l
        if (AL(ql,1)==p)|(AL(ql,2)==p)
            YB(p,p)=YB(p,p)+(AL(ql,6)*i);
        end    
    end
end
YB;
  
for pl=1:l
    for ql=1:l
        if (se(ql)==se(pl,:))&(re(ql)==re(pl,:))
            Ygb2(se(pl,:),re(pl,:))=(AL(ql,5)+AL(ql,6)*i);
            Ygb2(re(pl,:),se(pl,:))=(AL(ql,5)+AL(ql,6)*i);
        end
    end
end
Ygb2;
g2=real(Ygb2);
b2=imag(Ygb2);


if sz~=0
    YTT=[Y+YB+YZ];
else
    YTT=[Y+YB];
end


if sz~=0
    YTTNTT=[YTTNT+YB+YZ];
else
    YTTNTT=[YTTNT+YB];
end

YTTNTT;

BL=YTT(g+1:nb,:);
YLG=BL(:,1:g);
YLL=BL(:,g+1:nb);
YLLI=inv(YLL);
FLG=-((YLLI)*(YLG));

FLGM=abs(FLG);
FLGA=angle(FLG);

for p=1:nb
    for q=1:nb
        YTTM(p,q)=abs(YTT(p,q));
        YTTA(p,q)=angle(YTT(p,q));
    end
end
YTTM;
YTTA;
YTTMA=abs(YTTNTT);
YTTAA=angle(YTTNTT);

%formation of G and B matrix

for p=1:nb
    for q=1:nb
        G(p,q)=real(YTT(p,q));
        B(p,q)=imag(YTT(p,q));
    end
end
rd=AD(:,1);
rg=Ag(:,1);

	
Qmax=zeros(g,1);
for mg=1:g
    if (Ag(mg,1)==rg(mg))
            Qmax(mg)=Qmax(mg)+Ag(mg,3);
    end
end
% Qmax=Qmax/100;

Qmin=zeros(g,1);
for mg=1:g
    if (Ag(mg,1)==rg(mg))
            Qmin(mg)=Qmin(mg)+Ag(mg,4);
    end
end
% Qmin=Qmin/100;
 
Ps=zeros(nb,1);
for mg=1:g
    if (Ag(mg,1)==rg(mg))
        Ps(mg)=Ps(mg)+Ag(mg,2);
    end
end

for md=1:pq
    if (AD(md,1)==rd(md))
        Ps(rd(md))=Ps(rd(md))-AD(md,2);
    end
end
% Ps=Ps/100;

Qs=zeros(nb,1);
for md=1:pq
    if (AD(md,1)==rd(md))
        Qs(rd(md))=Qs(rd(md))-AD(md,3);
    end
end
% Qs=Qs/100;


VM=ones(nb,1);
for mg=1:g
    if (Ag(mg,1)==rg(mg))
        VM(mg)=Ag(mg,5);
    end
end
VM;
DL=zeros(nb,1);
VMvr=1.0;
DLvr=0*conv;

ku=9;mu=31;


PMref=0.4;
QMref=0.25;

Xcr=0.1;
VM0=1.0;

DLcr=-atan(PMref/abs(QMref));
VMcr=(Xcr/VM0)*sqrt(PMref*PMref+QMref*QMref);

Pmks=-PMref;
Pbbs=0;
Qmks=-QMref;

zcr=(.1*i);
zvr=(.1*i);

yvr=1/zvr;
ycr=1/zcr;

gcr=real(ycr);
bcr=imag(ycr);
gvr=real(yvr);
bvr=imag(yvr);

Gkm=-gcr;
Bkm=-bcr;
Gmm=gcr;
Bmm=(bcr);
Gkk=gcr+gvr;
Bkk=bcr+bvr;

Gmk=Gkm;
Bmk=Bkm;

Gvr=-gvr;
Bvr=-bvr;

itn=0.0;

while(itn<=itr) 
  
P=zeros(nb,1);
Q=zeros(nb,1);
for p=1:nb
   for q=1:nb
      P(p)=P(p)+VM(p)*VM(q)*YTTM(p,q)*cos(DL(p)-DL(q)-YTTA(p,q));
      Q(p)=Q(p)+VM(p)*VM(q)*YTTM(p,q)*sin(DL(p)-DL(q)-YTTA(p,q)); 
   end
end
	P;
	Q;
 
Pk=VM(ku)*VM(ku)*gvr+VM(ku)*VMcr*(Gkm*cos(DL(ku)-DLcr)+Bkm*sin(DL(ku)-DLcr))+VM(ku)*VMvr*(Gvr*cos(DL(ku)-DLvr)+Bvr*sin(DL(ku)-DLvr));
Pm=VM(mu)*VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
Qk=-VM(ku)*VM(ku)*bvr+VM(ku)*VMcr*(Gkm*sin(DL(ku)-DLcr)-Bkm*cos(DL(ku)-DLcr))+VM(ku)*VMvr*(Gvr*sin(DL(ku)-DLvr)-Bvr*cos(DL(ku)-DLvr));
Qm=VM(mu)*VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));

Pmk=VM(mu)*VM(mu)*Gmm+VM(ku)*VM(mu)*(Gkm*cos(DL(mu)-DL(ku))+Bkm*sin(DL(mu)-DL(ku)))+VM(mu)*VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
Qmk=-VM(mu)*VM(mu)*Bmm+VM(mu)*VM(ku)*(Gkm*sin(DL(mu)-DL(ku))-Bkm*cos(DL(mu)-DL(ku)))+VM(mu)*VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));

Pcr=VMcr*VMcr*Gmm+VM(ku)*VMcr*(Gkm*cos(DLcr-DL(ku))+Bkm*sin(DLcr-DL(ku)))+VM(mu)*VMcr*(Gmm*cos(DLcr-DL(mu))+Bmm*sin(DLcr-DL(mu)));
Pvr=-VMvr*VMvr*Gvr+VMvr*VM(ku)*(Gvr*cos(DLvr-DL(ku))+Bvr*sin(DLvr-DL(ku)));
Pbb=Pcr+Pvr;

Qcr=-VMcr*VMcr*Bmm+VM(ku)*VMcr*(Gkm*sin(DLcr-DL(ku))-Bkm*cos(DLcr-DL(ku)))+VM(mu)*VMcr*(Gmm*sin(DLcr-DL(mu))-Bmm*cos(DLcr-DL(mu)));
Qvr=VMvr*VMvr*Bvr+VMvr*VM(ku)*(Gvr*sin(DLvr-DL(ku))-Bvr*cos(DLvr-DL(ku)));


JPk=VM(ku)*(Gvr*cos(DL(ku)-DLvr)+Bvr*sin(DL(ku)-DLvr));
JQk=VM(ku)*VMvr*(Gvr*sin(DL(ku)-DLvr)-Bvr*cos(DL(ku)-DLvr));
JPbb=-2*VMvr*gvr+VM(ku)*(Gvr*cos(DLvr-DL(ku))+Bvr*sin(DLvr-DL(ku)));
Pbb=Pcr+Pvr;

P(ku)=P(ku)+Pk;
P(mu)=P(mu)+Pm;
Q(ku)=Q(ku)+Qk;
Q(mu)=Q(mu)+Qm;
P;
Q;
for p=1:nb
   DP(p)=	Ps(p)-P(p);
end
DPA=DP(:,2:nb);
for p=1:nb
   DQ(p)=	Qs(p)-Q(p);
end
DQA=DQ(:,g+1:nb);
DPmk=Pmks-Pmk;
DQmk=Qmks-Qmk;
DPbb=Pbbs-Pbb;
DPMP=[DPA';DPmk;DPbb];
DPMQ=[DQA';DQmk];

DPQABS=[abs(DPMP);abs(DPMQ)];

DPQ=[DPMP;DPMQ];

if (DPQABS<=tol)
   disp('THE SOLUTION IS ------------------------')
    
    minsing=min(svd(J))
    Pmk
    Qmk
    VMcr
    DLcrD=DLcr/conv
    VMvr
    DLvrD=DLvr/conv
    Pbb
    Pcr
    Pvr
    Qcr
    Qvr
    VM
    DLD=DLD'
    P
    Q
    Ploss=zeros(1,1);
    for p=1:nb
        for jp=1
            Ploss=Ploss+P(p,jp);
        end
    end
    Ploss
    Qloss=zeros(1,1);
    for p=1:nb
        for jp=1
            Qloss=Qloss+Q(p,jp);
        end
    end
    Qloss
    
    for p=1:nb
   for q=1:nb
      if (q~=p)&(YTTM(p,q)>0.001)
         PF(p,q)=-VM(p)*VM(p)*YTTM(p,q)*cos(YTTA(p,q))+VM(p)*VM(q)*YTTM(p,q)*cos(DL(p)-DL(q)-YTTA(p,q));
         QF(p,q)=VM(p)*VM(p)*YTTM(p,q)*sin(YTTA(p,q))+VM(p)*VM(q)*YTTM(p,q)*sin(DL(p)-DL(q)-YTTA(p,q))-VM(p)*VM(p)*b2(p,q);     
      p;,q;,PF(p,q);,QF(p,q);
     end
   end
end



for p=1:nb
       VC(p)=(VM(p)*cos(DL(p))+VM(p)*sin(DL(p))*i);
   end

fls=zeros(ld,1);
for j2=1:ld
   for i2=1:g
      fls(j2)=fls(j2)+(FLG(j2,i2)*(VC(i2)/VC(j2+g)));
      Lj(j2)=abs(1-(fls(j2)));
   end
end
Lj=Lj'

break;
else
    
%formation of jacobian matrix
%diagonal elements of J1 
for p1=2:nb
   HD(p1,p1)=-Q(p1)-(B(p1,p1)*VM(p1)*VM(p1));
end 
HD;
%off-diagonal elements of J1
%HO=zeros(nb,nb);
for p1=2:nb
   for q1=2:nb
      if q1~=p1
         HO(p1,q1)=VM(p1)*VM(q1)*(G(p1,q1)*sin(DL(p1)-DL(q1))-B(p1,q1)*cos(DL(p1)-DL(q1)));
      end
   end
end
HO;
H=[HD]+[HO];
%diagonal elements of J2
for p2=g+1:nb
   ND(p2,p2)=(P(p2)/VM(p2))+(G(p2,p2)*VM(p2));
end 
ND;
%off-diagonal elements of J2
%NO=zeros(nb,nb);
for p1=2:nb
   for q2=g+1:nb
      if q2~=p1
         NO(p1,q2)=VM(p1)*(G(p1,q2)*cos(DL(p1)-DL(q2))+B(p1,q2)*sin(DL(p1)-DL(q2)));
     end
   end
end
NO;
N=[ND]+[NO];
%diagonal elements of J3
for p2=g+1:nb
   MD(p2,p2)=P(p2)-(G(p2,p2)*VM(p2)*VM(p2));
end 
MD;
% off-diagonal elements OF J3
%MO=zeros(nb,nb);
for p2=g+1:nb
   for q1=2:nb
      if q1~=p2
          MO(p2,q1)=-VM(p2)*VM(q1)*(G(p2,q1)*cos(DL(p2)-DL(q1))+B(p2,q1)*sin(DL(p2)-DL(q1)));
       end
   end
end
MO;
M=[MD]+[MO];
%diagonal elements of J4
for p2=g+1:nb
   LD(p2,p2)=(Q(p2)/VM(p2))-(B(p2,p2)*VM(p2));
end 
LD;
% off-diagonal elements of J4
%LO=zeros(nb,nb);
for p2=g+1:nb
   for q2=g+1:nb
      if q2~=p2
         LO(p2,q2)=VM(p2)*(G(p2,q2)*sin(DL(p2)-DL(q2))-B(p2,q2)*cos(DL(p2)-DL(q2)));
       end
   end
end
LO;
L=[LD]+[LO];
J1=H(2:nb,2:nb);
J2=N(2:nb,g+1:nb);
J3=M(g+1:nb,2:nb);
J4=L(g+1:nb,g+1:nb);
J=[J1 J2;J3 J4];

PkDLk=VM(ku)*VMcr*(-Gkm*sin(DL(ku)-DLcr)+Bkm*cos(DL(ku)-DLcr))+VM(ku)*VMvr*(-Gvr*sin(DL(ku)-DLvr)+Bvr*cos(DL(ku)-DLvr));
PkDLCR=VM(ku)*VMcr*(Gkm*sin(DL(ku)-DLcr)-Bkm*cos(DL(ku)-DLcr));
PkDLVR=VM(ku)*VMvr*(Gvr*sin(DL(ku)-DLvr)-Bvr*cos(DL(ku)-DLvr));
PmDLm=VM(mu)*VMcr*(-Gmm*sin(DL(mu)-DLcr)+Bmm*cos(DL(mu)-DLcr));
PmDLCR=VM(mu)*VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));
PmkDLk=VM(ku)*VM(mu)*(Gkm*sin(DL(mu)-DL(ku))-Bkm*cos(DL(mu)-DL(ku)));
PmkDLm=VM(ku)*VM(mu)*(-Gkm*sin(DL(mu)-DL(ku))+Bkm*cos(DL(mu)-DL(ku)))+VM(mu)*VMcr*(-Gmm*sin(DL(mu)-DLcr)+Bmm*cos(DL(mu)-DLcr));
PmkDLCR=VM(mu)*VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));
PBBDLk=VM(ku)*VMcr*(Gkm*sin(DLcr-DL(ku))-Bkm*cos(DLcr-DL(ku)))+VM(ku)*VMvr*(Gvr*sin(DLvr-DL(ku))-Bvr*cos(DLvr-DL(ku)));
PBBDLm=VM(mu)*VMcr*(Gmm*sin(DLcr-DL(mu))-Bmm*cos(DLcr-DL(mu)));
PBBDLCR=VM(mu)*VMcr*(-Gmm*sin(DLcr-DL(mu))+Bmm*cos(DLcr-DL(mu)))+VM(ku)*VMcr*(-Gkm*sin(DLcr-DL(ku))+Bkm*cos(DLcr-DL(ku)));
PBBDLVR=VM(ku)*VMvr*(-Gvr*sin(DLvr-DL(ku))+Bvr*cos(DLvr-DL(ku)));


PkVk=2*VM(3)*gvr+VMcr*(Gkm*cos(DL(ku)-DLcr)+Bkm*sin(DL(ku)-DLcr))+VMvr*(Gvr*cos(DL(ku)-DLvr)+Bvr*sin(DL(ku)-DLvr));
PkVCR=VM(ku)*(Gkm*cos(DL(ku)-DLcr)+Bkm*sin(DL(ku)-DLcr));
PkVVR=VM(ku)*(Gvr*cos(DL(ku)-DLvr)+Bvr*sin(DL(ku)-DLvr));
PmVm=VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
PmVCR=VM(mu)*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
PmkVk=VM(mu)*(Gkm*cos(DL(mu)-DL(ku))+Bkm*sin(DL(mu)-DL(ku)));
PmkVm=2*VM(mu)*Gmm+VM(ku)*(Gkm*cos(DL(mu)-DL(ku))+Bkm*sin(DL(mu)-DL(ku)))+VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
PmkVCR=VM(mu)*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
PBBVk=VMcr*(Gkm*cos(DLcr-DL(ku))+Bkm*sin(DLcr-DL(ku)))+VMvr*(Gvr*cos(DLvr-DL(ku))+Bvr*sin(DLvr-DL(ku)));
PBBVm=VMcr*(Gmm*cos(DLcr-DL(mu))+Bmm*sin(DLcr-DL(mu)));
PBBVCR=2*VMcr*Gmm+VM(ku)*(Gkm*cos(DLcr-DL(ku))+Bkm*sin(DLcr-DL(ku)))+VM(mu)*(Gmm*cos(DLcr-DL(mu))+Bmm*sin(DLcr-DL(mu)));


QkDLk=VM(ku)*VMcr*(Gkm*cos(DL(ku)-DLcr)+Bkm*sin(DL(ku)-DLcr))+VM(ku)*VMvr*(Gvr*cos(DL(ku)-DLvr)+Bvr*sin(DL(ku)-DLvr));
QkDLCR=VM(ku)*VMcr*(-Gkm*cos(DL(ku)-DLcr)-Bkm*sin(DL(ku)-DLcr));
QkDLVR=VM(ku)*VMvr*(-Gvr*cos(DL(ku)-DLvr)-Bvr*sin(DL(ku)-DLvr));
QmDLm=VM(mu)*VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
QmDLCR=VM(mu)*VMcr*(-Gmm*cos(DL(mu)-DLcr)-Bmm*sin(DL(mu)-DLcr));
QmkDLk=VM(ku)*VM(mu)*(-Gkm*cos(DL(mu)-DL(ku))-Bkm*sin(DL(mu)-DL(ku)));
QmkDLm=VM(ku)*VM(mu)*(Gkm*cos(DL(mu)-DL(ku))+Bkm*sin(DL(mu)-DL(ku)))+VM(mu)*VMcr*(Gmm*cos(DL(mu)-DLcr)+Bmm*sin(DL(mu)-DLcr));
QmkDLCR=VM(mu)*VMcr*(-Gmm*cos(DL(mu)-DLcr)-Bmm*sin(DL(mu)-DLcr));


QkVk=-2*VM(ku)*bvr+VMcr*(Gkm*sin(DL(ku)-DLcr)-Bkm*cos(DL(ku)-DLcr))+VMvr*(Gvr*sin(DL(ku)-DLvr)-Bvr*cos(DL(ku)-DLvr));
QkVCR=VM(ku)*(Gkm*sin(DL(ku)-DLcr)-Bkm*cos(DL(ku)-DLcr));
QmVm=VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));
QmVCR=VM(mu)*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));
QmkVk=VM(mu)*(Gkm*sin(DL(mu)-DL(ku))-Bkm*cos(DL(mu)-DL(ku)));
QmkVm=-2*VM(mu)*Bmm+VM(ku)*(Gkm*sin(DL(mu)-DL(ku))-Bkm*cos(DL(mu)-DL(ku)))+VMcr*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));
QmkVCR=VM(mu)*(Gmm*sin(DL(mu)-DLcr)-Bmm*cos(DL(mu)-DLcr));

JM1=zeros(nb+1,nb+1);
for im1=1:nb-1
    for jm1=1:nb-1
        JM1(im1,jm1)=JM1(im1,jm1)+J1(im1,jm1);
    end
end
JM1;
JM2=zeros(nb+1,(nb-g)+1);
for im2=1:nb-1
    for jm2=1:nb-g
        JM2(im2,jm2)=JM2(im2,jm2)+J2(im2,jm2);
    end
end
JM2;
JM3=zeros((nb-g)+1,nb+1);
for im3=1:nb-g
    for jm3=1:nb-1
        JM3(im3,jm3)=JM3(im3,jm3)+J3(im3,jm3);
    end
end
JM3;
JM4=zeros(nb-g+1,nb-g+1);
for im4=1:nb-g
    for jm4=1:nb-g
        JM4(im4,jm4)=JM4(im4,jm4)+J4(im4,jm4);
    end
end
JM4;

JM1(ku-1,ku-1)=JM1(ku-1,ku-1)+PkDLk;
JM1(ku-1,nb)=JM1(ku-1,nb)+PkDLCR;
JM1(ku-1,nb+1)=JM1(ku-1,nb+1)+PkDLVR;
JM1(mu-1,mu-1)=JM1(mu-1,mu-1)+PmDLm;
JM1(mu-1,nb)=JM1(mu-1,nb)+PmDLCR;
JM1(nb,ku-1)=JM1(nb,ku-1)+PmkDLk;
JM1(nb,mu-1)=JM1(nb,mu-1)+PmkDLm;
JM1(nb,nb)=JM1(nb,nb)+PmkDLCR;
JM1(nb+1,ku-1)=JM1(nb+1,ku-1)+PBBDLk;
JM1(nb+1,mu-1)=JM1(nb+1,mu-1)+PBBDLm;
JM1(nb+1,nb)=JM1(nb+1,nb)+PBBDLCR;
JM1(nb+1,nb+1)=JM1(nb+1,nb+1)+PBBDLVR;

JM2(ku-1,ku-g)=JM2(ku-1,ku-g)+PkVk;
JM2(ku-1,nb-g+1)=JM2(ku-1,nb-g+1)+PkVCR;
JM2(mu-1,mu-g)=JM2(mu-1,mu-g)+PmVm;
JM2(mu-1,nb-g+1)=JM2(mu-1,nb-g+1)+PmVCR;
JM2(nb,ku-g)=JM2(nb,ku-g)+PmkVk;
JM2(nb,mu-g)=JM2(nb,mu-g)+PmkVm;
JM2(nb,nb-g+1)=JM2(nb,nb-g+1)+PmkVCR;
JM2(nb+1,ku-g)=JM2(nb+1,ku-g)+PBBVk;
JM2(nb+1,mu-g)=JM2(nb+1,mu-g)+PBBVm;
JM2(nb+1,nb-g+1)=JM2(nb+1,nb-g+1)+PBBVCR;

JM3(ku-g,ku-1)=JM3(ku-g,ku-1)+QkDLk;
JM3(ku-g,nb)=JM3(ku-g,nb)+QkDLCR;
JM3(ku-g,nb+1)=JM3(ku-g,nb+1)+QkDLVR;
JM3(mu-g,mu-1)=JM3(mu-g,mu-1)+QmDLm;
JM3(mu-g,nb)=JM3(mu-g,nb)+QmDLCR;
JM3(nb-g+1,ku-1)=JM3(nb-g+1,ku-1)+QmkDLk;
JM3(nb-g+1,mu-1)=JM3(nb-g+1,mu-1)+QmkDLm;
JM3(nb-g+1,nb)=JM3(nb-g+1,nb)+QmkDLCR;

JM4(ku-g,ku-g)=JM4(ku-g,ku-g)+QkVk;
JM4(ku-g,nb-g+1)=JM4(ku-g,nb-g+1)+QkVCR;
JM4(mu-g,mu-g)=JM4(mu-g,mu-g)+QmVm;
JM4(mu-g,nb-g+1)=JM4(mu-g,nb-g+1)+QmVCR;
JM4(nb-g+1,ku-g)=JM4(nb-g+1,ku-g)+QmkVk;
JM4(nb-g+1,mu-g)=JM4(nb-g+1,mu-g)+QmkVm;
JM4(nb-g+1,nb-g+1)=JM4(nb-g+1,nb-g+1)+QmkVCR;

JM1;
JM2;
JM3;
JM4;
JM=[JM1 JM2;JM3 JM4];

PV=zeros((2*nb-1-g+3),1);
PV(ku-1,1)=JPk;
PV(nb+1,1)=JPbb;
PV(nb+1+ku-g,1)=JQk;
PV;
JM((1:2*nb-1-g+3),nb+1+ku-g:nb+1+ku-g)=PV;

JM;
ddv=inv(JM)*[DPQ];

for p1=2:nb
   DL(p1)=(DL(p1)+ddv(p1-1));
   DLD(p1)=(DL(p1)/conv);
end
for p22=nb+2:2*nb-1-g+2
if p22~=nb+1+ku-g
    VM(p22-(nb-g+1))=VM(p22-(nb-g+1))+ddv(p22);
    
end
 end
end

DLcr=DLcr+ddv(nb,1);
DLvr=DLvr+ddv(nb+1,1);
DL;
DLD;
VM;
VMvr=VMvr+ddv(nb+1+ku-g);
VMcr=VMcr+ddv(2*nb-1-g+3);
itn=itn+1
end

t=etime(clock,t0)
