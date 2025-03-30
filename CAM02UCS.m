function [Japos,aaposM,baposM,h,Mapos]=CAM02UCS(SPD,Rx)

persistent SourcefileTM30

 if isempty(SourcefileTM30)
    % load data
    filenameTM30 = fullfile('source_TM30_20.csv');
    SourcefileTM30 = csvread(filenameTM30);

 end

xbar_10=SourcefileTM30(:,104);
ybar_10=SourcefileTM30(:,105);
zbar_10=SourcefileTM30(:,106);

% calculate coordinates CAM02UCS - under 1 light source only
% input is needed: SPD

% Tristimulus values 
littlek=100/sum(SPD.*ybar_10);
X=littlek*sum(SPD.*Rx.*xbar_10);
Y=littlek*sum(SPD.*Rx.*ybar_10);
Z=littlek*sum(SPD.*Rx.*zbar_10);

% Adopted white in test illuminant (like the old Xn,Yb,Zn): Xw, Yw, Zw
littlek=100/sum(SPD.*ybar_10);
Xw=littlek*sum(SPD.*xbar_10);
Yw=littlek*sum(SPD.*ybar_10);
Zw=littlek*sum(SPD.*zbar_10);

% Background in test conditions: Yb (=20 CRICAM02UCS, Fairchild)
Yb=20;
% Reference white in reference illuminant: Xwr=Ywr=Zwr=100, which are fixed in the model
Xwr=100;
Ywr=100;
Zwr=100;
% Luminance of test adapting field (cd/m2): LA (=100 CRICAM02UCS) (=60 Luo
% developments) LA=100 (TM-30-18)
LA=100;
% SURROUND PARAMETERS (Luo developments)
% Surface colors: surround=average, La=60, Sr=1, F=1.0, c=0.69, Nc=1.0
% self-luminous (display): surround=dim, La=20, Sr=0.15, F=0.9, c=0.59, Nc=0.9
% projection: surround=dark, La=30, Sr=0, F=0.8, c=0.535, Nc=0.8
Sr=1;
F=1;
c=0.69; 
Nc=1;
% D degree of adaptation. between 0-1. D=1 to discount illuminant (Fairchild,CRICAM02UCS)  
% D=F*(1-(1/3.6)*exp((-LA-42)/92));
% change D for mixed adaptation
D=1;

% Calculate RGB via CAT02 
MCAT02=[0.7328 0.4296 -0.1624;-0.7036 1.6975 0.0061;0.0030 0.0136 0.9834];

% Use Royer's inverse MCAT02 with 6 decimal points. 28 July 2021
inv_MCAT02=[1.096124 -0.278869 0.182745; 0.454369 0.473533 0.072098; -0.009628 -0.005698 1.015326];

RwGwBw=MCAT02*[Xw;Yw;Zw];
Rw=RwGwBw(1);
Gw=RwGwBw(2);
Bw=RwGwBw(3);
DR=D*Yw/Rw+1-D;
DG=D*Yw/Gw+1-D;
DB=D*Yw/Bw+1-D;
k=1/(5*LA+1);
FL=0.2*k^4*(5*LA)+0.1*(1-k^4)^2*((5*LA)^(1/3));
n=Yb/Yw;
z=1.48+sqrt(n);
Nbb=0.725*(1/n)^0.2;
Ncb=Nbb;
Rwc=DR*Rw;
Gwc=DG*Gw;
Bwc=DB*Bw;
MHPE=[0.38971 0.68898 -0.07868;-0.22981 1.18340 0.04641;0 0 1];

RaposwGaposwBaposw=MHPE*inv_MCAT02*[Rwc;Gwc;Bwc];
Raposw=RaposwGaposwBaposw(1);
Gaposw=RaposwGaposwBaposw(2);
Baposw=RaposwGaposwBaposw(3);
Raposaw=400*(((FL*Raposw/100)^0.42)/((FL*Raposw/100)^0.42+27.13))+0.1;
Gaposaw=400*(((FL*Gaposw/100)^0.42)/((FL*Gaposw/100)^0.42+27.13))+0.1;
Baposaw=400*(((FL*Baposw/100)^0.42)/((FL*Baposw/100)^0.42+27.13))+0.1;
Aw=(2*Raposaw+Gaposaw+Baposaw/20-0.305)*Nbb;

% Calculate (sharpened) cone responses (transfer CMFs to sharper sensors)
RGB=MCAT02*[X;Y;Z];
R=RGB(1);
G=RGB(2);
B=RGB(3);

% Calculate the corresponding (sharpened) cone response
Rc=DR*R;
Gc=DG*G;
Bc=DB*B;

% calculate the Hunt-Pointer-Estevez response
RaposGaposBapos=MHPE*inv_MCAT02*[Rc;Gc;Bc];
Rapos=RaposGaposBapos(1);
Gapos=RaposGaposBapos(2);
Bapos=RaposGaposBapos(3);

% post-adaptation cone response 
if Rapos<0
  Raposa=-400*(((-FL*Rapos/100)^0.42)/((-FL*Rapos/100)^0.42+27.13))+0.1;
else
  Raposa=400*(((FL*Rapos/100)^0.42)/((FL*Rapos/100)^0.42+27.13))+0.1;  
end

if Gapos<0
  Gaposa=-400*(((-FL*Gapos/100)^0.42)/((-FL*Gapos/100)^0.42+27.13))+0.1;
else
  Gaposa=400*(((FL*Gapos/100)^0.42)/((FL*Gapos/100)^0.42+27.13))+0.1;
end

if Bapos<0
  Baposa=-400*(((-FL*Bapos/100)^0.42)/((-FL*Bapos/100)^0.42+27.13))+0.1;
else
  Baposa=400*(((FL*Bapos/100)^0.42)/((FL*Bapos/100)^0.42+27.13))+0.1;
 end

% Added _ after a and b, so that fit_ellipse function can use a and b as
% output parameteres 
a_=Raposa-12*Gaposa/11+Baposa/11;
b_=(Raposa+Gaposa-2*Baposa)/9;

% calculate (h) hue angle 
if (a_>0) && (b_>0) 
    h=atand(b_/a_);
elseif (a_>0) && (b_<0) 
    h=atand(b_/a_)+360;  
else
    h=atand(b_/a_)+180;  
end

% et eccentricity, H hue composition
if (h<20.14)
   hapos=h+360;
else
   hapos=h;
end

if (hapos>=20.14) && (hapos<90)
    et=(cos(hapos*pi/180+2)+3.8)/4;
    H=(100*(hapos-20.14)/0.8)/(((hapos-20.14)/0.8)+((90-hapos)/0.7));
elseif (hapos>=90) && (hapos<164.25)
    et=(cos(hapos*pi/180+2)+3.8)/4;
    H=100+(100*(hapos-90)/0.7)/(((hapos-90)/0.7)+((164.25-hapos)/1));    
elseif (hapos>=164.25) && (hapos<237.53)
    et=(cos(hapos*pi/180+2)+3.8)/4;
    H=200+(100*(hapos-164.25)/1)/(((hapos-164.25)/1)+((237.53-hapos)/1.2));     
elseif (hapos>=237.53) && (hapos<380.14)
    et=(cos(hapos*pi/180+2)+3.8)/4;
    H=300+(100*(hapos-237.53)/1.2)/(((hapos-237.53)/1.2)+((380.14-hapos)/0.8));
else
    hapos
    warning('et and H not defined')
    
end 
    
% achromatic response A
A=(2*Raposa+Gaposa+Baposa/20-0.305)*Nbb;  
% Lightiness J
J=100*(A/Aw)^(c*z);
% brightness Q
Q=(4/c)*sqrt(J/100)*(Aw+4)*FL^0.25;
% chroma C, colourfulness M, saturation s 
t=(50000/13*Nc*Ncb*et*sqrt(a_^2+b_^2))/(Raposa+Gaposa+(21/20)*Baposa);
C=(t^0.9)*sqrt(J/100)*(1.64-0.29^n)^0.73;
M=C*FL^0.25;
s=100*sqrt(M/Q);


% CIECAM02UCS coordinates 
Japos=(1.7*J)/(1+0.007*J);
% Ronnie Luo prefers (1/0.228) instead of 43.86 (CIE)
Mapos=(1/0.0228)*log(1+0.0228*M);
aaposM=Mapos*cosd(h);
baposM=Mapos*sind(h);

