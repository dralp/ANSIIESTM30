function [R_f, R_g, CCT, Duv,TM30_struct]=TM_30_20(SPDin)
%[R_f, R_g, T_x, Duv_test,TM30_struct]=TM_30_20(SPDin)
%
% Calculate IES-TM-30-20, Rf, Rg, CCT, and Duv
% IES TM-30: https://store.ies.org/product/tm-30-20-ies-method-for-evaluating-light-source-color-rendition/
%
% User INPUT needed: SPDin (spectrum of test light source, wavelength in nm in first column, spectral distribution in second)
% E.g., For equal-energy radiator: SPDin=[(360:830)',ones(471,1)];
%Example:
%
%  SPDin=[(360:830)',ones(471,1)]
% 
%  [R_f, R_g, T_x, Duv_test,TM30_struct]=TM_30_20(SPDin)




persistent SourcefileTM30


%load data file
if isempty(SourcefileTM30)
    % load data
    filenameTM30 = fullfile('source_TM30_20.csv');
    SourcefileTM30 = csvread(filenameTM30);

end

if nargin()<1
    SPDin=loadStandardLEDspectra;
end

SPDin=cropSPD(SPDin,380,780);

SPDin=interpSPD1onSPD2(SPDin,[(380:780)',ones(401,1)]);

Stest1=SPDin(:,2);

wavelength=SourcefileTM30(:,100);
%CIE 2 and 10 degree observer
xbar=SourcefileTM30(:,101);
ybar=SourcefileTM30(:,102);
zbar=SourcefileTM30(:,103);
xbar_10=SourcefileTM30(:,104);
ybar_10=SourcefileTM30(:,105);
zbar_10=SourcefileTM30(:,106);
% constants from CIE15:2004
S_0=SourcefileTM30(:,107);
S_1=SourcefileTM30(:,108);
S_2=SourcefileTM30(:,109);
% Rescale test SPD so SPDxVlamd=100
Stest=Stest1*(100/sum(Stest1.*ybar_10));

%% generate the reference illuminant

% calculate test SPD CCT and Duv (Tt) according to Ohno CCT Duv Leukos
[T_x,Duv_test]=CCT_Duv(Stest1);
T_t=T_x;
if T_t <= 4000
    % Planckian radiation constants per IES TM-30-20 Excel 2.04
    c2=0.014388;
    c1=3.7415*10^-16;

    % Le,lamnda numerator
    for kk=1:401
        Le_num(kk)=(c1/((wavelength(kk,1)*0.000000001)^5)/(exp(c2/((wavelength(kk,1)*0.000000001)*T_t))-1));
    end

    % Normalize for Y10=100
    Sref=Le_num'*(100/sum(Le_num'.*ybar_10));

elseif T_t >= 5000

    % Daylight
    if T_t <= 7000
        x_D=((-4.6070*10^9)/T_t^3)+((2.9678*10^6)/T_t^2)+((0.09911*10^3)/T_t)+0.244063;
        y_D=-3*(x_D^2)+2.87*x_D-0.275;
    else
        x_D=((-2.0064*10^9)/T_t^3)+((1.9018*10^6)/T_t^2)+((0.24748*10^3)/T_t)+0.23704;
        y_D=-3*(x_D^2)+2.87*x_D-0.275;
    end
    M_1=(-1.3515-1.7703*x_D+5.9114*y_D)/(0.0241+0.2562*x_D-0.7341*y_D);
    M_2=(0.03-31.4424*x_D+30.0717*y_D)/(0.0241+0.2562*x_D-0.7341*y_D);
    Sref=S_0+(M_1*S_1)+(M_2*S_2);

    % Rescale ref SPD so SPDxVlamd=100
    Sref=Sref*(100/sum(Sref.*ybar_10));

else
    % a mix of Plackian and Daylight

    % Planckian radiation constants per IES TM-30-20 Excel 2.04
    c2=0.014388;
    c1=3.7415*10^-16;

    % Le,lamnda numerator
    for kk=1:401
        Le_num(kk)=(c1/((wavelength(kk,1)*0.000000001)^5)/(exp(c2/((wavelength(kk,1)*0.000000001)*T_t))-1));
    end

    % Rescale ref SPD so SPDxVlamd=100
    Sref_P=Le_num'*(100/sum(Le_num'.*ybar_10));


    % Daylight
    if T_t <= 7000
        x_D=((-4.6070*10^9)/T_t^3)+((2.9678*10^6)/T_t^2)+((0.09911*10^3)/T_t)+0.244063;
        y_D=-3*(x_D^2)+2.87*x_D-0.275;
    else
        x_D=((-2.0064*10^9)/T_t^3)+((1.9018*10^6)/T_t^2)+((0.24748*10^3)/T_t)+0.23704;
        y_D=-3*(x_D^2)+2.87*x_D-0.275;
    end
    M_1=(-1.3515-1.7703*x_D+5.9114*y_D)/(0.0241+0.2562*x_D-0.7341*y_D);
    M_2=(0.03-31.4424*x_D+30.0717*y_D)/(0.0241+0.2562*x_D-0.7341*y_D);
    Sref_D=S_0+(M_1*S_1)+(M_2*S_2);

    % normalize for Y10=100
    Sref_D2=Sref_D*(100/sum(Sref_D.*ybar_10));

    % mixed reference illuminant
    Sref=((5000-T_t)/1000)*Sref_P+(1-((5000-T_t)/1000))*Sref_D2;

    % Rescale ref SPD so SPDxVlamd=100
    Sref=Sref*(100/sum(Sref.*ybar_10));
end

CCT=round(T_x);
Duv=round(Duv_test,4);


%%

% reference illuminant tristimulus values
k_r=100/sum(Sref.*ybar_10);
X_10_r=k_r*sum(Sref.*xbar_10);
Y_10_r=k_r*sum(Sref.*ybar_10);
Z_10_r=k_r*sum(Sref.*zbar_10);

% calculate deltaE for 99 test samples

for iii=1:99
    Rx=SourcefileTM30(:,iii);

    % calculate coordinates for reference illuminant
    SPD1=Sref;
    [Japos_ref(iii),aaposM_ref(iii),baposM_ref(iii)]=CAM02UCS(SPD1,Rx);


    % calculate coordinates for test SPD
    SPD1=Stest;
    [Japos_test(iii),aaposM_test(iii),baposM_test(iii)]=CAM02UCS(SPD1,Rx);

    % calculate colour difference in CIECAM02UCS
    deltaE_CAM02(:,iii)=sqrt((Japos_test(iii)-Japos_ref(iii))^2+(aaposM_test(iii)-aaposM_ref(iii))^2+(baposM_test(iii)-baposM_ref(iii))^2);
    ab_test(iii,:)=[aaposM_test(iii),baposM_test(iii)];
    ab_ref(iii,:)=[aaposM_ref(iii),baposM_ref(iii)];
end

% Rf fidelity index
R_f_apos=100-6.73*mean(deltaE_CAM02);
R_f=10*log(exp(R_f_apos/10)+1);

% Define hue bins
huebin(1,1)=polyshape([0 1000 1000],[0 0 414.21]);
huebin(2,1)=polyshape([0 1000 1000],[0 414.21 1000]);
huebin(3,1)=polyshape([0 1000 414.21],[0 1000 1000]);
huebin(4,1)=polyshape([0 414.21 0],[0 1000 1000]);
huebin(5,1)=polyshape([0 -414.21 0],[0 1000 1000]);
huebin(6,1)=polyshape([0 -1000 -414.21],[0 1000 1000]);
huebin(7,1)=polyshape([0 -1000 -1000],[0 414.21 1000]);
huebin(8,1)=polyshape([0 -1000 -1000],[0 0 414.21]);
huebin(9,1)=polyshape([0 -1000 -1000],[0 0 -414.21]);
huebin(10,1)=polyshape([0 -1000 -1000],[0 -414.21 -1000]);
huebin(11,1)=polyshape([0 -1000 -414.21],[0 -1000 -1000]);
huebin(12,1)=polyshape([0 -414.21 0],[0 -1000 -1000]);
huebin(13,1)=polyshape([0 0 414.21],[0 -1000 -1000]);
huebin(14,1)=polyshape([0 414.21 1000],[0 -1000 -1000]);
huebin(15,1)=polyshape([0 1000 1000],[0 -1000 -414.21]);
huebin(16,1)=polyshape([0 1000 1000],[0 -414.21 0]);

% Identify hue bins for each sample (a-b coordinate)

% we could simplyfi the code by doing a loop like below. Having numbered
% variables is a bit of an 
% for n=1:16
%     huebin_ID{n}=find(isinterior(huebin(n,1),ab_ref(:,1),ab_ref(:,2))==1);
% end


huebin_ID_1= find(isinterior(huebin(1,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_2= find(isinterior(huebin(2,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_3= find(isinterior(huebin(3,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_4= find(isinterior(huebin(4,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_5= find(isinterior(huebin(5,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_6= find(isinterior(huebin(6,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_7= find(isinterior(huebin(7,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_8= find(isinterior(huebin(8,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_9= find(isinterior(huebin(9,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_10=find(isinterior(huebin(10,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_11=find(isinterior(huebin(11,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_12=find(isinterior(huebin(12,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_13=find(isinterior(huebin(13,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_14=find(isinterior(huebin(14,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_15=find(isinterior(huebin(15,1),ab_ref(:,1),ab_ref(:,2))==1);
huebin_ID_16=find(isinterior(huebin(16,1),ab_ref(:,1),ab_ref(:,2))==1);


% for n=1:16
%     for hh=1:numel(huebin_ID_1)
%         hb_test{n}(hh,:)=[ab_test((huebin_ID{n}(hh,1)),1),ab_test((huebin_ID{n}(hh,1)),2)];
%         hb_ref(hh,:)=[ab_ref((huebin_ID{n}(hh,1)),1),ab_ref((huebin_ID{n}(hh,1)),2)];
%     end
% end


% write sample a-b values to a matrix. later to average
for hh=1:numel(huebin_ID_1)
    hb_1_test(hh,:)=[ab_test((huebin_ID_1(hh,1)),1),ab_test((huebin_ID_1(hh,1)),2)];
    hb_1_ref(hh,:)=[ab_ref((huebin_ID_1(hh,1)),1),ab_ref((huebin_ID_1(hh,1)),2)];
end
for hh=1:numel(huebin_ID_2)
    hb_2_test(hh,:)=[ab_test((huebin_ID_2(hh,1)),1),ab_test((huebin_ID_2(hh,1)),2)];
    hb_2_ref(hh,:)=[ab_ref((huebin_ID_2(hh,1)),1),ab_ref((huebin_ID_2(hh,1)),2)];
end
for hh=1:numel(huebin_ID_3)
    hb_3_test(hh,:)=[ab_test((huebin_ID_3(hh,1)),1),ab_test((huebin_ID_3(hh,1)),2)];
    hb_3_ref(hh,:)=[ab_ref((huebin_ID_3(hh,1)),1),ab_ref((huebin_ID_3(hh,1)),2)];
end
for hh=1:numel(huebin_ID_4)
    hb_4_test(hh,:)=[ab_test((huebin_ID_4(hh,1)),1),ab_test((huebin_ID_4(hh,1)),2)];
    hb_4_ref(hh,:)=[ab_ref((huebin_ID_4(hh,1)),1),ab_ref((huebin_ID_4(hh,1)),2)];
end
for hh=1:numel(huebin_ID_5)
    hb_5_test(hh,:)=[ab_test((huebin_ID_5(hh,1)),1),ab_test((huebin_ID_5(hh,1)),2)];
    hb_5_ref(hh,:)=[ab_ref((huebin_ID_5(hh,1)),1),ab_ref((huebin_ID_5(hh,1)),2)];
end
for hh=1:numel(huebin_ID_6)
    hb_6_test(hh,:)=[ab_test((huebin_ID_6(hh,1)),1),ab_test((huebin_ID_6(hh,1)),2)];
    hb_6_ref(hh,:)=[ab_ref((huebin_ID_6(hh,1)),1),ab_ref((huebin_ID_6(hh,1)),2)];
end
for hh=1:numel(huebin_ID_7)
    hb_7_test(hh,:)=[ab_test((huebin_ID_7(hh,1)),1),ab_test((huebin_ID_7(hh,1)),2)];
    hb_7_ref(hh,:)=[ab_ref((huebin_ID_7(hh,1)),1),ab_ref((huebin_ID_7(hh,1)),2)];
end
for hh=1:numel(huebin_ID_8)
    hb_8_test(hh,:)=[ab_test((huebin_ID_8(hh,1)),1),ab_test((huebin_ID_8(hh,1)),2)];
    hb_8_ref(hh,:)=[ab_ref((huebin_ID_8(hh,1)),1),ab_ref((huebin_ID_8(hh,1)),2)];
end
for hh=1:numel(huebin_ID_9)
    hb_9_test(hh,:)=[ab_test((huebin_ID_9(hh,1)),1),ab_test((huebin_ID_9(hh,1)),2)];
    hb_9_ref(hh,:)=[ab_ref((huebin_ID_9(hh,1)),1),ab_ref((huebin_ID_9(hh,1)),2)];
end
for hh=1:numel(huebin_ID_10)
    hb_10_test(hh,:)=[ab_test((huebin_ID_10(hh,1)),1),ab_test((huebin_ID_10(hh,1)),2)];
    hb_10_ref(hh,:)=[ab_ref((huebin_ID_10(hh,1)),1),ab_ref((huebin_ID_10(hh,1)),2)];
end
for hh=1:numel(huebin_ID_11)
    hb_11_test(hh,:)=[ab_test((huebin_ID_11(hh,1)),1),ab_test((huebin_ID_11(hh,1)),2)];
    hb_11_ref(hh,:)=[ab_ref((huebin_ID_11(hh,1)),1),ab_ref((huebin_ID_11(hh,1)),2)];
end
for hh=1:numel(huebin_ID_12)
    hb_12_test(hh,:)=[ab_test((huebin_ID_12(hh,1)),1),ab_test((huebin_ID_12(hh,1)),2)];
    hb_12_ref(hh,:)=[ab_ref((huebin_ID_12(hh,1)),1),ab_ref((huebin_ID_12(hh,1)),2)];
end
for hh=1:numel(huebin_ID_13)
    hb_13_test(hh,:)=[ab_test((huebin_ID_13(hh,1)),1),ab_test((huebin_ID_13(hh,1)),2)];
    hb_13_ref(hh,:)=[ab_ref((huebin_ID_13(hh,1)),1),ab_ref((huebin_ID_13(hh,1)),2)];
end
for hh=1:numel(huebin_ID_14)
    hb_14_test(hh,:)=[ab_test((huebin_ID_14(hh,1)),1),ab_test((huebin_ID_14(hh,1)),2)];
    hb_14_ref(hh,:)=[ab_ref((huebin_ID_14(hh,1)),1),ab_ref((huebin_ID_14(hh,1)),2)];
end
for hh=1:numel(huebin_ID_15)
    hb_15_test(hh,:)=[ab_test((huebin_ID_15(hh,1)),1),ab_test((huebin_ID_15(hh,1)),2)];
    hb_15_ref(hh,:)=[ab_ref((huebin_ID_15(hh,1)),1),ab_ref((huebin_ID_15(hh,1)),2)];
end
for hh=1:numel(huebin_ID_16)
    hb_16_test(hh,:)=[ab_test((huebin_ID_16(hh,1)),1),ab_test((huebin_ID_16(hh,1)),2)];
    hb_16_ref(hh,:)=[ab_ref((huebin_ID_16(hh,1)),1),ab_ref((huebin_ID_16(hh,1)),2)];
end





% average a-b values for each hue bin

% The long lines of code could be done like this
% for n=1:16
%     gamut_test(:,n)=[mean(hb_test{n}(:,1)),mean(hb_test{n}(:,2))];
%     gamut_ref(:,n)=[mean(hb_ref{n}(:,1)),mean(hb_ref{n}(:,2))];
% end

gamut_test=[mean(hb_1_test(:,1)),mean(hb_1_test(:,2));mean(hb_2_test(:,1)),mean(hb_2_test(:,2));mean(hb_3_test(:,1)),mean(hb_3_test(:,2));mean(hb_4_test(:,1)),mean(hb_4_test(:,2));mean(hb_5_test(:,1)),mean(hb_5_test(:,2));mean(hb_6_test(:,1)),mean(hb_6_test(:,2));mean(hb_7_test(:,1)),mean(hb_7_test(:,2));mean(hb_8_test(:,1)),mean(hb_8_test(:,2));mean(hb_9_test(:,1)),mean(hb_9_test(:,2));mean(hb_10_test(:,1)),mean(hb_10_test(:,2));mean(hb_11_test(:,1)),mean(hb_11_test(:,2));mean(hb_12_test(:,1)),mean(hb_12_test(:,2));mean(hb_13_test(:,1)),mean(hb_13_test(:,2));mean(hb_14_test(:,1)),mean(hb_14_test(:,2));mean(hb_15_test(:,1)),mean(hb_15_test(:,2));mean(hb_16_test(:,1)),mean(hb_16_test(:,2))];
gamut_ref=[mean(hb_1_ref(:,1)),mean(hb_1_ref(:,2));mean(hb_2_ref(:,1)),mean(hb_2_ref(:,2));mean(hb_3_ref(:,1)),mean(hb_3_ref(:,2));mean(hb_4_ref(:,1)),mean(hb_4_ref(:,2));mean(hb_5_ref(:,1)),mean(hb_5_ref(:,2));mean(hb_6_ref(:,1)),mean(hb_6_ref(:,2));mean(hb_7_ref(:,1)),mean(hb_7_ref(:,2));mean(hb_8_ref(:,1)),mean(hb_8_ref(:,2));mean(hb_9_ref(:,1)),mean(hb_9_ref(:,2));mean(hb_10_ref(:,1)),mean(hb_10_ref(:,2));mean(hb_11_ref(:,1)),mean(hb_11_ref(:,2));mean(hb_12_ref(:,1)),mean(hb_12_ref(:,2));mean(hb_13_ref(:,1)),mean(hb_13_ref(:,2));mean(hb_14_ref(:,1)),mean(hb_14_ref(:,2));mean(hb_15_ref(:,1)),mean(hb_15_ref(:,2));mean(hb_16_ref(:,1)),mean(hb_16_ref(:,2))];
% Rg gamut index
R_g=100*(area(polyshape(gamut_test))/area(polyshape(gamut_ref)));

%Rcs,hj Local chroma shift
R_cs_h(1,:)=((cosd(11.25)*(mean(hb_1_test(:,1))-mean(hb_1_ref(:,1)))/(sqrt(mean(hb_1_ref(:,1))^2+mean(hb_1_ref(:,2))^2)))+(sind(11.25)*(mean(hb_1_test(:,2))-mean(hb_1_ref(:,2)))/(sqrt(mean(hb_1_ref(:,1))^2+mean(hb_1_ref(:,2))^2))));
R_cs_h(2,:)=((cosd(33.75)*(mean(hb_2_test(:,1))-mean(hb_2_ref(:,1)))/(sqrt(mean(hb_2_ref(:,1))^2+mean(hb_2_ref(:,2))^2)))+(sind(33.75)*(mean(hb_2_test(:,2))-mean(hb_2_ref(:,2)))/(sqrt(mean(hb_2_ref(:,1))^2+mean(hb_2_ref(:,2))^2))));
R_cs_h(3,:)=((cosd(56.25)*(mean(hb_3_test(:,1))-mean(hb_3_ref(:,1)))/(sqrt(mean(hb_3_ref(:,1))^2+mean(hb_3_ref(:,2))^2)))+(sind(56.25)*(mean(hb_3_test(:,2))-mean(hb_3_ref(:,2)))/(sqrt(mean(hb_3_ref(:,1))^2+mean(hb_3_ref(:,2))^2))));
R_cs_h(4,:)=((cosd(78.75)*(mean(hb_4_test(:,1))-mean(hb_4_ref(:,1)))/(sqrt(mean(hb_4_ref(:,1))^2+mean(hb_4_ref(:,2))^2)))+(sind(78.75)*(mean(hb_4_test(:,2))-mean(hb_4_ref(:,2)))/(sqrt(mean(hb_4_ref(:,1))^2+mean(hb_4_ref(:,2))^2))));
R_cs_h(5,:)=((cosd(101.25)*(mean(hb_5_test(:,1))-mean(hb_5_ref(:,1)))/(sqrt(mean(hb_5_ref(:,1))^2+mean(hb_5_ref(:,2))^2)))+(sind(101.25)*(mean(hb_5_test(:,2))-mean(hb_5_ref(:,2)))/(sqrt(mean(hb_5_ref(:,1))^2+mean(hb_5_ref(:,2))^2))));
R_cs_h(6,:)=((cosd(123.75)*(mean(hb_6_test(:,1))-mean(hb_6_ref(:,1)))/(sqrt(mean(hb_6_ref(:,1))^2+mean(hb_6_ref(:,2))^2)))+(sind(123.75)*(mean(hb_6_test(:,2))-mean(hb_6_ref(:,2)))/(sqrt(mean(hb_6_ref(:,1))^2+mean(hb_6_ref(:,2))^2))));
R_cs_h(7,:)=((cosd(146.25)*(mean(hb_7_test(:,1))-mean(hb_7_ref(:,1)))/(sqrt(mean(hb_7_ref(:,1))^2+mean(hb_7_ref(:,2))^2)))+(sind(146.25)*(mean(hb_7_test(:,2))-mean(hb_7_ref(:,2)))/(sqrt(mean(hb_7_ref(:,1))^2+mean(hb_7_ref(:,2))^2))));
R_cs_h(8,:)=((cosd(168.75)*(mean(hb_8_test(:,1))-mean(hb_8_ref(:,1)))/(sqrt(mean(hb_8_ref(:,1))^2+mean(hb_8_ref(:,2))^2)))+(sind(168.75)*(mean(hb_8_test(:,2))-mean(hb_8_ref(:,2)))/(sqrt(mean(hb_8_ref(:,1))^2+mean(hb_8_ref(:,2))^2))));
R_cs_h(9,:)=((cosd(191.25)*(mean(hb_9_test(:,1))-mean(hb_9_ref(:,1)))/(sqrt(mean(hb_9_ref(:,1))^2+mean(hb_9_ref(:,2))^2)))+(sind(191.25)*(mean(hb_9_test(:,2))-mean(hb_9_ref(:,2)))/(sqrt(mean(hb_9_ref(:,1))^2+mean(hb_9_ref(:,2))^2))));
R_cs_h(10,:)=((cosd(213.75)*(mean(hb_10_test(:,1))-mean(hb_10_ref(:,1)))/(sqrt(mean(hb_10_ref(:,1))^2+mean(hb_10_ref(:,2))^2)))+(sind(213.75)*(mean(hb_10_test(:,2))-mean(hb_10_ref(:,2)))/(sqrt(mean(hb_10_ref(:,1))^2+mean(hb_10_ref(:,2))^2))));
R_cs_h(11,:)=((cosd(236.25)*(mean(hb_11_test(:,1))-mean(hb_11_ref(:,1)))/(sqrt(mean(hb_11_ref(:,1))^2+mean(hb_11_ref(:,2))^2)))+(sind(236.25)*(mean(hb_11_test(:,2))-mean(hb_11_ref(:,2)))/(sqrt(mean(hb_11_ref(:,1))^2+mean(hb_11_ref(:,2))^2))));
R_cs_h(12,:)=((cosd(258.75)*(mean(hb_12_test(:,1))-mean(hb_12_ref(:,1)))/(sqrt(mean(hb_12_ref(:,1))^2+mean(hb_12_ref(:,2))^2)))+(sind(258.75)*(mean(hb_12_test(:,2))-mean(hb_12_ref(:,2)))/(sqrt(mean(hb_12_ref(:,1))^2+mean(hb_12_ref(:,2))^2))));
R_cs_h(13,:)=((cosd(281.25)*(mean(hb_13_test(:,1))-mean(hb_13_ref(:,1)))/(sqrt(mean(hb_13_ref(:,1))^2+mean(hb_13_ref(:,2))^2)))+(sind(281.25)*(mean(hb_13_test(:,2))-mean(hb_13_ref(:,2)))/(sqrt(mean(hb_13_ref(:,1))^2+mean(hb_13_ref(:,2))^2))));
R_cs_h(14,:)=((cosd(303.75)*(mean(hb_14_test(:,1))-mean(hb_14_ref(:,1)))/(sqrt(mean(hb_14_ref(:,1))^2+mean(hb_14_ref(:,2))^2)))+(sind(303.75)*(mean(hb_14_test(:,2))-mean(hb_14_ref(:,2)))/(sqrt(mean(hb_14_ref(:,1))^2+mean(hb_14_ref(:,2))^2))));
R_cs_h(15,:)=((cosd(326.25)*(mean(hb_15_test(:,1))-mean(hb_15_ref(:,1)))/(sqrt(mean(hb_15_ref(:,1))^2+mean(hb_15_ref(:,2))^2)))+(sind(326.25)*(mean(hb_15_test(:,2))-mean(hb_15_ref(:,2)))/(sqrt(mean(hb_15_ref(:,1))^2+mean(hb_15_ref(:,2))^2))));
R_cs_h(16,:)=((cosd(348.75)*(mean(hb_16_test(:,1))-mean(hb_16_ref(:,1)))/(sqrt(mean(hb_16_ref(:,1))^2+mean(hb_16_ref(:,2))^2)))+(sind(348.75)*(mean(hb_16_test(:,2))-mean(hb_16_ref(:,2)))/(sqrt(mean(hb_16_ref(:,1))^2+mean(hb_16_ref(:,2))^2))));

%Rhs,hj Local hue shift
R_hs_h(1,:)=-(sind(11.25)*(mean(hb_1_test(:,1))-mean(hb_1_ref(:,1)))/(sqrt(mean(hb_1_ref(:,1))^2+mean(hb_1_ref(:,2))^2)))+(cosd(11.25)*(mean(hb_1_test(:,2))-mean(hb_1_ref(:,2)))/(sqrt(mean(hb_1_ref(:,1))^2+mean(hb_1_ref(:,2))^2)));
R_hs_h(2,:)=-(sind(33.75)*(mean(hb_2_test(:,1))-mean(hb_2_ref(:,1)))/(sqrt(mean(hb_2_ref(:,1))^2+mean(hb_2_ref(:,2))^2)))+(cosd(33.75)*(mean(hb_2_test(:,2))-mean(hb_2_ref(:,2)))/(sqrt(mean(hb_2_ref(:,1))^2+mean(hb_2_ref(:,2))^2)));
R_hs_h(3,:)=-(sind(56.25)*(mean(hb_3_test(:,1))-mean(hb_3_ref(:,1)))/(sqrt(mean(hb_3_ref(:,1))^2+mean(hb_3_ref(:,2))^2)))+(cosd(56.25)*(mean(hb_3_test(:,2))-mean(hb_3_ref(:,2)))/(sqrt(mean(hb_3_ref(:,1))^2+mean(hb_3_ref(:,2))^2)));
R_hs_h(4,:)=-(sind(78.75)*(mean(hb_4_test(:,1))-mean(hb_4_ref(:,1)))/(sqrt(mean(hb_4_ref(:,1))^2+mean(hb_4_ref(:,2))^2)))+(cosd(78.75)*(mean(hb_4_test(:,2))-mean(hb_4_ref(:,2)))/(sqrt(mean(hb_4_ref(:,1))^2+mean(hb_4_ref(:,2))^2)));
R_hs_h(5,:)=-(sind(101.25)*(mean(hb_5_test(:,1))-mean(hb_5_ref(:,1)))/(sqrt(mean(hb_5_ref(:,1))^2+mean(hb_5_ref(:,2))^2)))+(cosd(101.25)*(mean(hb_5_test(:,2))-mean(hb_5_ref(:,2)))/(sqrt(mean(hb_5_ref(:,1))^2+mean(hb_5_ref(:,2))^2)));
R_hs_h(6,:)=-(sind(123.75)*(mean(hb_6_test(:,1))-mean(hb_6_ref(:,1)))/(sqrt(mean(hb_6_ref(:,1))^2+mean(hb_6_ref(:,2))^2)))+(cosd(123.75)*(mean(hb_6_test(:,2))-mean(hb_6_ref(:,2)))/(sqrt(mean(hb_6_ref(:,1))^2+mean(hb_6_ref(:,2))^2)));
R_hs_h(7,:)=-(sind(146.25)*(mean(hb_7_test(:,1))-mean(hb_7_ref(:,1)))/(sqrt(mean(hb_7_ref(:,1))^2+mean(hb_7_ref(:,2))^2)))+(cosd(146.25)*(mean(hb_7_test(:,2))-mean(hb_7_ref(:,2)))/(sqrt(mean(hb_7_ref(:,1))^2+mean(hb_7_ref(:,2))^2)));
R_hs_h(8,:)=-(sind(168.75)*(mean(hb_8_test(:,1))-mean(hb_8_ref(:,1)))/(sqrt(mean(hb_8_ref(:,1))^2+mean(hb_8_ref(:,2))^2)))+(cosd(168.75)*(mean(hb_8_test(:,2))-mean(hb_8_ref(:,2)))/(sqrt(mean(hb_8_ref(:,1))^2+mean(hb_8_ref(:,2))^2)));
R_hs_h(9,:)=-(sind(191.25)*(mean(hb_9_test(:,1))-mean(hb_9_ref(:,1)))/(sqrt(mean(hb_9_ref(:,1))^2+mean(hb_9_ref(:,2))^2)))+(cosd(191.25)*(mean(hb_9_test(:,2))-mean(hb_9_ref(:,2)))/(sqrt(mean(hb_9_ref(:,1))^2+mean(hb_9_ref(:,2))^2)));
R_hs_h(10,:)=-(sind(213.75)*(mean(hb_10_test(:,1))-mean(hb_10_ref(:,1)))/(sqrt(mean(hb_10_ref(:,1))^2+mean(hb_10_ref(:,2))^2)))+(cosd(213.75)*(mean(hb_10_test(:,2))-mean(hb_10_ref(:,2)))/(sqrt(mean(hb_10_ref(:,1))^2+mean(hb_10_ref(:,2))^2)));
R_hs_h(11,:)=-(sind(236.25)*(mean(hb_11_test(:,1))-mean(hb_11_ref(:,1)))/(sqrt(mean(hb_11_ref(:,1))^2+mean(hb_11_ref(:,2))^2)))+(cosd(236.25)*(mean(hb_11_test(:,2))-mean(hb_11_ref(:,2)))/(sqrt(mean(hb_11_ref(:,1))^2+mean(hb_11_ref(:,2))^2)));
R_hs_h(12,:)=-(sind(258.75)*(mean(hb_12_test(:,1))-mean(hb_12_ref(:,1)))/(sqrt(mean(hb_12_ref(:,1))^2+mean(hb_12_ref(:,2))^2)))+(cosd(258.75)*(mean(hb_12_test(:,2))-mean(hb_12_ref(:,2)))/(sqrt(mean(hb_12_ref(:,1))^2+mean(hb_12_ref(:,2))^2)));
R_hs_h(13,:)=-(sind(281.25)*(mean(hb_13_test(:,1))-mean(hb_13_ref(:,1)))/(sqrt(mean(hb_13_ref(:,1))^2+mean(hb_13_ref(:,2))^2)))+(cosd(281.25)*(mean(hb_13_test(:,2))-mean(hb_13_ref(:,2)))/(sqrt(mean(hb_13_ref(:,1))^2+mean(hb_13_ref(:,2))^2)));
R_hs_h(14,:)=-(sind(303.75)*(mean(hb_14_test(:,1))-mean(hb_14_ref(:,1)))/(sqrt(mean(hb_14_ref(:,1))^2+mean(hb_14_ref(:,2))^2)))+(cosd(303.75)*(mean(hb_14_test(:,2))-mean(hb_14_ref(:,2)))/(sqrt(mean(hb_14_ref(:,1))^2+mean(hb_14_ref(:,2))^2)));
R_hs_h(15,:)=-(sind(326.25)*(mean(hb_15_test(:,1))-mean(hb_15_ref(:,1)))/(sqrt(mean(hb_15_ref(:,1))^2+mean(hb_15_ref(:,2))^2)))+(cosd(326.25)*(mean(hb_15_test(:,2))-mean(hb_15_ref(:,2)))/(sqrt(mean(hb_15_ref(:,1))^2+mean(hb_15_ref(:,2))^2)));
R_hs_h(16,:)=-(sind(348.75)*(mean(hb_16_test(:,1))-mean(hb_16_ref(:,1)))/(sqrt(mean(hb_16_ref(:,1))^2+mean(hb_16_ref(:,2))^2)))+(cosd(348.75)*(mean(hb_16_test(:,2))-mean(hb_16_ref(:,2)))/(sqrt(mean(hb_16_ref(:,1))^2+mean(hb_16_ref(:,2))^2)));




%Rf,hj Local color fidelity
R_f_h_apos(1,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_1));
R_f_h(1,:)=10*log(exp(R_f_h_apos(1,:)/10)+1);
R_f_h_apos(2,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_2));
R_f_h(2,:)=10*log(exp(R_f_h_apos(2,:)/10)+1);
R_f_h_apos(3,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_3));
R_f_h(3,:)=10*log(exp(R_f_h_apos(3,:)/10)+1);
R_f_h_apos(4,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_4));
R_f_h(4,:)=10*log(exp(R_f_h_apos(4,:)/10)+1);
R_f_h_apos(5,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_5));
R_f_h(5,:)=10*log(exp(R_f_h_apos(5,:)/10)+1);
R_f_h_apos(6,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_6));
R_f_h(6,:)=10*log(exp(R_f_h_apos(6,:)/10)+1);
R_f_h_apos(7,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_7));
R_f_h(7,:)=10*log(exp(R_f_h_apos(7,:)/10)+1);
R_f_h_apos(8,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_8));
R_f_h(8,:)=10*log(exp(R_f_h_apos(8,:)/10)+1);
R_f_h_apos(9,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_9));
R_f_h(9,:)=10*log(exp(R_f_h_apos(9,:)/10)+1);
R_f_h_apos(10,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_10));
R_f_h(10,:)=10*log(exp(R_f_h_apos(10,:)/10)+1);
%are the lines below here correct? ... should huebin_ID_1 be huebin_ID_11
R_f_h_apos(11,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_11));
R_f_h(11,:)=10*log(exp(R_f_h_apos(11,:)/10)+1);
R_f_h_apos(12,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_12));
R_f_h(12,:)=10*log(exp(R_f_h_apos(12,:)/10)+1);
R_f_h_apos(13,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_13));
R_f_h(13,:)=10*log(exp(R_f_h_apos(13,:)/10)+1);
R_f_h_apos(14,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_14));
R_f_h(14,:)=10*log(exp(R_f_h_apos(14,:)/10)+1);
R_f_h_apos(15,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_15));
R_f_h(15,:)=10*log(exp(R_f_h_apos(15,:)/10)+1);
R_f_h_apos(16,:)=100-6.73*mean(deltaE_CAM02(huebin_ID_16));
R_f_h(16,:)=10*log(exp(R_f_h_apos(16,:)/10)+1);

TM30_struct.R_f_h=R_f_h;
TM30_struct.R_f_h_apos=R_f_h_apos;
TM30_struct.R_hs_h=R_hs_h;
TM30_struct.ab_test=ab_test;
TM30_struct.ab_ref=ab_ref;
TM30_struct.Japos_test=Japos_test';
TM30_struct.Japos_ref=Japos_ref';

% TM_30_CVG;

end



function SPDout=cropSPD(SPDin,startWavelength,endWavelength)
% SPDout=cropSPD(SPDin,startWavelength,endWavelength)
%
% Crops a SPD to within a given range or single value,
% Note: Does not interpolate
% ===========================
% Input                     :
% ===========================
% input1, short description data type, size, other info
% ei.
% SPDin          : Spectral power distribution with wavelength in first column and power in secound column
% startWavelength: lower wavelength of crop interval or single wavelength
%                   (giving only startWavelength will give the nearest single value SPD)
% endWavelength  : higher wavelength for end of interval
% ===========================
% Output                    :
% ===========================
% SPDout        :The cropped SPD
%
%
% ===========================
% Code status               :
% ===========================
% Stage		status
% Developing	[ongoing]
% Test		     []
% Release		[-]
%

% ===========================
% Authors
% ===========================
% Anders Thorseth, DTU Fotonik, 2014
%

if nargin<3

    [~,idx]=min(abs(SPDin(:,1)-startWavelength));
    SPDout(1,1)=SPDin(idx,1);
    SPDout(1,2)=SPDin(idx,2);
else

    if startWavelength>endWavelength
        error('start wavelength cannot be larger that endWavelength')
    end

    lambda = SPDin(:,1);
    power  = SPDin(:,2);

    cropPattern = and((lambda>=startWavelength),(lambda<=endWavelength));

    SPDout(:,1) = lambda(cropPattern);
    SPDout(:,2) = power(cropPattern);

end
end

function SPDout=interpSPD1onSPD2(SPD1,SPD2,interpMethod)
%SPD1out=interpSPD1onSPD2(SPD1,SPD2,interpMethod)
%
%Anders Thorseth, DTU Fotonik 2015
%

if nargin<3
    interpMethod='linear';
end

SPDout(:,2) = interp1(SPD1(:,1),SPD1(:,2),SPD2(:,1),interpMethod,0);
SPDout(:,1) = SPD2(:,1);

end

