function [T_x_Parabolic,T_x_Triangulation_corr,v_Tx,paraPoly]=Planck_Table(u_t,v_t)


% Generate Planckian (blackbody) source 
% look-up table for CCT and Duv

% Load table from IES-TM-30 online calc. 

persistent Table_Planck

 if isempty(Table_Planck)
    % load data
    filenamePlanck = fullfile('source_PlanckTable.csv');
    Table_Planck = csvread(filenamePlanck);

 end



for ttt=1:size(Table_Planck,1)
 % distance between planck u,v and test u,v
    Table_Planck(ttt,4)=sqrt((u_t-Table_Planck(ttt,2))^2+(v_t-Table_Planck(ttt,3))^2);
    Table_Planck(ttt,5)=ttt;
    
end 
        
%find minumum value in Planck table
[minimum,m_m]=min(Table_Planck);
m_m_location=m_m(1,4);
      
%if chromaticity is outside of boundaries choose the min or max value
if m_m_location>m_m(1,2)
    m_m_location=m_m(1,2)-1;
elseif m_m_location<m_m(1,1)+1
    m_m_location=m_m(1,1)+1;
else
end

% find T_x by triangulation 
l_m=((Table_Planck(m_m_location+1,2)-Table_Planck(m_m_location-1,2))^2+(Table_Planck(m_m_location+1,3)-Table_Planck(m_m_location-1,3))^2)^(1/2);
x_m=(Table_Planck(m_m_location-1,4)^2-Table_Planck(m_m_location+1,4)^2+l_m^2)/(2*l_m);
T_x_Triangulation=Table_Planck(m_m_location-1,1)+(Table_Planck(m_m_location+1,1)-Table_Planck(m_m_location-1,1))*x_m/l_m;

% find Duv
v_Tx=Table_Planck(m_m_location-1,3)+(Table_Planck(m_m_location+1,3)-Table_Planck(m_m_location-1,3))*x_m/l_m;

if v_t-v_Tx >=0 
Duv_test=sqrt(Table_Planck(m_m_location-1,4)^2-x_m^2);  
else
Duv_test=sqrt(Table_Planck(m_m_location-1,4)^2-x_m^2)*-1;    
end

% CCT parabolic method 
X_poly=(Table_Planck(m_m_location+1,1)-Table_Planck(m_m_location,1))*(Table_Planck(m_m_location-1,1)-Table_Planck(m_m_location+1,1))*(Table_Planck(m_m_location,1)-Table_Planck(m_m_location-1,1));
a_poly=(Table_Planck(m_m_location-1,1)*(Table_Planck(m_m_location+1,4)-Table_Planck(m_m_location,4))+Table_Planck(m_m_location,1)*(Table_Planck(m_m_location-1,4)-Table_Planck(m_m_location+1,4))+Table_Planck(m_m_location+1,1)*(Table_Planck(m_m_location,4)-Table_Planck(m_m_location-1,4)))*X_poly^-1;
b_poly=-(Table_Planck(m_m_location-1,1)^2*(Table_Planck(m_m_location+1,4)-Table_Planck(m_m_location,4))+Table_Planck(m_m_location,1)^2*(Table_Planck(m_m_location-1,4)-Table_Planck(m_m_location+1,4))+Table_Planck(m_m_location+1,1)^2*(Table_Planck(m_m_location,4)-Table_Planck(m_m_location-1,4)))*X_poly^-1; 
c_poly=-((Table_Planck(m_m_location-1,4)*(Table_Planck(m_m_location+1,1)-Table_Planck(m_m_location,1))*Table_Planck(m_m_location,1)*Table_Planck(m_m_location+1,1))+(Table_Planck(m_m_location,4)*(Table_Planck(m_m_location-1,1)-Table_Planck(m_m_location+1,1))*Table_Planck(m_m_location-1,1)*Table_Planck(m_m_location+1,1))+(Table_Planck(m_m_location+1,4)*(Table_Planck(m_m_location,1)-Table_Planck(m_m_location-1,1))*Table_Planck(m_m_location-1,1)*Table_Planck(m_m_location,1)))*X_poly^-1; 
%T_x closest point on planckian radiation
T_x_Parabolic=-b_poly/(2*a_poly);

% Linear shift according to Ohno CCT calculation in TM-30 v2.04
T_x_Triangulation_corr= T_x_Triangulation + (T_x_Parabolic - T_x_Triangulation)* (sqrt(Table_Planck(m_m_location-1,4)^2-x_m^2))*(1/0.002);

paraPoly.a=a_poly;
paraPoly.b=b_poly;
paraPoly.c=c_poly;
