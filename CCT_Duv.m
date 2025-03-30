function [CCT_test,Duv_test]=CCT_Duv(Stest)

persistent SourcefileTM30

 if isempty(SourcefileTM30)
    % load data
    filenameTM30 = fullfile('source_TM30_20.csv');
    SourcefileTM30 = csvread(filenameTM30);

 end



wavelength=SourcefileTM30(:,100);
%CIE 2 and 10 degree observer
xbar=SourcefileTM30(:,101);
ybar=SourcefileTM30(:,102);
zbar=SourcefileTM30(:,103);


% CCT and Duv calculation from Ohno (2014). Practical use and calculation of CCT and Duv. Leukos, 10(1), 47-55.
% combined method. Triangular (Duv<0.002) and parabolic solution

% Tristimulus values - 2 degree observer 
X_CCT_t=sum(Stest.*xbar);
Y_CCT_t=sum(Stest.*ybar);
Z_CCT_t=sum(Stest.*zbar);

% 1960 u v coordinates 
u_t=(4*X_CCT_t)/(X_CCT_t+15*Y_CCT_t+3*Z_CCT_t);
v_t=(6*Y_CCT_t)/(X_CCT_t+15*Y_CCT_t+3*Z_CCT_t);

L_fp=sqrt((u_t-0.292)^2+(v_t-0.24)^2);
a_CCT=acos((u_t-0.292)/L_fp);
L_bb=(-0.00616793*a_CCT^6)+(0.0893944*a_CCT^5)+(-0.5179722*a_CCT^4)+(1.5317403*a_CCT^3)+(-2.4243787*a_CCT^2)+(1.925865*a_CCT)-0.471106;
Duv_test_initial=L_fp-L_bb;

% Find Duv. if |Duv| <0.002 triangular, otherwise parabolic 

if abs(Duv_test_initial) < 0.002
    % triangular solution 
    
    % New plack table function: 
    [~,T_x_Triangulation_corr]=Planck_Table(u_t,v_t);
    
    % CCT 
    T_x=T_x_Triangulation_corr;
      Duv_test=Duv_test_initial;
else
    % parabolic solution 
    
    % New plack table function: 
    [T_x_Parabolic,~,v_Tx,paraPoly]=Planck_Table(u_t,v_t);

    % CCT 
    T_x=T_x_Parabolic;
    
    % Duv for parabolic
    if v_t-v_Tx >=0 
    Duv_test=(paraPoly.a*(T_x)^2)+(paraPoly.b*T_x)+paraPoly.c;  
    else
    Duv_test=-((paraPoly.a*(T_x)^2)+(paraPoly.b*T_x)+paraPoly.c);    
    end

end

% CCT of the test light source 
CCT_test=T_x;
% 1931 x y coordinates of the test light source 
x_test=(X_CCT_t)/(X_CCT_t+Y_CCT_t+Z_CCT_t);
y_test=(Y_CCT_t)/(X_CCT_t+Y_CCT_t+Z_CCT_t);
% 1976 u' v' coordinates of the test light source 
uprime_test=(4*X_CCT_t)/(X_CCT_t+15*Y_CCT_t+3*Z_CCT_t);
vprime_test=(9*Y_CCT_t)/(X_CCT_t+15*Y_CCT_t+3*Z_CCT_t);
