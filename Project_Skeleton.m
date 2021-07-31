%% Project Skeleton Code %%

clear; clc; close all;

%% initial parameters

%link lengths, units of metres
r1 = 7.8*10^(-2); % o2o3
r2 = 2.5*10^(-2); % o2a2
r3 = 13.8*10^(-2); % o3B
r5 = 4.75*10^(-2); % BC
r7 = 17.1*10^(-2); % o4o3

syms t theta2(t) theta3(t) r4(t) theta5(t) r6(t); %establishes variables are all functions of t

% link 2 movement
dtheta2 = 1800; % deg/s from 300 rpm
ddtheta2 = 0; % constant angular velocity
theta2(t) = t*dtheta2; % theta2 value 

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, using LCEs
% Found using link closure equations as detailed in the report

theta3(t) = atand((r2*sind(theta2(t)))/(r2*cosd(theta2(t))-r1))+180;

r4(t)=(r2*cosd(theta2(t))-r1)/cosd(theta3(t));

theta5(t) = acosd((r7+r3*cosd(theta3(t))) / r5);

r6(t) = - r5*sind(theta5(t)) + r3*sind(theta3(t));

% Hint: Check if the angle needs to be adjusted to its true value
% Hint: Check this for all other angles too

%% Derivative equations of kinematic vars (d/dt) 
syms dtheta3(t) ddtheta3(t) dr4(t) dtheta5(t) ddtheta5(t) dr6(t) ddr6(t)

dtheta3(t) = diff(theta3(t), t); % differentiates theta3_eqn with respect to t
ddtheta3(t) = diff(dtheta3(t), t); % differentiates dtheta3_eqn with respect to t, similar processes for the rest of the variables

dr4(t) = diff(r4(t), t);
ddr4(t) = diff(dr4(t), t);

dtheta5(t) = diff(theta5(t), t);
ddtheta5(t) = diff(dtheta5(t), t);

dr6(t) = diff(r6(t), t);
ddr6(t) = diff(dr6, t);

%% Calculate Values %%
tscale = 100; % scale factor for time (100 means hundredths, etc.
end_time = 60; % number of time increments to calculate

% preallocation of arrays for the first-order values
theta3_array = zeros(end_time, 1);
theta5_array = zeros(end_time, 1);

r4_array = zeros(end_time, 1);
r6_array = zeros(end_time, 1);

% preallocation of arrays for first-derivative values
dtheta3_array = zeros(end_time, 1);
dtheta5_array = zeros(end_time, 1);

dr4_array = zeros(end_time, 1);
dr6_array = zeros(end_time, 1);

% preallocation of arrays for the second-derivative values
t_array = zeros(end_time, 1);
ddtheta3_array = zeros(end_time, 1);
ddtheta5_array = zeros(end_time, 1);

ddr4_array = zeros(end_time, 1);
ddr6_array = zeros(end_time, 1);

for s = 1:1:end_time % calculates values for each variable from of t=0.01 to end_time/tscale seconds
    t_array(s) = (s/tscale);
    
%first order values calculated first
    theta3_array(s) = subs(theta3(t), t, s/tscale);
    theta5_array(s) = subs(theta5(t), t, s/tscale);
    
    r4_array(s) = subs(r4(t), t, s/tscale);
    r6_array(s) = subs(r6(t), t, s/tscale);
    
% first-order derivative terms
    dtheta3_array(s) = subs(dtheta3(t), t, s/tscale);
    dtheta5_array(s) = subs(dtheta5(t), t, s/tscale);
    
    dr4_array(s) = subs(dr4(t), t, s/tscale);
    dr6_array(s) = subs(dr6(t), t, s/tscale);
    
% second-order derivative term calculations
    ddtheta3_array(s) = subs(ddtheta3(t), t, s/tscale);
    ddtheta5_array(s) = subs(ddtheta5(t), t, s/tscale);
    
    ddr4_array(s) = subs(ddr4(t), t, s/tscale);
    ddr6_array(s) = subs(ddr6(t), t, s/tscale);
end

%% Plot vars;

% Plot all desired deliverables. 
% Plotting first-order values first again

% theta3, theta5 plotted against time
figure (1);
plot(t_array, theta3_array, t_array, theta5_array);
grid on;
title('$\theta_3$, $\theta_5$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('\theta_3, \theta_5 (degrees)');
legend('\theta_3', '\theta_5');

% r4 and r6 plotted against time
figure (2);
plot(t_array, r4_array, t_array, r6_array);
grid on;
title('$R_4$, $R_6$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('R_4, R_6 (m)');
legend('R_4', 'R_6');

% Plotting first-derivative values
figure (3);
plot(t_array, dtheta3_array, t_array, dtheta5_array);
grid on;
title('d$\theta_3$, d$\theta_5$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d\theta_3, d\theta_5 (degrees/s)');
legend('d\theta_3', 'd\theta_5');

% dr4 and dr6 plotted against time
figure (4);
plot(t_array, dr4_array, t_array, dr6_array);
grid on;
title('$dR_4$, d$R_6$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('dR_4, dR_6 (m/s)');
legend('dR_4', 'dR_6');

% Plotting second-derivative values

figure (5);
plot(t_array, ddtheta3_array, t_array, ddtheta5_array);
grid on;
title('$d^2\theta_3$, $d^2\theta_5$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d^2\theta_3, d^2\theta_5 (degrees/s^2)');
legend('d^2\theta_3', 'd^2\theta_5');

figure (6);
plot(t_array, ddr4_array, t_array, ddr6_array);
grid on;
title('$d^2R_4$, $d^2R_6$ vs $time$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d^2R_4, d^2R_6 (m/s^2)');
legend('d^2R_4', 'd^2R_6');

% *****************************************************
%% Part 2 - Force and Moment Calculation

syms theta2 theta3(theta2) theta5(theta2) ddtheta3(theta2) ddtheta5(theta2) % redefining variables used in part 1 to be functions of theta2 instead of t

theta3(theta2) = atand((r2*sind(theta2))/(r2*cosd(theta2)-r1))+180;
theta5(theta2) = acosd((r7+r3*cosd(theta3(theta2))) / r5);

ddtheta3(theta2) = diff(theta3(theta2), theta2, 2);
ddtheta5(theta2) = diff(theta5(theta2), theta2, 2);

%%initial parameters:

dtheta2 = 1800; % theta2 dot
ddtheta2 = 0; % theta2 doble-dot - second derivative

rho = 2.7/1000; % density, kg/cm3
r = 0.25; % radius, cm

m2 = pi*(r^2)*rho*r2*100 ; % link 2, O2A2, kg
m3 = pi*(r^2)*rho*r3*100 ; % link 3, O3B, kg
m5 = pi*(r^2)*rho*r5*100 ; % link 5, BC, kg
m4 = 5/1000 ; % slider 4, kg
m6 = 5/1000 ; % slider 6, kg

%Formula: I = 1/12*m*L^2;
IG5 = m5*(r5^2)/12; %Moment of Inertia, link 5, kg*m^2
IG3 = m3*(r3^2)/12; %MOI, link 3, kg*m^2
IG2 = m2*(r2^2)/12; %MOI, link 2, kg*m^2

M12_list = [];
theta2_list = [];
Fs_list = [];  % shaking force
Fsx_list = [];
Fsy_list = [];
As_list = []; % direction of a shaking force
Ms_list =[]; % Shaking moment

%List of Forces and the Angles at which forces act%
F16x_list = [];
F16_list = [];
F16_alpha = [];

F56x_list = [];
F56y_list = [];
F56_list = [];
F56_alpha = [];

F35x_list = [];
F35y_list = [];
F35_list = [];
F35_alpha = [];

F34_list = [];
F34_alpha = [];

F13x_list = [];
F13y_list = [];
F13_list = [];
F13_alpha = [];

F24x_list = [];
F24y_list = [];
F24_list = [];
F24_alpha = [];

F12x_list = [];
F12y_list = [];
F12_list = [];
F12_alpha = [];

for theta2 = 0:1:360

    % kinematic variables are caculated based on loop eqn
    syms RG2(theta2) V2(theta2) A2(theta2) RG3(theta2) V3(theta2) A3(theta2)
    syms RG4(theta2) V4(theta2) A4(theta2) RG5(theta2) V5(theta2) A5(theta2) RG6(theta2) V6(theta2) A6(theta2)
    % kinematic variables are caculated based on loop eqn
    RG2(theta2) = [r2/2*cosd(theta2) , r2/2*sind(theta2)];
    V2(theta2) = diff(RG2(theta2),theta2);
    A2(theta2) = diff(RG2(theta2),theta2,2);
    
    
    RG3(theta2) = [r1 + r3/2*cosd(theta3(theta2)) , r3/2*sind(theta3(theta2))];
    V3(theta2) = diff(RG3(theta2),theta2);
    A3(theta2) = diff(RG3(theta2),theta2,2);
    
    RG4(theta2) = [r2*cosd(theta2) , r2*sind(theta2)];
    V4(theta2) = diff(RG4(theta2),theta2);
    A4(theta2) = diff(RG4(theta2),theta2,2);
    
    RG5(theta2) = [r1-r7+r5/2*cosd(theta5(theta2)) , r6+r5/2*sind(theta5(theta2))];
    V5(theta2) = diff(RG5(theta2),theta2);
    A5(theta2) = diff(RG5(theta2),theta2,2);
    
    RG6(theta2) = [r1-r7 , r6];
    V6(theta2) = diff(RG6(theta2),theta2);
    A6(theta2) = diff(RG6(theta2),theta2,2);

    ag2x = A2(theta2)*cosd(theta2);
    ag2y = A2(theta2)*sind(theta2);
    ag3x = A3(theta2)*cosd(theta2);
    ag3y = A3(theta2)*sind(theta2);
    ag4x = A4(theta2)*cosd(theta2);
    ag4y = A4(theta2)*sind(theta2);
    ag5x = A5(theta2)*cosd(theta2);
    ag5y = A5(theta2)*sind(theta2);
    ag6x = A6(theta2)*cosd(theta2);
    ag6y = A6(theta2)*sind(theta2);
    
    B = get_ma_vector(m2,m3,m4,m5,m6,ag2x,ag2y,ag3x,ag3y,ag4x,ag4y,ag5x,ag5y,ag6x,ag6y,IG3,IG5,ddtheta3(theta2),ddtheta5(theta2));
    
    A = get_A_matrix(theta5(theta2),r5,theta3(theta2),r3,theta2,r2);

    x = A\ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
  % Collecting force magnitudes from x vector and adding to respective lists:
    
    F16x = x(1);
    F16x_list = [F16x_list; F16x];
    F56x = x(2);
    F56x_list = [F56x_list; F56x];
    F56y = x(3);
    F56y_list = [F56y_list; F56y];
    F35x = x(4);
    F35x_list = [F35x_list; F35x];
    F35y = x(5);
    F35y_list = [F35y_list; F35y];
    F34 = x(6);
    F34_list = [F34_list; F34];
    F13x = x(7);
    F13x_list = [F13x_list; F13x];
    F13y = x(8);
    F13y_list = [F13y_list; F13y];
    F24x = x(9);
    F24x_list = [F24x_list; F24x];
    F24y = x(10);
    F24y_list = [F24y_list; F24y];
    F12x = x(11);
    F12x_list = [F12x_list; F12x];
    F12y = x(12);
    F12y_list = [F12y_list; F12y];
    M12 = x(13);
    M12_list = [M12_list; M12];
    
    %Shaking Force Equations from kinetic analysis
    Fsx = -F16x - F12x - F13x
    Fsx_list = [Fsx_list; Fsx];
    Fsy = -F12y - F13y
    Fsy_list = [Fsy_list; Fsy];
    
    % Magnitudes of all forces: 
    % Atan is defined on [-pi/2; pi/2]. 
    % This if clause will help to adjust the value of the angle 
    % to its true value:	

    
    % Directions of all forces:  
    fx = Fsx;
    fy = Fsy;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    Fs_alpha = [Fs_alpha; alpha_f];
    
    fx = F16x;
    fy = 0;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F16_alpha = [F16_alpha; alpha_f];

    fx = F56x;
    fy = F56y;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F56_alpha = [F56_alpha; alpha_f];
    
    fx = F35x;
    fy = F35y;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F35_alpha = [F35_alpha; alpha_f];
    
    F34_alpha = [F34_alpha; theta3(theta2)];
    
    fx = F13x;
    fy = F13y;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F13_alpha = [F13_alpha; alpha_f];
    
    fx = F24x;
    fy = F24y;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F24_alpha = [F24_alpha; alpha_f];

    fx = F12x;
    fy = F12y;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    F12_alpha = [F12_alpha; alpha_f];

  
    % Collecting the values of theta2:
    theta2_list = [theta2_list, theta2];
    
    %Finding force length from the x and y components:
    F16_list = F16x_list;
    F56_list = sqrt((F56x_list.^2)+(F56y_list.^2));
    F35_list = sqrt((F35x_list.^2)+(F35y_list.^2));
    F13_list = sqrt((F13x_list.^2)+(F13y_list.^2));
    F24_list = sqrt((F24x_list.^2)+(F24y_list.^2));
    F12_list = sqrt((F12x_list.^2)+(F12y_list.^2));
    Fs_list = sqrt((Fsx_list.^2)+(Fsy_list.^2));
    
    
end


% Regular and Polar plots:
% Might have to transpose the Force vectors for polar plot. Do so if needed
% Polar plot only works with radians so will have to do it accordingly

figure (9);
plot(theta2_list,M12_list);
grid on;
title('M_{12} vs \theta_2');
xlabel('\theta_2   unit: degree');
ylabel('M12   unit: N-m');


% Convert degrees to the radians
theta2_rad = deg2rad(theta2_list);

figure (10);
polarplot(Fs_alpha,Fs_list);
grid on;
title('F_s polar plot');

figure (11);
polarplot(F16_alpha,F16_list);
grid on;
title('F_{16} polar plot');

figure (12);
polarplot(F56_alpha,F56_list);
grid on;
title('F_{56} polar plot');

figure (13);
polarplot(F35_alpha,F35_list);
grid on;
title('F_{35} polar plot');

figure (14);
polarplot(F34_alpha,F34_list);
grid on;
title('F_{34} polar plot');

figure (15);
polarplot(F13_alpha,F13_list);
grid on;
title('F_{13} polar plot');

figure (16);
polarplot(F24_alpha,F24_list);
grid on;
title('F_{24} polar plot');

figure (17);
polarplot(F12_alpha,F12_list);
grid on;
title('F_{12} polar plot');



