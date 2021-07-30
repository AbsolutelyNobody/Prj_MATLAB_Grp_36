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

for s = 1:1:end_time % calculates values over the course of t=0.01s to 0.90s
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
title('$time$ vs $\theta_3$, $\theta_5$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('\theta_3, \theta_5 (degrees)');
legend('\theta_3', '\theta_5');

% r4 and r6 plotted against time
figure (2);
plot(t_array, r4_array, t_array, r6_array);
grid on;
title('$time$ vs $R_4$, $R_6$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('R_4, R_6 (m)');
legend('R_4', 'R_6');

% Plotting first-derivative values
figure (3);
plot(t_array, dtheta3_array, t_array, dtheta5_array);
grid on;
title('$time$ vs d$\theta_3$, d$\theta_5$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d\theta_3, d\theta_5 (degrees/s)');
legend('d\theta_3', 'd\theta_5');

% dr4 and dr6 plotted against time
figure (4);
plot(t_array, dr4_array, t_array, dr6_array);
grid on;
title('$time$ vs d$R_4$, d$R_6$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('dR_4, dR_6 (m/s)');
legend('dR_4', 'dR_6');

% Plotting second-derivative values

figure (5);
plot(t_array, ddtheta3_array, t_array, ddtheta5_array);
grid on;
title('$time$ vs $d^2\theta_3$, $d^2\theta_5$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d^2\theta_3, d^2\theta_5 (degrees/s^2)');
legend('d^2\theta_3', 'd^2\theta_5');

figure (6);
plot(t_array, ddr4_array, t_array, ddr6_array);
grid on;
title('$time$ vs $d^2R_4$, $d^2R_6$', 'Interpreter','latex');
xlabel('time (s)');
ylabel('d^2R_4, d^2R_6 (m/s^2)');
legend('d^2R_4', 'd^2R_6');
