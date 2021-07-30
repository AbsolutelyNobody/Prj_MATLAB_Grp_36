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

r6(t) = r5*sind(theta5(t)) - r3*sind(theta3(t));

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
end_time = 20; % number of time increments to calculate

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
t_array = zeros(end_time, 1)
ddtheta3_array = zeros(end_time, 1);
ddtheta5_array = zeros(end_time, 1);

ddr4_array = zeros(end_time, 1);
ddr6_array = zeros(end_time, 1);

for s = 1:1:end_time % calculates values over the course of t=0.01s to 0.90s
    t_array(s) = (s/tscale)
    
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
% 
% figure (1)
% plot(t, theta3_array)
% grid on;
% title('$\theta_3$ vs $\theta_2$', 'Interpreter','latex')
% xlabel('\theta_2   unit: degree')
% ylabel('\theta_3   unit: degree')
