function ma = get_ma_vector(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);

% ENTER YOUR CODE HERE %
m6*ag5x=F35x-F65x
m6*ag5y=F35y-F65y
IG5*a5=-(F35x+F65x)*sind(theta5)*r5/2 + (F35y+F65y)*cosd(theta5)*r5/2

m3*ag3x=F13x+F43*cosd(theta3)-F53x
m3*ag3y=F13y+F43*sind(theta3)-F53y
IG3*a3=(F13x+F53x)*sind(theta3)*r3/2 - (F13y+F53y)*cosd(theta3)*r3/2

m4*ag4x=F24x-F34*cosd(theta3)
m4*ag4y=F24y-F34*sind(theta3)

m2*ag2x=F12x-F42x
m2*ag2y=F12y-F42y
IG2*a2= 0 =(F12x+F42x)*sind(theta2)*r2/2 - (F12y+F42y)*cosd(theta2)*r2/2 + M12

%shaking force and moment
Fsx=-F61x-F21x-F31x
Fsy=-F21y-F31y
Ms=-M21-r1*F31y+r6*F61x %only for r6 upwards. don't know what to do about when it is belox the x axis, maybe use vector instead?

ma = [ % mi * a_g_ix; ...
       % mi * a_g_iy; ... etc
       % ENTER YOUR CODE HERE %
    ];

end
