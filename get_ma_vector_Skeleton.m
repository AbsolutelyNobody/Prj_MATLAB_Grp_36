function ma = get_ma_vector(m2,m3,m4,m5,m6,ag2x,ag2y,ag3x,ag3y,ag4x,ag4y,ag5x,ag5y,ag6x,ag6y,IG3,IG5,ddtheta3,ddtheta5);

%Below are the force equations used to get the quantities listed on the
%left side (that are used in the overall ma vector) for reference

% m6*ag6x=F16x+F56x
% m6*ag6y=F56y
% 
% m5*ag5x=F35x-F65x
% m5*ag5y=F35y-F65y
% IG5*ddtheta5(t)=-(F35x+F65x)*sind(theta5)*r5/2 + (F35y+F65y)*cosd(theta5)*r5/2
% 
% m3*ag3x=F13x+F43*cosd(theta3)-F53x
% m3*ag3y=F13y+F43*sind(theta3)-F53y
% IG3*ddtheta3(t)=(F13x+F53x)*sind(theta3)*r3/2 - (F13y+F53y)*cosd(theta3)*r3/2
% 
% m4*ag4x=F24x-F34*cosd(theta3)
% m4*ag4y=F24y-F34*sind(theta3)
% 
% m2*ag2x=F12x-F42x
% m2*ag2y=F12y-F42y
% IG2*ddtheta2= 0 =(F12x+F42x)*sind(theta2)*r2/2 - (F12y+F42y)*cosd(theta2)*r2/2 + M12

a = m6*ag6x;
b = m6*ag6y;
c = m5*ag5x;
d = m5*ag5y;
e = IG5*ddtheta5;
f = m3*ag3x;
g = m3*ag3y;
h = IG3*ddtheta3;
i = m4*ag4x;
j = m4*ag4y;
k = m2*ag2x;
l = m2*ag2y;
m = 0;

ma = [a, b, c, d, e, f, g, h, i, j, k, l, m];
end
