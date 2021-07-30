function ma = get_ma_vector(m2,m3,m4,m5,m6,ag2x,ag2y,ag3x,ag3y,ag4x,ag4y,ag5x,ag5y,ag6x,ag6y,IG3,IG5,ddtheta3(t),ddtheta5(t));

m6*ag6x=F16x+F56x
m6*ag6y=F56y

m5*ag5x=F35x-F65x
m5*ag5y=F35y-F65y
IG5*ddtheta5(t)=-(F35x+F65x)*sind(theta5)*r5/2 + (F35y+F65y)*cosd(theta5)*r5/2

m3*ag3x=F13x+F43*cosd(theta3)-F53x
m3*ag3y=F13y+F43*sind(theta3)-F53y
IG3*ddtheta3(t)=(F13x+F53x)*sind(theta3)*r3/2 - (F13y+F53y)*cosd(theta3)*r3/2

m4*ag4x=F24x-F34*cosd(theta3)
m4*ag4y=F24y-F34*sind(theta3)

m2*ag2x=F12x-F42x
m2*ag2y=F12y-F42y
IG2*ddtheta2= 0 =(F12x+F42x)*sind(theta2)*r2/2 - (F12y+F42y)*cosd(theta2)*r2/2 + M12

%shaking force and moment
Fsx=-F61x-F21x-F31x
Fsy=-F21y-F31y
Ms=-M21-r1*F31y+r6*F61x %only for r6 upwards. don't know what to do about when it is belox the x axis, maybe use vector instead?

ma = [ m6*ag6x ; 
       m6*ag6y ;
       m5*ag5x ;
       m5*ag5y ;
       IG5*ddtheta5(t) ;
       m3*ag3x ;
       m3*ag3y ;
       IG3*ddtheta3(t) ;
       m4*ag4x ;
       m4*ag4y ;
       m2*ag2x ;
       m2*ag2y ;
       IG2*ddtheta2 ;
    ];

end
