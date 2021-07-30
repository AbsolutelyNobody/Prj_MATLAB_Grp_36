function A = get_A_matrix(theta5,r5,theta3,r3,theta2,r2)


A = [ 1,1,0,0,0,0,0,0,0,0,0,0,0;
      0,0,1,0,0,0,0,0,0,0,0,0,0;
      0,-1,0,1,0,0,0,0,0,0,0,0,0;
      0,0,-1,0,1,0,0,0,0,0,0,0,0;
      0,-sind(theta5)*r5/2,cosd(theta5)*r5/2,-sind(theta5)*r5/2,cosd(theta5)*r5/2,0,0,0,0,0,0,0,0;
      0,0,0,-1,0,cosd(theta3),1,0,0,0,0,0,0;
      0,0,0,0,-1,sind(theta3),0,1,0,0,0,0,0;
      0,0,0,sind(theta3)*r3/2,-cosd(theta3)*r3/2,0,sind(theta3)*r3/2,-cosd(theta3)*r3/2,0,0,0,0,0;
      0,0,0,0,0,-cosd(theta3),0,0,1,0,0,0,0;
      0,0,0,0,0,-sind(theta3),0,0,0,1,0,0,0;
      0,0,0,0,0,0,0,0,-1,0,1,0,0;
      0,0,0,0,0,0,0,0,0,-1,0,1,0;
      0,0,0,0,0,0,0,0,sind(theta2)*r2/2,-cosd(theta2)*r2/2,sind(theta2)*r2/2,-cosd(theta2)*r2/2,1;
    ];

end
