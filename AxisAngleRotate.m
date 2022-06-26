function [VnewX,VnewY,VnewZ]=AxisAngleRotate(X,Y,Z,direction,alpha)
%euler angles transformation

%rotationsangle
theta=deg2rad(alpha);

%rotate around x-axis
if direction==X
    VnewX=X;
    VnewY=Y*cos(theta)-Z*sin(theta);
    VnewZ=Y*sin(theta)+Z*cos(theta);
end 
%rotate around y-axis
if direction==Y
    VnewX=X*cos(theta)+Z*sin(theta);
    VnewY=Y;
    VnewZ=Z*cos(theta)-X*sin(theta);
end 
%rotate around z-axis
if direction==Z
    VnewX=X*cos(theta)-Y*sin(theta);
    VnewY=X*sin(theta)+Y*cos(theta);
    VnewZ=Z;
end 

end 