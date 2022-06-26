function [X2 Y2 Z2] = drehnung(X,Y,Z,theta)
Vver=-Z;
Vhori=Y;
V=Vver+Vhori;
k=Vver/norm(Vver);
vrot =cos(theta).*V+((1-cos(theta)).*dot(k,V).*k)+sin(theta).*cross(k,V);
Y2 = vrot-Vver;
Z2 = -Vver;
X2 = cross(Y2,Z2);
end 