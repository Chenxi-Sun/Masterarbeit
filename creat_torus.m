function [x,y,z]=creat_torus(R,r0,u,v)
u=deg2rad(u);
v=deg2rad(v);
x=(R+r0.*cos(v)).*cos(u);
y=r0*sin(v)
z=(R+r0.*cos(v)).*sin(u);