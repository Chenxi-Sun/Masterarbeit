t=0:pi/50:6*pi;   % to have one complete round
r = 4;            % radius
h = 8;            % height
x = r * sin(t);
y = r * cos(t);
z = h/(2*pi) * t;   
plot3(x,y,z,'LineWidth',2)

hold on 
gm = multicylinder(4,24)
model = createpde;
model.Geometry = gm
pdegplot(model,"CellLabels","on")
% r = 2; % radius of cylinder
% h = 2; % height of cylinder
% circumference_pnts = 100; % number of points in the circumference
% theta = linspace(0, 2*pi, circumference_pnts); % angle to compute 'x' and 'y'
% x = repmat(r*cos(theta),2,1); % compute coordinates and put in appropriate form
% y = repmat(r*sin(theta),2,1); % compute coordinates and put in appropriate form
% z = [zeros(1,circumference_pnts);h*ones(1,circumference_pnts)]; % array of 'z' values: first row is basis of cylinder and second row is top of cylinder
% surf(x,y,z)