function [neu]=drehen(xaxis,yaxis,zaxis,alf,bett,gammel)

% xtest=[1 0 0]'; ytest=[0 1 0]'; ztest=[0 0 1]';

NNnew = [xaxis(1,:);yaxis(1,:);zaxis(1,:)]';

% alf = 0; bett = 180; gammel = 0; 
eulerdreh = [cosd(alf)*cosd(gammel)-sind(alf)*cosd(bett)*sind(gammel) sind(alf)*cosd(gammel)+cosd(alf)*cosd(bett)*sind(gammel) sind(bett)*sind(gammel); 
    -cosd(alf)*sind(gammel)-sind(alf)*cosd(bett)*cosd(gammel) -sind(alf)*sind(gammel)+cosd(alf)*cosd(bett)*cosd(gammel) sind(bett)*cosd(gammel);    
    sind(alf)*sind(bett) -cosd(alf)*sind(bett) cosd(bett)];

neu = (NNnew*eulerdreh'*inv(NNnew))*NNnew;

%%
% figure
% quiver3(M(1,1),M(1,2),M(1,3), xaxis(1,1), xaxis(1,2), xaxis(1,3),'r','Linewidth',2)
% hold on
% quiver3(M(1,1),M(1,2),M(1,3), yaxis(1,1), yaxis(1,2), yaxis(1,3),'g','Linewidth',2)
% hold on 
% quiver3(M(1,1),M(1,2),M(1,3), zaxis(1,1), zaxis(1,2), zaxis(1,3),'b','Linewidth',2)
% hold on
% axis equal
% 
% figure
% quiver3(M(1,1),M(1,2),M(1,3),neu(1,1), neu(2,1), neu(3,1),'r','Linewidth',2)
% hold on
% quiver3(M(1,1),M(1,2),M(1,3),neu(1,2), neu(2,2), neu(3,2),'g','Linewidth',2)
% hold on 
% quiver3(M(1,1),M(1,2),M(1,3),neu(1,3), neu(2,3), neu(3,3),'b','Linewidth',2)
% 
% axis equal



% 
% alf=radtodeg(o2(1,1));
% bett=radtodeg(o2(1,2));
% gammel=radtodeg(o2(1,3));

% alf=pi/2;
% bett=0;
% gammel=0;

%[alf,bett,gammel] = EulerAngles_nic(xaxis1(1,:),yaxis1(1,:),zaxis1(1,:),[1 0 0],[0 1 0],[0 0 1]);

% NNold=[1 0 0; 0 1 0; 0 0 1];
% NNnew=[1 0 0; 0 0 1; 0 -1 0];
% % 
% % Nold.x=[1 0 0];
% % Nold.y=[0 1 0];
% % Nold.z=[0 0 1];
% % 
% % Nnew.x=[0 0 -1];
% % Nnew.y=[0 -1 0];
% % Nnew.z=[-1 0 0];
% 
% %dreh=[cos(bett) 0 sin(bett); 0 1 0; -sin(bett) 0 cos(bett)]
% 
% 
% neu=(NNnew'*dreh*inv(NNnew'))
% 
% 
% % [aba2,baba2,gaba2] = EulerAngles_nic(Nnew.x,Nnew.y,Nnew.z,Nold.x,Nold.y,Nold.z);
% % [aba,baba,gaba] = EulerAngles_nic(Nold.x,Nold.y,Nold.z,Nnew.x,Nnew.y,Nnew.z);
% [aba,baba,gaba] = EulerAngles_nic(NNold,NNnew);
% % [alf,bett,gammel] = EulerAngles(Nold,Nnew);
% %alf=radtodeg(alf);
% 
% alf = radtodeg(aba);
% bett = radtodeg(baba);
% gammel = radtodeg(gaba);
% 
% 
% 
% % eulerdreh = [cos(gammel)*cos(alf)-cos(bett)*sin(alf)*sin(gammel) cos(gammel)*sin(alf)+cos(bett)*cos(alf)*sin(gammel) sin(gammel)*sin(bett);
% %     -sin(alf)*cos(alf)-cos(bett)*sin(alf)*cos(gammel) -sin(gammel)*sin(alf)+cos(bett)*cos(alf)*cos(gammel) cos(gammel)*sin(bett);
% %     sin(bett)*sin(alf) -sin(bett)*cos(alf) cos(bett)];
% 
% eulerdreh = [cosd(alf)*cosd(gammel)-sind(alf)*cosd(bett)*sind(gammel) sind(alf)*cosd(gammel)+cosd(alf)*cosd(bett)*sind(gammel) sind(bett)*sind(gammel); 
%     -cosd(alf)*sind(gammel)-sind(alf)*cosd(bett)*cosd(gammel) -sind(alf)*sind(gammel)+cosd(alf)*cosd(bett)*cosd(gammel) sind(bett)*cosd(gammel);    
%     sind(alf)*sind(bett) -cosd(alf)*sind(bett) cosd(bett)];
% % 
% zyzdreh = [-sind(alf)*sind(gammel)+cosd(alf)*cosd(bett)*cosd(gammel) cosd(alf)*sind(gammel)+sind(alf)*cosd(bett)*cosd(gammel) -sind(bett)*cosd(gammel);
%     -sind(alf)*cosd(gammel)-cosd(alf)*cosd(bett)*sind(gammel) cosd(alf)*cosd(gammel)-sind(alf)*cosd(bett)*sind(gammel) sind(bett)*sind(gammel);
%     cosd(alf)*sind(bett) sind(alf)*sind(bett) cosd(bett)];
% 
% % xaxis1(2,:)'
% % (inv(eulerdreh)*dreh*eulerdreh)*vtest'
% M(:,1)=eulerdreh'*xtest';
% M(:,2)=eulerdreh'*ytest';
% M(:,3)=eulerdreh'*ztest'

%(eulerdreh*dreh*inv(eulerdreh))*vtest'%xaxis1(1,:)'
% dreh*eulerdreh*xaxis1(1,:)'


% zyzdreh*vtest'%xaxis1(1,:)'