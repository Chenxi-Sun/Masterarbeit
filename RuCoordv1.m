function [o1,o2,r]=RuCoordv1(R,X1,Y1,Z1,X2,Y2,Z2)

%X,Y,Z are the Coordinate axes of the spin labels and R is the Distance vector. o are the euler angles. 
for k1=1:size(R,1)
r(k1,:)=norm(R(k1,:));
x1(k1,:)=X1(k1,:)./norm(X1(k1,:));
y1(k1,:)=Y1(k1,:)./norm(Y1(k1,:)); 
z1(k1,:)=Z1(k1,:)./norm(Z1(k1,:));
x2(k1,:)=X2(k1,:)./norm(X2(k1,:));
y2(k1,:)=Y2(k1,:)./norm(Y2(k1,:)); 
z2(k1,:)=Z2(k1,:)./norm(Z2(k1,:));


Z0(k1,:)=R(k1,:)./(norm(R(k1,:))); 
a=dot(Z0(k1,:),Z1(k1,:));
 if a<0
     Z0(k1,:)=-Z0(k1,:);
 end

% eps=10.^-9;
% if abs(a)>=1-eps
%     X0=X1; Y0=Y1; Z0=Z1;
% else
    L(k1,:)=cross(Z0(k1,:),Z1(k1,:));
    X0(k1,:)=L(k1,:)./sqrt(dot(L(k1,:),L(k1,:)));  %sqrt(dot(L(k1,:),L(k1,:)))=norm(L(k1,:))
    Y0(k1,:)=cross(Z0(k1,:),X0(k1,:));
%  end

[o1(k1,1),o1(k1,2),o1(k1,3)]=eulers(X0(k1,:),Y0(k1,:),Z0(k1,:),x1(k1,:),y1(k1,:),z1(k1,:));
[o2(k1,1),o2(k1,2),o2(k1,3)]=eulers(X0(k1,:),Y0(k1,:),Z0(k1,:),x2(k1,:),y2(k1,:),z2(k1,:));

end
% function [A,B,C]=eulers(X1,Y1,Z1,X2,Y2,Z2)
% eps=10.^-9;
% b=dot(Z1,Z2);
% B=angle(b+1i.*sqrt(1-b.^2));
% 
% if abs(b)>=1-eps
%     A=0; C=angle(dot(X2,X1)+1i.*dot(X2,Y1));
% else
%     L=cross(Z1,Z2); l=L./sqrt(dot(L,L));
%     A=angle(dot(l,X1)+1i.*dot(l,Y1));
% %     C=-angle(dot(l,X2)+1i.*dot(l,Y2));
%      yy=cross(Z2,l);
%          ax=dot(l,X2); ay=dot(yy,X2);
%      C=angle(ax+1i.*ay);
% end
% % gr=pi/180;
% A=A./gr;
% B=B./gr;
% C=C./gr;
