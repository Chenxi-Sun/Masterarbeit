clc
clear all
close all

wlarmor = larmorfrq('1H',1200);
sigma   = 0.1;      %linienbreite der Gausskurven
frq     = (wlarmor-10:0.01:wlarmor+10);
tau     = 0.3;      %in us

aiso    = 3;
T       = 1; 

theta = (0:90)*pi/180;				%generiert alle Winkel von 0 bis pi/2
A     = aiso+T*(3*((cos(theta)).^2)-1);		%berechnet Parameter A
B     = 3*T*sin(theta).*cos(theta);			%berechnet Parameter B

%Berechnung der Resonanzfrequenzen
wa = ((wlarmor+(A/2)).^2+B.^2/4).^0.5;
wb = ((wlarmor-(A/2)).^2+B.^2/4).^0.5;

%Intensitaetsverteilung fuer die erlaubten und verbotenen Uebergaenge
Iallowed   = abs(wlarmor^2-0.25*(wa-wb).^2)/wa.*wb;				
Iforbidden = abs(wlarmor^2-0.25*(wa+wb).^2)/wa.*wb;				


%ordnet jedem Wert wa und wb eine Gausskurven zu
lauf=1;
for lauf=1:length(wa)
gaussa(:,lauf)=1/(2*pi*sigma^2)^0.5*exp(-(frq-wa(lauf)).^2/(2*sigma^2));		
gaussb(:,lauf)=1/(2*pi*sigma^2)^0.5*exp(-(frq-wb(lauf)).^2/(2*sigma^2));		
end

%Gewichtung mit sin(theta) und Multiplikation mit der Intensitaetsverteilung fuer verbotene Uebergaenge
Y1 = sum(sin(theta).*Iforbidden.*gaussa,2);				
Y2 = sum(sin(theta).*Iforbidden.*gaussb,2);

%Gewichtung mit sin(theta) und Multiplikation mit der Intensitaetsverteilung fuer erlaubten Uebergaenge
Y3 = sum(sin(theta).*Iallowed.*gaussa,2);
Y4 = sum(sin(theta).*Iallowed.*gaussb,2);

%gesamtes Spektrum der verbotenen Uebergaenge
Yforb      = Y1+Y2;

%gesamtes Spektrum der verbotenen Uebergaenge mit blindspots
Yforbblind = (Y1+Y2).*sin((frq.'-wlarmor)*2*pi*tau).^2;

%gesamtes Spektrum der erlaubten Uebergaenge
Yall       = (Y3+Y4);

%gesamtes Spektrum der erlaubten Uebergaenge mit blindspots
Yallblind  = (Y3+Y4).*sin((frq.'-wlarmor)*2*pi*tau).^2;


figure
plot(frq-wlarmor,Yall/max(Yall))                  %beide normalen

%xlim([-6 6])
ylim([0 1.1])
xlabel('\nu / MHz')
% T=-(Aortho-Aparallel)/3;
% aiso=Aparallel+2*(Aortho-Aparallel)/3;