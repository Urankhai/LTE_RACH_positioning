function [ x,y ] = geo2gausskrueger( L, B, L0 )
%GEO2GAUSSKRUEGER Summary of this function goes here
% L longitude 
% B latitude (80S-80N)
% L0 longitude reference meridian, should be a multiple of 3deg, 
%    if eq. to -999 automatic zone is selected for the data


B0=floor(B);
if L0==-999
    L0=L-mod(L,3);
    if L-L0>=2
        L0=L0+3;
    end
end

b=(B-B0)*pi/180;
l=(L-L0)*pi/180;

sinbo=sin(B0*pi/180);
cosbo=cos(B0*pi/180);
sin2bo=sinbo^2;
sin4bo=sinbo^4;
sin6bo=sinbo^6;
sin8bo=sinbo^8;
sin10bo=sinbo^10;

Ro=6380000;
%y0=2*pi*Ro*cos(B*pi/180)
%y0=0;

%geodetic reference system 1980 wgs84
% Google Maps use WGS84 https://productforums.google.com/forum/?hl=en#!category-topic/maps/base-map-data/Y6VEhIHdcB4
A1=6378137;  A2=6356752.3142;
%Bessel
% A1=6377397.15508;   A2=6356078.96290;


E2=1-((A2^2)/(A1^2));

E4=E2^2;
E6=E2^3;
E8=E2^4;
E10=E2^5;

e2s2bo=E2*sin2bo;
le2s2bo=1-E2*sin2bo;
le2=1-E2;

x10=(A1*cosbo)/(1-E2*sin2bo)^.5;
x11=(-A1*(1-E2)*sinbo)  /(1-E2*sin2bo )^1.5;
x30=A1*cosbo*(1-2*sin2bo+E2*sin4bo)/(6*(1-E2)*(1-E2*sin2bo)^.5);
x12=(-A1*(1-E2)*cosbo*(1+2*E2*sin2bo)) / (2*(1-E2*sin2bo)^2.5);
x31=((-A1*sinbo)/(6*(1-E2)*(1-E2*sin2bo )^1.5))*(5-E2-6*sin2bo*(1+E2)+3*E2*sin4bo*(3+E2)-4*E4*sin6bo);
x13=((A1*(1-E2)*sinbo)/(6*(1-E2*sin2bo )^3.5))*(1-9*E2+2*E2*sin2bo*(5-3*E2)+4*E4*sin4bo);
x50=((A1*cosbo)/(120*(1-E2)^3*(1-E2*sin2bo )^0.5))*(5-E2-4*sin2bo*(7+4*E2)+2*sin4bo*(12+43*E2+13*E4)-4*E2*sin6bo*(18+25*E2+3*E4)+E4*sin8bo*(77+39*E2)-28*E6*sin10bo);
x32=((-A1*cosbo)/(120*(1-E2)*(1-E2*sin2bo)^2.5))*(5-E2-2*sin2bo*(9+4*E2+E4)+15*E2*sin4bo*(3+E2)-2*E4*sin6bo*(23+3*E2)+16*E6*sin8bo);
x14=((A1*(1-E2)*cosbo)/(24*(1-E2*sin2bo )^4.5))*( 1-9*E2+36*E2*sin2bo*(1-2*E2)+12*E4*sin4bo*(5-2*E2)+8*E6*sin6bo );

y01=(A1*le2)/(le2s2bo^1.5);
y20=(A1*cosbo*sinbo)/(2*(le2s2bo)^0.5);
y02=(3*A1*E2*le2*cosbo*sinbo)/(2*(le2s2bo)^2.5);
y21=(A1*(1-2*sin2bo+E2*sin4bo))/(2*(le2s2bo)^2.5);
y03=(A1*E2*le2)/(2*(le2s2bo)^3.5)  *  (1-2*sin2bo*(1-2*E2)-3*E2*sin2bo);
y40=((A1*cosbo*sinbo)/(24*le2*le2s2bo^0.5)) * (5-E2-6*sin2bo*(1+E2)+3*(E2^1.5)*sin4bo*(3+E2)-4*E4*(sinbo^5)  );
y22=((-A1*cosbo*sinbo)/(4*le2s2bo^2.5)) * (4-3*E2-2*E2*sin2bo+E4*sin4bo);
y04=((-A1*E2*le2*cosbo*sinbo)/(8*le2s2bo^4.5)) * (4-15*E2+2*E2*sin2bo*(11-10*E2)+9*E2*sin4bo);
y41=(A1)/(24*le2^2*le2s2bo^1.5)  * ( 5-E2-4*sin2bo*(7+4*E2)+2*sin4bo*(12+43*E2+13*E4) -4*E2*sin6bo*(18+25*E2+3*E4)+E4*sin8bo*(77+39*E2)-28*E6*sin10bo);
y23=(-A1)/(12*le2s2bo^3.5)  *  (4-3*E2-4*sin2bo*(2-4*E2+3*E4)-2*E2*sin4bo*(2-5*E2)-4*E4*sin6bo+E6*sin8bo);
y05=(-A1*E2*le2)/(40*le2s2bo^5.5)  *  (4-15*E2-4*sin2bo*(2-32*E2+45*E4)-2*E2*sin4bo*(38-181*E2+60*E4)-4*E4*sin6bo*(41-34*E2)-27*(E2^2.5)*(sinbo^5)) ;

xlb=x10*l+x11*l*b+x30*l^3+x12*l*b^2+x31*l^3*b+x13*l*b^3+x50*l^5+x32*l^3*b^2+x14*l*b^4;

ylb=y01*b+y20*l^2+y02*b^2+y21*l^2*b+y03*b^3+y40*l^4+y22*l^2*b^2+y04*b^4+y41*l^4*b+y23*l^2*b^3+y05*b^5 ;

y0=A1*( (B0*pi/180) * (1-E2/4-3*E4/64-5*E6/256-175*E8/16384)-(3/8)*E2*(1+E2/4+15*E4/128-35*E6/512)*sin(2*B0*pi/180)+(15/256)*E4*(1+3*E2/4+35*E4/64)*sin(4*B0*pi/180)-(35/3072)*E6*(1+5*E2/4)*sin(6*B0*pi/180)+(315/131072)*E8*sin(8*B0*pi/180));


x=xlb+1000000*L0/3+500000;
y=y0+ylb;