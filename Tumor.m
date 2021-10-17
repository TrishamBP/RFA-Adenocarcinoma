%% Applied Finite Element Methods in Engineering.
%% Trisham Bharat Patil-688519827
%% HW-6/Final Assignment.
%% Data Inpust Task.(DIT)
clear all
clc
syms Source % S is the source value which can be inputted by the user.
theta=0.5;
time=0; dt=0.2;
ngp=2; % Input of how many gauss points do we want.
[Gauss, WGauss] = GPoints(ngp);
[fname,fpath]=uigetfile('*.txt','Please select the Quad input file');
filename=sprintf('%s%s',fpath,fname);
fid = fopen(filename, 'r');
% Setting up the vector space.
nn = fscanf(fid,'%d',1);
ne = fscanf(fid,'%d',1);
nT1 = fscanf(fid,'%d',1);
nSrc = fscanf(fid,'%d',1);
x = zeros(nn,1);  y = x; 
BC1 = zeros(nT1,2);
src = zeros(nSrc,2); % Heat Source.
Source=input('Please enter the value of the source term:  ');

for i = 1:nn
    node_num = fscanf(fid, '%d',1);
    x(node_num) = fscanf(fid, '%f',1);
    y(node_num) = fscanf(fid, '%f',1);
end
for L = 1:ne 
    elem_num = fscanf(fid, '%d',1);
    El(elem_num,:) = fscanf(fid, '%d',4);
    mp_num(elem_num) = fscanf(fid, '%d',1);% Reading in the MP number.
end
for i = 1:nT1 
    bdy_num(i)=fscanf(fid, '%d',1);
end
for i = 1:nSrc 
    Loc(i) = fscanf(fid, '%d',1);
    S(Loc) = Source;  
end    
fclose(fid)
% Code to look at the geometry.
coords = zeros(nn,3); 
Elcolor = ['b';'m';'r';'g';'c';'y'];
coords(:,1) = x;coords(:,2) = y;coords(:,3) = 0;
for L = 1:ne 
    patch('Faces',El(L,1:4),'Vertices',coords,'FaceVertexCData',0,'FaceColor',Elcolor(mp_num(L)));
end
%% Data Analysis Task.(DAT).
Lhs = zeros(nn,nn);
Rhs= zeros(nn,nn);
R=zeros(nn,1);xL = zeros(4,1); yL = xL;SL = xL; D=zeros(6,1);
% Applying the IC.
T=zeros(nn,1);
vol_sp=[0.8;0.6;0.5;0.6;1.0;0.6];
conductivity=[0.002;0.003;0.003;0.006;0.004;0.006];

for L = 1:ne
    xL= x(El(L,:)); 
    yL= y(El(L,:));
    K= conductivity(mp_num(L));
    c=vol_sp(mp_num(L));
    D = K/c;
    for gp1 = 1:ngp
        Xi = Gauss(gp1); 
        wt1 = WGauss(gp1);
        for gp2 = 1:ngp
            Eta = Gauss(gp2); 
            wt2 = WGauss(gp2);
            [N, dNdx, dNdy, DJ] = gp2d(xL,yL,Xi,Eta);W = N; dWdx = dNdx; dWdy = dNdy;
            for i = 1:4                    
                iG = El(L,i);
                for j=1:4
                    jG=El(L,j);
                    Term_1=W(i)*N(j);
                    Term_2=(D*((dWdx(i)*dNdx(j))+(dWdy(i)*dNdy(j))))+(W(i)*N(j));
                    Lhs(iG,jG)=Lhs(iG,jG)+wt1*wt2*DJ*(Term_1+((theta*dt)*Term_2));
                    Rhs(iG,jG)=Rhs(iG,jG)+wt1*wt2*DJ*(Term_1-(((1-theta)*dt)*Term_2));
                end
            end
        end
    end
end

Converge=1.0e-4; Diff=1;
count=0; increment=1;
%% This part of the code will plot for t=1.5.
while time<=1.5
    time=time+dt;  
    Rside=Rhs*T;
    Rside=Rside+R;
    
    % Calculating the type 1 BC.
    for i=1:nT1
        Loc=bdy_num(i);
        Lhs(Loc,:)=0;
        Lhs(Loc,Loc)=1;
        Rside(Loc)=1.0;
    end  
    Tnew=Lhs\Rside;
    
    Diff=max(abs(Tnew-T))
    
    T=Tnew;
    %% Data Output Task (DOT).
    coords(:,3) = Tnew;
    patch('Faces',El,'Vertices',coords,'FaceVertexCData',Tnew,'FaceColor','interp');
end
fid = fopen('HW6Temperatures-c.txt','w');
fprintf(fid,'  These are the Temperature values \n');
fprintf(fid,'Temperature     X        Y \n');
for i = 1:nn
    fprintf(fid,'   %d      %8.4f    %8.4f \n',x(i),y(i),T(i));
end
%% This commented out code will give the time to reach steady state.
% while Diff>Converge
%     time=time+dt; count=count+1;
%     Rside=Rhs*T;
%     Rside=Rside+R;
%     
%     % Calculating the type 1 BC.
%     for i=1:nT1
%         Loc=bdy_num(i);
%         Lhs(Loc,:)=0;
%         Lhs(Loc,Loc)=1;
%         Rside(Loc)=1.0;
%     end  
%     Tnew=Lhs\Rside;
%     
%     Diff=max(abs(Tnew-T))
%     
%     T=Tnew;
%     %% Data Output Task (DOT).
%     coords(:,3) = Tnew;
%     patch('Faces',El,'Vertices',coords,'FaceVertexCData',Tnew,'FaceColor','interp');
% end

%% List of functions.
function [N,DNX,DNY,DJ] = gp2d(X,Y,XI,ETA)
% GAUSS POINT 2 DIMENSIONAL LINEAR ELEMENTS
% THIS PROGRAM WAS LAST UPDATED ON 10-27-83
% BY JOHN M. SULLIVAN, JR.   THE PROGRAM PERFORMS
% GAUSS QUADRATURE. LINEAR IN (XI) AND (ETA) DIRECTIONS
%
%      *----------*
%      %4        3%
%      %          %   NODE ORDERING IS CCW, 
%      %1        2%    FIRST NODE LOCATION
%      *----------*     IS NOT CRITICAL.
%
% STANDARD SERINDEPITY FINITE FORMULATIONS
N = zeros(4,1);    DJA = zeros(4,1);
DNX = zeros(4,1);  DNY = zeros(4,1);
DNE = zeros(4,1);  DNXI = zeros(4,1);
%
% BASIS FUNCTIONS FOR UNKNOWNS (I.E. TEMPERATURE)
N(1) = 0.25*(1.-XI)*(1.-ETA);  N(2) = 0.25*(1.+XI)*(1.-ETA);
N(3) = 0.25*(1.+XI)*(1.+ETA);  N(4) = 0.25*(1.-XI)*(1.+ETA);
%
% DERIVATIVES OF BASIS FUNCTIONS W.R.T. XI DIRECTION
DNXI(1) = -0.25*(1.-ETA);
DNXI(2) = 0.25*(1.-ETA);
DNXI(3) = 0.25*(1+ETA);
DNXI(4) = -0.25*(1.+ETA);
%
% DERIVATIVES OF BASIS FUNCTIONS W.R.T. ETA DIRECTION
DNE(1) = -0.25*(1.-XI);
DNE(2) = -0.25*(1.+XI);
DNE(3) = 0.25*(1+XI);
DNE(4) = 0.25*(1-XI);
%
% JACOBIAN
for I=1:4
    DJA(1)=X(I)*DNXI(I) + DJA(1);
    DJA(2)=Y(I)*DNXI(I) + DJA(2);
    DJA(3)=X(I)*DNE(I)  + DJA(3);
    DJA(4)=Y(I)*DNE(I)  + DJA(4);
end
DJ = DJA(1)*DJA(4) - DJA(2)*DJA(3);
%
% DERIVATIVES W.R.T. X AND Y
for I=1:4
    DNX(I) = (DJA(4)*DNXI(I)-DJA(2)*DNE(I))/DJ;
    DNY(I) = (-DJA(3)*DNXI(I) + DJA(1)*DNE(I))/DJ;
end
end
function [X, W]= GPoints(n)
%% This function receives the number of desired Gauss Points and
%   returns the specific locations "X" and weights "W"
%   John Sullivan  2014

X = zeros(n,1); W=zeros(n,1);
switch n
    case 1
        W(1) = 	2.0;
        X(1) = 0.0;
    case 2
        W(1) = 1; W(2) = 1;
        X(1) = -0.5773502691896257; X(2) = -X(1);
    case 3
        W(2) = 	0.8888888888888888;
        X(2) = 0.0000;
        W(1) = 0.5555555555555556; W(3)=W(1);
        X(1) = -0.7745966692414834; X(3) = -X(1);
    case 4
        W(1) =	0.6521451548625461; W(2) = W(1);
        X(1) = -0.3399810435848563; X(2) = -X(1);
        W(3) = 	0.3478548451374538; W(4) = W(3);
        X(3) = -0.8611363115940526; X(4) = -X(3);
    case 5
        W(3) =	0.5688888888888889;
        X(3) = 0.0000000000000000;
        W(2) =	0.4786286704993665; W(4) = W(2);
        X(2) = -0.5384693101056831; X(4) = -X(2);
        W(1) = 0.2369268850561891;  W(5) = W(1);
        X(1) = -0.9061798459386640; X(5) = -X(1);
    case 6
        W(1)= 0.3607615730481386; X(1)= 0.6612093864662645;
        W(2)= 0.3607615730481386; X(2)= -0.6612093864662645;
        W(3)= 0.4679139345726910; X(3)= -0.2386191860831969;
        W(4)= 0.4679139345726910; X(4)= 0.2386191860831969;
        W(5)= 0.1713244923791704; X(5)= -0.9324695142031521;
        W(6)= 0.1713244923791704; X(6)= 0.9324695142031521;
    case 7
        W(1)= 0.4179591836734694; X(1)=	0.0000000000000000;
        W(2)= 0.3818300505051189; X(2)=	0.4058451513773972;
        W(3)= 0.3818300505051189; X(3)=	-0.4058451513773972;
        W(4)= 0.2797053914892766; X(4)=	-0.7415311855993945;
        W(5)= 0.2797053914892766; X(5)=	0.7415311855993945;
        W(6)= 0.1294849661688697; X(6)=	-0.9491079123427585;
        W(7)= 0.1294849661688697; X(7)=	0.9491079123427585;
    case 8
        W(1) = 0.3626837833783620; X(1)= -0.1834346424956498;
        W(2) = 0.3626837833783620; X(2)= 0.1834346424956498;
        W(3) = 0.3137066458778873; X(3)= -0.5255324099163290;
        W(4) = 0.3137066458778873; X(4)= 0.5255324099163290;
        W(5) = 0.2223810344533745; X(5)= -0.7966664774136267;
        W(6) = 0.2223810344533745; X(6)= 0.7966664774136267;
        W(7) = 0.1012285362903763; X(7)= -0.9602898564975363;
        W(8) = 0.1012285362903763; X(8)= 0.9602898564975363;
    case 9
        W(1) = 0.3302393550012598; X(1) = 0.0000000000000000;
        W(2) = 0.1806481606948574; X(2) = -0.8360311073266358;
        W(3) = 0.1806481606948574; X(3) = 0.8360311073266358;
        W(4) = 0.0812743883615744; X(4) = -0.9681602395076261;
        W(5) = 0.0812743883615744; X(5) = 0.9681602395076261;
        W(6) = 0.3123470770400029; X(6) = -0.3242534234038089;
        W(7) = 0.3123470770400029; X(7) = 0.3242534234038089;
        W(8) = 0.2606106964029354; X(8) = -0.6133714327005904;
        W(9) = 0.2606106964029354; X(9) = 0.6133714327005904;
    case 10
        W(1) = 0.2955242247147529; X(1) = -0.1488743389816312;
        W(2) = 0.2955242247147529; X(2) = 0.1488743389816312;
        W(3) = 0.2692667193099963; X(3) = -0.4333953941292472;
        W(4) = 0.2692667193099963; X(4) = 0.4333953941292472;
        W(5) = 0.2190863625159820; X(5) = -0.6794095682990244;
        W(6) = 0.2190863625159820; X(6) = 0.6794095682990244;
        W(7) = 0.1494513491505806; X(7) = -0.8650633666889845;
        W(8) = 0.1494513491505806; X(8) = 0.8650633666889845;
        W(9) = 0.0666713443086881; X(9) = -0.9739065285171717;
        W(10) = 0.0666713443086881; X(10) = 0.9739065285171717;
    case 11
        W(1) = 0.2729250867779006; X(1) = 0.0000000000000000;
        W(2) = 0.2628045445102467; X(2) = -0.2695431559523450;
        W(3) = 0.2628045445102467; X(3) = 0.2695431559523450;
        W(4) = 0.2331937645919905; X(4) = -0.5190961292068118;
        W(5) = 0.2331937645919905; X(5) = 0.5190961292068118;
        W(6) = 0.1862902109277343; X(6) = -0.7301520055740494;
        W(7) = 0.1862902109277343; X(7) = 0.7301520055740494;
        W(8) = 0.1255803694649046; X(8) = -0.8870625997680953;
        W(9) = 0.1255803694649046; X(9) = 0.8870625997680953;
        W(10) = 0.0556685671161737; X(10) = -0.9782286581460570;
        W(11) = 0.0556685671161737; X(11) = 0.9782286581460570;
    otherwise
        warning('Invalid passing parameters');
end
end
