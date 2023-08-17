% Body + NoseFins + TailFins + NoseFins
%% Variation of Moment coefficient w.r.t. Elevator deflection
format long g

g = 9.81;
Cl = [];
Alpha = [];
C_m = [];

%Balance load matrix
CM = [63.080043	0.144499	-0.206795	1.354260	1.630051	4.275882
    -0.123649	150.309342	0.592082	-0.725847	0.030856	0.393628
    0.024098	-0.689773	151.831777	0.096597	-0.571799	-4.414537
    0.152944	-2.334107	0.037781	77.595997	0.445712	5.841784
    -0.006030	0.114642	-0.574072	-0.065800	79.176337	0.322466
    0.047415	0.466131	0.099431	0.208017	0.190935	44.877349];

% Model reference parameters
% Sref(sqm)	        Chord (m)	    Span (m)
Sref = 0.009677;	c = 0.111;	    b = 0.111;

% Moment reference point (nose)
% X (m)		Y (m)		Z (m)
% 0	        0	        0

% Balance center location
%X (m)		Y (m)		Z (m)
X =-0.465;	Y = 0;	    Z = 0;

% Nowind data
%     Ax		    N1		    N2		    S1		    S2		    Rm
NW = [-0.000285	-0.000793	-0.000714	0.001382	0.001304	-0.000896 % -10
      0.000770	0.000858	0.000866	0.000727	-0.000826	0.000939  % -05
      -0.000060	0.000720	-0.000535	-0.000702	0.000175	-0.000938 % 00
      0.000584	0.000924	0.001203	0.000854	-0.001032	0.000897  % 05
      0.001288	0.001094	0.000372	0.000063	-0.000004	0.001433];% 10

Pitchi = 3.947917;

for i = 1:3
if i == 1
% Data for Speed 1
%       Dyn. head	Elevator	Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [956.613636	-10     	0.000000	0.000000	0.008969	0.014374	-0.016207	0.003729	0.000641	-0.000555
       962.477273	-05     	0.000000	0.000000	0.007567	0.015164	-0.010113	0.005263	-0.003927	0.001079
       980.386364	 00      	0.000000	0.000000	0.006000	0.011205	-0.004648	0.002760	-0.002287	-0.000951
       963.318182	 05      	0.000000	0.000000	0.007266	0.006149	0.005136	0.005586	-0.005078	0.000798
       951.579545	 10      	0.000000	0.000000	0.009792	0.000265	0.013822	0.005092	-0.004629	0.001546];


end
if i == 2
%Data for Speed 2
%       Dyn. head	Elevator	Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [1484.375000	-10     	0.000000	0.000000	0.013694	0.024866	-0.026298	0.004252	0.000996	-0.000635
       1497.863636	-05     	0.000000	0.000000	0.011547	0.021348	-0.014914	0.008450	-0.005916	0.001035
       1503.238636	 00      	0.000000	0.000000	0.009064	0.016918	-0.007061	0.004863	-0.003809	-0.001007
       1492.897727	 05      	0.000000	0.000000	0.010794	0.009277	0.006837	0.008420	-0.007463	0.000730
       1514.727273	 10      	0.000000	0.000000	0.014524	0.000673	0.020373	0.008270	-0.007440	0.001494];


end
if i == 3
% Data for Speed 3
%       Dyn. head	Pitch		Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [2155.170455	-10     	0.000000	0.000000	0.019466	0.037130	-0.038417	0.004791	0.001493	-0.000757
       2148.931818	-05     	0.000000	0.000000	0.016331	0.034577	-0.025182	0.011087	-0.007623	0.000859
       2155.477273	 00     	0.000000	0.000000	0.012750	0.024123	-0.010274	0.007539	-0.005805	-0.001135
       2150.988636	 05     	0.000000	0.000000	0.015385	0.012956	0.008959	0.012206	-0.010416	0.000700
       2152.750000	 10     	0.000000	0.000000	0.019810	0.001425	0.027528	0.011871	-0.010392	0.001382];

end

Le = 0.884;
% Transformation Matrix from CG to Body center
x = 0.2126; % Xnp --> 0.605 (Neutral point from nose)
y = 0;
z = 0;
TM =[-1,  0,  0,  0,  0,  0
      0,  1,  0,  0,  0,  0
      0,  0,  -1, 0,  0,  0
      0,  -z, y,  1,  0,  0
      z,  0,  -x, 0,  1,  0
      -y, x,  0,  0,  0,  1];

    for j = 1:5
        aj = [DSi(j,5) ; DSi(j,6) ; DSi(j,7) ; DSi(j,8) ; DSi(j,9) ; DSi(j,10)];
        % Convert the obtained normalized voltage signals to kg
        Aj = CM * (aj-(NW(j,1:6))');
        % Calculate the forces and moments about the balance center
        fm = [Aj(1) ; (Aj(4)+Aj(5)) ; (Aj(2)+Aj(3)) ; Aj(6) ; (Aj(2)-Aj(3))*0.065 ; (Aj(4)-Aj(5))*0.065].*g;
        % Transform the forces and moments to the body axis
        FM = [Aj(1) ; (Aj(4)+Aj(5)) ; (Aj(2)+Aj(3)) ; Aj(6) ; (Aj(2)-Aj(3))*0.065 ; (Aj(4)-Aj(5))*0.065].*g; 
        % Transform the forces and moments to the C.G of the flight vehicle
        FMcg = TM * FM;
        % The variation of longitudinal force and moment coefficients
        Cf = (1/(DSi(j,1)*Sref)).*[FMcg(1:3)];
        Cm = (1/(DSi(j,1)*Sref*Le)).*[FMcg(4:6)];
        
        Ele = (DSi(j,2));
        % Variation of Aerodynamic force coefficients with angle of attack
        Aero_TM = [ sind(Pitchi)   0   -cosd(Pitchi)
                   -cosd(Pitchi)   0   -sind(Pitchi)
                   0               1   0];
        
        Aero_Coeff = Aero_TM*Cf;
        
        % Recording The Data
        Alpha = [Alpha; Ele'];
        Cl  = [Cl; Aero_Coeff'];
        C_m = [C_m; Cm'];

    end
end

        A1= Alpha(1:5);
        B1= Cl(1:5,1);
        D1= C_m(1:5,2);

        A2= Alpha(6:10);
        B2= Cl(6:10,1);
        D2= C_m(6:10,2);

        A3= Alpha(11:15);
        B3= Cl(11:15,1);
        D3= C_m(11:15,2);
        
figure(1)        
        hold on
        plot(A1,B1,'r-*',A2,B2,'g--o',A3,B3,'b-+')   
        hold off       
        xlabel('\delta_e')
        ylabel('Cl')
        title('Cl vs \delta_e')
        legend('V1=40m/s','V2=50m/s','V3=60m/s')
        grid on

figure(2)
        hold on
        plot(A1,D1,'r-*',A2,D2,'g--o',A3,D3,'b-+')   
        hold off
        xlabel('\delta_e')
        ylabel('Cm')
        title('Cm vs \delta_e')
        legend('V1=40m/s','V2=50m/s','V3=60m/s')
        grid on

% Calculating Slope of Moment Coefficient
%For V1
        cc1=0;
        Cm_Alpha1i = 0;
        Cl_Alpha1i = 0;
        for l = 2:1:5
        D_Alpha1 =  A1(l)-A1(l-1) ;
        D_Cl1 = B1(l)-B1(l-1);
        D_Cm1 = D1(l)-D1(l-1);
        Cm_alpha1 = D_Cm1/D_Alpha1;
        Cl_alpha1 = D_Cl1/D_Alpha1;
        Cl_Alpha1i = Cl_Alpha1i + Cl_alpha1;
        Cm_Alpha1i = Cm_Alpha1i + Cm_alpha1;
        cc1 = cc1+1;
        end        
        F_Cm_Alpha1 = Cm_Alpha1i/cc1;
        F_Cl_Alpha1 = Cl_Alpha1i/cc1;
        disp('For V1 = 40 m/s')
        disp('Average value of Cl_delta_e')
        disp(F_Cl_Alpha1)
        disp('Average value of Cm_delta_e')
        disp(F_Cm_Alpha1)

%For V2
        cc2=0;
        Cm_Alpha2i = 0;
        Cl_Alpha2i = 0;
        for l = 2:1:5
        D_Alpha2 =  A2(l)-A2(l-1) ;
        D_Cl2 = B2(l)-B2(l-1);
        D_Cm2 = D2(l)-D2(l-1);
        Cm_alpha2 = D_Cm2/D_Alpha2;
        Cl_alpha2 = D_Cl2/D_Alpha2;
        Cl_Alpha2i = Cl_Alpha2i + Cl_alpha2;
        Cm_Alpha2i = Cm_Alpha2i + Cm_alpha2;
        cc2 = cc2+1;
        end
        F_Cm_Alpha2 = Cm_Alpha2i/cc2;
        F_Cl_Alpha2 = Cl_Alpha2i/cc2;
        disp('For V2 = 50 m/s')
        disp('Average value of Cl_delta_e')
        disp(F_Cl_Alpha2)
        disp('Average value of Cm_delta_e')
        disp(F_Cm_Alpha2)

%For V3
        cc3=0;
        Cm_Alpha3i = 0;
        Cl_Alpha3i = 0;
        for l = 2:1:5
        D_Alpha3 =  A3(l)-A3(l-1) ;
        D_Cl3 = B3(l)-B3(l-1);
        D_Cm3 = D3(l)-D3(l-1);
        Cm_alpha3 = D_Cm3/D_Alpha3;
        Cl_alpha3 = D_Cl3/D_Alpha3;
        Cl_Alpha3i = Cl_Alpha3i + Cl_alpha3;
        Cm_Alpha3i = Cm_Alpha3i + Cm_alpha3;
        cc3 = cc3+1;
        end
        F_Cm_Alpha3 = Cm_Alpha3i/cc3;
        F_Cl_Alpha3 = Cl_Alpha3i/cc3;
        disp('For V3 = 60 m/s')
        disp('Average value of Cl_delta_e')
        disp(F_Cl_Alpha3)
        disp('Average value of Cm_delta_e')
        disp(F_Cm_Alpha3)