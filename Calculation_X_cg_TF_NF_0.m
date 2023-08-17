% 6694	5/8/2010 10:50:12 AM
% Bomb Model-II Test with Body + NoseFins + TailFins + delta 0
%% Variation of Moment coefficient w.r.t. Center of gravity 
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
    0.047415	0.466131	0.099431	0.208017	0.190935	44.877349 ];

% Model reference parameters
% Sref(sqm)	        Chord (m)	    Span (m)
Sref = 0.009677;	c = 0.111;	    b = 0.111;

% Moment reference point (nose)
% X (m)		Y (m)		Z (m)
% 0	        0	        0

% Balance center location
%X (m)		Y (m)		Z (m)
X =-0.465;	Y = 0;	    Z = 0;
X_bc = 0.465 ; % distance from Balence Center to Nose of missile

% Nowind data
%     Ax		    N1		    N2		    S1		    S2		    Rm
NW = [-0.000060	    0.000720	-0.000535	-0.000702	0.000175	-0.000938];

for i = 1:3
if i == 1
% Data for Speed 1
DyHi = 980.386364;
%       Dyn. head	Pitch		Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [ 980.386364	-9.947917	0.000000	0.000000	0.005824	-0.018730	-0.003349	-0.000120	-0.000445	-0.000447
        980.386364	-8.052083	0.000000	0.000000	0.005779	-0.018470	0.000780	0.000486	-0.001005	-0.000446
        980.386364	-6.052083	0.000000	0.000000	0.005848	-0.015727	0.002675	0.001113	-0.001597	-0.000533
        980.386364	-4.052083	0.000000	0.000000	0.005919	-0.010716	0.002267	0.001595	-0.002072	-0.000610
        980.386364	-2.052083	0.000000	0.000000	0.006174	-0.005618	0.001225	0.002030	-0.002352	-0.000731
        980.386364	-0.052083	0.000000	0.000000	0.006382	0.000178	-0.000973	0.002233	-0.002247	-0.000805
        980.386364	0.052083	0.000000	0.000000	0.006418	0.000622	-0.001537	0.002333	-0.002409	-0.000828
        980.386364	1.947917	0.000000	0.000000	0.006301	0.005986	-0.003291	0.002541	-0.002314	-0.000872
        980.386364	3.947917	0.000000	0.000000	0.006000	0.011205	-0.004648	0.002760	-0.002287	-0.000951
        980.386364	5.947917	0.000000	0.000000	0.005786	0.016277	-0.005346	0.003033	-0.002387	-0.001062
        980.386364	7.947917	0.000000	0.000000	0.005650	0.020372	-0.004799	0.003231	-0.002487	-0.001154
        980.386364	9.958333	0.000000	0.000000	0.005792	0.020911	-0.001107	0.003629	-0.002522	-0.001206
        980.386364	11.947917	0.000000	0.000000	0.005782	0.019842	0.004121	0.003848	-0.002633	-0.001215
        980.386364	13.947917	0.000000	0.000000	0.005569	0.018813	0.009425	0.004079	-0.002668	-0.001259
        980.386364	15.947917	0.000000	0.000000	0.005211	0.017681	0.014868	0.004138	-0.002628	-0.001260
        980.386364	17.947917	0.000000	0.000000	0.004741	0.017236	0.020059	0.004069	-0.002389	-0.001278
        980.386364	19.947917	0.000000	0.000000	0.004227	0.016561	0.025673	0.003951	-0.002322	-0.001245
        980.386364	21.947917	0.000000	0.000000	0.003654	0.016103	0.031224	0.003849	-0.002502	-0.001226
        980.386364	23.947917	0.000000	0.000000	0.003050	0.016164	0.036243	0.003657	-0.002624	-0.001240
        980.386364	25.947917	0.000000	0.000000	0.002483	0.016457	0.040852	0.004120	-0.003231	-0.001296];%20
end
if i == 2
%Data for Speed 2
DyHi = 1503.238636;
%       Dyn. head	Pitch		Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [ 1503.238636	-9.947917	0.000000	0.000000	0.008672	-0.029165	-0.004914	-0.000077	-0.000414	-0.000311
        1503.238636	-8.052083	0.000000	0.000000	0.008730	-0.029110	0.001760	0.000991	-0.001441	-0.000276
        1503.238636	-6.052083	0.000000	0.000000	0.008875	-0.024977	0.004733	0.001999	-0.002528	-0.000371
        1503.238636	-4.052083	0.000000	0.000000	0.009019	-0.017260	0.004091	0.002831	-0.003249	-0.000504
        1503.238636	-2.052083	0.000000	0.000000	0.009410	-0.009170	0.002185	0.003521	-0.003766	-0.000601
        1503.238636	-0.052083	0.000000	0.000000	0.009866	-0.000294	-0.001190	0.003980	-0.003718	-0.000728
        1503.238636	0.052083	0.000000	0.000000	0.009806	0.000650	-0.002243	0.004068	-0.003786	-0.000868
        1503.238636	1.947917	0.000000	0.000000	0.009634	0.008807	-0.004948	0.004488	-0.003806	-0.000846
        1503.238636	3.947917	0.000000	0.000000	0.009064	0.016918	-0.007061	0.004863	-0.003809	-0.001007
        1503.238636	5.947917	0.000000	0.000000	0.008705	0.024919	-0.008314	0.005317	-0.004052	-0.001171
        1503.238636	7.947917	0.000000	0.000000	0.008391	0.031412	-0.007677	0.005684	-0.004145	-0.001335
        1503.238636	9.947917	0.000000	0.000000	0.008504	0.032596	-0.002311	0.006242	-0.004257	-0.001430
        1503.238636	11.947917	0.000000	0.000000	0.008532	0.031188	0.005413	0.006732	-0.004278	-0.001513
        1503.238636	13.947917	0.000000	0.000000	0.008216	0.029086	0.013746	0.006902	-0.004357	-0.001548
        1503.238636	15.947917	0.000000	0.000000	0.007643	0.027572	0.021870	0.007058	-0.004498	-0.001615
        1503.238636	17.947917	0.000000	0.000000	0.006965	0.026547	0.030059	0.006797	-0.003870	-0.001538
        1503.238636	19.947917	0.000000	0.000000	0.006100	0.025577	0.038565	0.006591	-0.003663	-0.001523
        1503.238636	21.947917	0.000000	0.000000	0.005224	0.024703	0.047349	0.006408	-0.003872	-0.001488
        1503.238636	23.947917	0.000000	0.000000	0.004296	0.024660	0.055195	0.006309	-0.004261	-0.001527
        1503.238636	25.947917	0.000000	0.000000	0.003412	0.024635	0.062857	0.006864	-0.005023	-0.001578];
end
if i == 3
% Data for Speed 3
DyHi = 2155.477273;
%       Dyn. head	Pitch		Yaw		    Roll		Ax		    N1		    N2		    S1		    S2		    Rm
DSi = [ 2155.477273	-9.947917	0.000000	0.000000	0.012224	-0.043297	-0.006181	-0.000401	0.000118	-0.000044
        2155.477273	-8.052083	0.000000	0.000000	0.012382	-0.043300	0.003434	0.001114	-0.001542	0.000027
        2155.477273	-6.052083	0.000000	0.000000	0.012492	-0.037048	0.007539	0.002741	-0.003228	-0.000118
        2155.477273	-4.052083	0.000000	0.000000	0.012852	-0.026257	0.006930	0.003982	-0.004439	-0.000348
        2155.477273	-2.052083	0.000000	0.000000	0.013389	-0.014235	0.003861	0.005388	-0.005434	-0.000540
        2155.477273	-0.052083	0.000000	0.000000	0.014060	-0.001266	-0.001222	0.006043	-0.005342	-0.000732
        2155.477273	0.052083	0.000000	0.000000	0.013911	0.000759	-0.003306	0.006368	-0.005640	-0.000822
        2155.477273	1.947917	0.000000	0.000000	0.013731	0.012064	-0.006754	0.006798	-0.005611	-0.000923
        2155.477273	3.947917	0.000000	0.000000	0.012750	0.024123	-0.010274	0.007539	-0.005805	-0.001135
        2155.477273	5.947917	0.000000	0.000000	0.012157	0.035754	-0.012153	0.008187	-0.006238	-0.001378
        2155.477273	7.947917	0.000000	0.000000	0.011812	0.045210	-0.011220	0.008632	-0.006314	-0.001616
        2155.477273	9.958333	0.000000	0.000000	0.011769	0.047922	-0.004543	0.009309	-0.006195	-0.001753
        2155.477273	11.947917	0.000000	0.000000	0.011876	0.046518	0.005962	0.009825	-0.005895	-0.001883
        2155.477273	13.947917	0.000000	0.000000	0.011479	0.042350	0.018673	0.010167	-0.006151	-0.001922
        2155.477273	15.947917	0.000000	0.000000	0.010666	0.040216	0.030442	0.010338	-0.006394	-0.002040
        2155.477273	17.937500	0.000000	0.000000	0.009717	0.038464	0.042449	0.010045	-0.005459	-0.001934
        2155.477273	19.947917	0.000000	0.000000	0.008543	0.036833	0.054952	0.009845	-0.005513	-0.001877
        2155.477273	21.947917	0.000000	0.000000	0.007323	0.035552	0.067531	0.008894	-0.005253	-0.001782
        2155.477273	23.947917	0.000000	0.000000	0.006033	0.035056	0.079428	0.009108	-0.006263	-0.001820
        2155.477273	25.947917	0.000000	0.000000	0.004720	0.034759	0.090736	0.009236	-0.007105	-0.001772];

end
Le = 0.884;
% Transformation Matrix from CG to Body center
for x = -0.4:0.1:0.4  % -0.2:0.01:-0.12
% Xnp --> 0.605 (x=-0.14)
y = 0;
z = 0;
TM =[-1,  0,  0,  0,  0,  0
      0,  1,  0,  0,  0,  0
      0,  0,  -1, 0,  0,  0
      0,  -z, y,  1,  0,  0
      z,  0,  -x, 0,  1,  0
      -y, x,  0,  0,  0,  1];

Alpha = [];
C_m = [];
    for j = 1:20
        aj = [DSi(j,5) ; DSi(j,6) ; DSi(j,7) ; DSi(j,8) ; DSi(j,9) ; DSi(j,10)];
        % Convert the obtained normalized voltage signals to kg
        Aj = CM * (aj-NW);
        % Calculate the forces and moments about the balance center.
        fm = [Aj(1) ; (Aj(4)+Aj(5)) ; (Aj(2)+Aj(3)) ; Aj(6) ; (Aj(2)-Aj(3))*0.065 ; (Aj(4)-Aj(5))*0.065].*g;
        % Transform the forces and moments to the body axis
        FM = [Aj(1) ; (Aj(4)+Aj(5)) ; (Aj(2)+Aj(3)) ; Aj(6) ; (Aj(2)-Aj(3))*0.065 ; (Aj(4)-Aj(5))*0.065].*g; 
        % Transform the forces and moments to the C.G of the flight vehicle
        FMcg = TM * FM;
        % The variation of longitudinal force and moment coefficients
        Cf = (1/(DyHi*Sref)).*[FMcg(1:3)];
        Cm = (1/(DyHi*Sref*Le)).*[FMcg(4:6)];

        Pitch = (DSi(j,2));
        % Variation of Aerodynamic force coefficients with angle of attack
        Aero_TM = [ sind(Pitch)   0   -cosd(Pitch)
                   -cosd(Pitch)   0   -sind(Pitch)
                   0              1   0];
  
        Aero_Coeff = Aero_TM*Cf;
        
        % Recording The Data
        Alpha = [Alpha; Pitch];        
        C_m = [C_m; Cm'];

 
    end
    % Iterating the x (from Balance center) value to obtain Cm_\alpha
   figure(i)
        A = Alpha(:,1);
        ct = C_m(:,2);
        hold on
        plot(A,ct)
        hold off
        xlabel('\alpha')
        ylabel('Cm')
        ylim([-6 6])
        title('Iterating the x to obtain Cm_\alpha')
        %legend('Xcg1=0.10','Xcg2=0.11','Xcg3=0.12','Xcg4=0.13','Xcg5=0.14','Xcg6=0.15','Xcg7=0.16','Xcg8=0.17','Xcg9=0.18')
        legend('Xcg1=-0.4','Xcg2=-0.3','Xcg3=-0.2','Xcg4=-0.1','Xcg5=0','Xcg6=0.1','Xcg7=0.2','Xcg8=0.3','Xcg9=0.4')
        grid on    
end
end
% Most horizontal line of the Output Graph is for x = -0.2 to -0.1
% Let Neutral Point be
X_bnp = 0.08; %distance between Balence Center and Neutral Point
disp(['distance between Balence Center and Neutral Point is--> ', num2str(X_bnp) ' m'])

X_np = X_bc - X_bnp; % distance of Balence Center from Nose of missile
disp(['distance of Neutral Point from Nose of missile is--> ', num2str(X_np) ' m'])

% Assuming Static Margine ( = 15%Le ) 
disp('Assuming Static Margine 15% of length of missile')
SM = 0.15*Le;
disp(['Static Margine is--> ', num2str(SM) ' m'])

% distance between Balence Center and C.G
X_bg = X_bnp + SM;
disp(['distance between Balence Center and C.G is--> ', num2str(X_bg) ' m'])

X_cg = X_np-SM;
disp(['distance of C.G. from Nose of missile is--> ', num2str(X_cg) ' m'])