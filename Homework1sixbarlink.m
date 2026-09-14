% six bar linkage
% Static Equilibrium, Velocity, and Acceleration Analysis

clc;
clear;

% define the joints
A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];

% Define the lengths of the bars
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lDE = norm(E - D);
lEF = norm(F - E);
lFG = norm(G - F);
lCE=norm(E-C);
lBE=norm(E-B);
lEG=norm(E-G);

% Weight of Each Link
g = 9.81;

MassAB = 23.311;   J_AB  = 1148.57 * 0.00064516;
MassBC = 56.403;   J_BC  = 15251.887 * 0.00064516;   % "BEC" in code -> use BC's real part
MassDE = 98.436;   J_DE  = 80331.584 * 0.00064516;   % "CD" in code -> use DE's real part
MassEF = 46.971;   J_EF  = 8861.866 * 0.00064516;
MassFG = 102.262;  J_FG  = 90039.888 * 0.00064516;

WAB  = [0 -MassAB*g 0];
WBEC = [0 -MassBC*g 0];
WCD  = [0 -MassDE*g 0];
WEF  = [0 -MassEF*g 0];
WFG  = [0 -MassFG*g 0];


% center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];

% Applied Force
AppliedForce = [0 200 0];

% Static Equilibrium Condition for Link AB
eqn1 = ForceA + ForceB + WAB == 0;
eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;
eqn4 = cross(B - S2, -ForceB) + cross(C - S2, ForceC) + cross(E - S2, ForceE) == 0;

eqn5 = -ForceC + ForceD + WCD == 0;
eqn6 = cross(C - S3, -ForceC) + cross(D - S3, ForceD) == 0;

eqn7 = -ForceE + ForceF + WEF == 0;
eqn8 = cross(E - S4, -ForceE) + cross(F - S4, ForceF) == 0;

eqn9  = -ForceF + ForceG + WFG + AppliedForce == 0;
eqn10 = cross(F - S5, -ForceF) + cross(G - S5, ForceG) == 0;

% Solving the 10 equations
eqnMatrix = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];

StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax(1) = double(StaticSolution.FAx);
Force_Ay(1) = double(StaticSolution.FAy);

Force_Bx(1) = double(StaticSolution.FBx);
Force_By(1) = double(StaticSolution.FBy);

Force_Cx(1) = double(StaticSolution.FCx);
Force_Cy(1) = double(StaticSolution.FCy);

Force_Dx(1) = double(StaticSolution.FDx);
Force_Dy(1) = double(StaticSolution.FDy);

Force_Ex(1) = double(StaticSolution.FEx);
Force_Ey(1) = double(StaticSolution.FEy);

Force_Fx(1) = double(StaticSolution.FFx);
Force_Fy(1) = double(StaticSolution.FFy);

Force_Gx(1) = double(StaticSolution.FGx);
Force_Gy(1) = double(StaticSolution.FGy);

Input_Torque = double(StaticSolution.Tin);

disp('Static equilibrium forces (N):')
fprintf('Force A: [%.4f, %.4f]\n', Force_Ax(1), Force_Ay(1));
fprintf('Force B: [%.4f, %.4f]\n', Force_Bx(1), Force_By(1));
fprintf('Force C: [%.4f, %.4f]\n', Force_Cx(1), Force_Cy(1));
fprintf('Force D: [%.4f, %.4f]\n', Force_Dx(1), Force_Dy(1));
fprintf('Force E: [%.4f, %.4f]\n', Force_Ex(1), Force_Ey(1));
fprintf('Force F: [%.4f, %.4f]\n', Force_Fx(1), Force_Fy(1));
fprintf('Force G: [%.4f, %.4f]\n', Force_Gx(1), Force_Gy(1));
fprintf('Input Torque (Nm): %.4f\n', Input_Torque);

%% Angular velocity calculations
% Loop 1: A-B-C-D-A
syms wBEC wCD
omega_AB  = [0 0 23.14];
omega_BEC = [0 0 wBEC];
omega_CD  = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loop1Solution = solve(eqn11,[wBEC wCD]);

% angular velocities from loop solution
angularVelocity_BEC = double(loop1Solution.wBEC);
angularVelocity_CD  = double(loop1Solution.wCD);

omegaBEC = [0 0 angularVelocity_BEC];
omegaCD  = [0 0 angularVelocity_CD];

% Loop 2: D-C-E-F-G-D
syms wEF wFG
omega_EF = [0 0 wEF];   % FIX: was mistakenly [0 0 WEF] (the weight vector)
omega_FG = [0 0 wFG];

eqn12 = cross(omegaCD,D-C) + cross(omegaBEC,E-C) + cross(omega_EF,F-E) + cross(omega_FG,G-F) == 0;

loop2solution = solve(eqn12,[wEF wFG]);

% Angular velocities from loop solution
angularVelocity_EF(1) = double(loop2solution.wEF);
angularVelocity_FG(1) = double(loop2solution.wFG);

omegaEF = [0 0 angularVelocity_EF];
omegaFG = [0 0 angularVelocity_FG];

%% Angular acceleration
% Loop 1: A-B-C-D-A
syms aBEC aCD
alpha_AB  = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD  = [0 0 aCD];

a_B_A = cross(alpha_AB,  B-A) + cross(omega_AB,  cross(omega_AB,  B-A));
a_C_B = cross(alpha_BEC, C-B) + cross(omegaBEC,  cross(omegaBEC,  C-B));
a_D_C = cross(alpha_CD,  D-C) + cross(omegaCD,   cross(omegaCD,   D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;
loop1AccSolution = solve(eqn13,[aBEC aCD]);

alphaBEC(1) = double(loop1AccSolution.aBEC);
alphaCD(1)  = double(loop1AccSolution.aCD);

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector  = [0 0 alphaCD];

% Loop 2: D-C-E-F-G-D
syms aEF aFG
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

a_C_D = cross(alphaCD_vector,  C-D) + cross(omegaCD,  cross(omegaCD,  C-D));
a_E_C = cross(alphaBEC_vector, E-C) + cross(omegaBEC, cross(omegaBEC, E-C));
a_F_E = cross(alpha_EF, F-E)        + cross(omegaEF,  cross(omegaEF,  F-E));
a_G_F = cross(alpha_FG, G-F)        + cross(omegaFG,  cross(omegaFG,  G-F));

eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;
loop2AccSolution = solve(eqn14,[aEF aFG]);


alphaEF(1) = double(loop2AccSolution.aEF);
alphaFG(1) = double(loop2AccSolution.aFG);

%% Velocity of joints B, C, E, F (A, D, G are fixed pivots -> zero velocity)
VA = [0 0 0];
VB = cross(omega_AB, B-A);
VC = VB + cross(omegaBEC, C-B);
VE = VB + cross(omegaBEC, E-B);
VF = VE + cross(omegaEF,  F-E);

%these should evaluate to ~[0 0 0] since D and G are fixed
VD_check = VC + cross(omegaCD, D-C);
VG_check = VF + cross(omegaFG, G-F);

%% Acceleration of joints B, C, E, F (A, D, G are fixed pivots -> zero acceleration)
AA = [0 0 0];
AB = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
AC = AB + cross(alphaBEC_vector, C-B) + cross(omegaBEC, cross(omegaBEC, C-B));
AE = AB + cross(alphaBEC_vector, E-B) + cross(omegaBEC, cross(omegaBEC, E-B));

alphaEF_vector = [0 0 alphaEF];
alphaFG_vector = [0 0 alphaFG];

AF = AE + cross(alphaEF_vector, F-E) + cross(omegaEF, cross(omegaEF, F-E));

% these should evaluate to ~[0 0 0] since D and G are fixed
AD_check = AC + cross(alphaCD_vector, D-C) + cross(omegaCD, cross(omegaCD, D-C));
AG_check = AF + cross(alphaFG_vector, G-F) + cross(omegaFG, cross(omegaFG, G-F));

%% Velocity of each link's center of mass (Vs1 - Vs5)
% Link 1 (AB) pivots about fixed A
Vs1 = cross(omega_AB, S1-A);
% Link 2 (BEC) - reference off known point B
Vs2 = VB + cross(omegaBEC, S2-B);
% Link 3 (CD) - reference off known point C
Vs3 = VC + cross(omegaCD, S3-C);
% Link 4 (EF) - reference off known point E
Vs4 = VE + cross(omegaEF, S4-E);
% Link 5 (FG) - reference off known point F
Vs5 = VF + cross(omegaFG, S5-F);

%% Acceleration of each link's center of mass (As1 - As5)
% Link 1 (AB) pivots about fixed A
As1 = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));
% Link 2 (BEC) - reference off known point B
As2 = AB + cross(alphaBEC_vector, S2-B) + cross(omegaBEC, cross(omegaBEC, S2-B));
% Link 3 (CD) - reference off known point C
As3 = AC + cross(alphaCD_vector, S3-C) + cross(omegaCD, cross(omegaCD, S3-C));
% Link 4 (EF) - reference off known point E
As4 = AE + cross(alphaEF_vector, S4-E) + cross(omegaEF, cross(omegaEF, S4-E));
% Link 5 (FG) - reference off known point F
As5 = AF + cross(alphaFG_vector, S5-F) + cross(omegaFG, cross(omegaFG, S5-F));

%% Display results
disp('Joint velocities [VB VC VE VF]:');
disp([VB; VC; VE; VF]);
disp('Joint accelerations [AB AC AE AF]:');
disp([AB; AC; AE; AF]);

disp('Center-of-mass velocities [Vs1 Vs2 Vs3 Vs4 Vs5]:');
disp([Vs1; Vs2; Vs3; Vs4; Vs5]);
disp('Center-of-mass accelerations [As1 As2 As3 As4 As5]:');
disp([As1; As2; As3; As4; As5]);

%NEwtons 2nd law


%Mass MOI

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NTin

%Define Forces

NForceA=[NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceT = [0 0 NTin];

%Equations for LinkAB
%sum of forces
eqn15=NForceA+NForceB+WAB == MassAB*As1;

%sum of moments=0
%sum of moments=0
eqn16 = cross(A-S1,NForceA) + cross(B-S1,NForceB) + NForceT == J_AB*alpha_AB;    % Link AB

% Equations for Link BEC
% sum of forces
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBC*As2;
% sum of moments about S2
eqn18 = cross(B-S2, -NForceB) + cross(C-S2, NForceC) + cross(E-S2, NForceE) == J_BC*alpha_BEC;

% Equations for Link CD
% sum of forces
eqn19 = -NForceC + NForceD + WCD == MassDE*As3;
% sum of moments about S3
eqn20 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) == J_DE*alpha_CD;

% Equations for Link EF
% sum of forces
eqn21 = -NForceE + NForceF + NForceA + WEF == MassEF*As4;
% sum of moments about S4
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) + cross(A-S4, NForceA) == J_EF*alpha_EF;

% Equations for Link FG
% sum of forces
eqn23 = -NForceF + NForceD + WFG == MassFG*As5;
% sum of moments about S5
eqn24 = cross(F-S5, -NForceF) + cross(G-S5, NForceD) == J_FG*alpha_FG;

% Collect equations
eqnNewton = [eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
dynamicSolution=solve(eqnNewton,[NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NTin])

% Example: substitute known symbolic parameters into the solution, then convert
symList = [aBEC, aCD, aEF, aFG, wBEC, wCD, wEF, wFG, FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, FEx, FEy, FFx, FFy, FGx, FGy, Tin];
valList = {alphaBEC, alphaCD, alphaEF, alphaFG, angularVelocity_BEC, angularVelocity_CD, angularVelocity_EF, angularVelocity_FG, Force_Ax, Force_Ay, Force_Bx, Force_By, Force_Cx, Force_Cy, Force_Dx, Force_Dy, Force_Ex, Force_Ey, Force_Fx, Force_Fy, Force_Gx, Force_Gy, Input_Torque};
dynamicSolutionSubs = structfun(@(s) subs(s, symList, valList), dynamicSolution, 'UniformOutput', false);

NFAx_val(1) = double(dynamicSolutionSubs.NFAx);
NFAy_val(1) = double(dynamicSolutionSubs.NFAy);
NFBx_val(1) = double(dynamicSolutionSubs.NFBx);
NFBy_val(1) = double(dynamicSolutionSubs.NFBy);
NFCx_val(1) = double(dynamicSolutionSubs.NFCx);
NFCy_val(1) = double(dynamicSolutionSubs.NFCy);
NFDx_val(1) = double(dynamicSolutionSubs.NFDx);
NFDy_val(1) = double(dynamicSolutionSubs.NFDy);
NFEx_val(1) = double(dynamicSolutionSubs.NFEx);
NFEy_val(1) = double(dynamicSolutionSubs.NFEy);
NFFx_val(1) = double(dynamicSolutionSubs.NFFx);
NFFy_val(1) = double(dynamicSolutionSubs.NFFy);
NTin_val(1) = double(dynamicSolutionSubs.NTin);

% Display extracted forces
disp(NForceA)
% Display extracted numeric reaction forces and torque
disp('Reaction forces at A [NFAx NFAy]:');
disp([NFAx_val NFAy_val]);
disp('Reaction forces at B [NFBx NFBy]:');
disp([NFBx_val NFBy_val]);
disp('Reaction forces at C [NFCx NFCy]:');
disp([NFCx_val NFCy_val]);
disp('Reaction forces at D [NFDx NFDy]:');
disp([NFDx_val NFDy_val]);
disp('Reaction forces at E [NFEx NFEy]:');
disp([NFEx_val NFEy_val]);
disp('Reaction forces at F [NFFx NFFy]:');
disp([NFFx_val NFFy_val]);
disp('Input torque (NTin):');
disp(NTin_val);

%Circle Intersection Technique

%Joint coordinates defined
%Length of links also defined
%% 

% Compute the initial angle of the input link
Initial_Theta = atan2(B(2)-A(2),B(1)-A(1));
disp('Initial angle of input link (Initial_Theta): ');
disp(Initial_Theta)

if(Initial_Theta<0)
    inputAngle=2*pi+Initial_Theta;
else
    inputAngle=Initial_Theta;
end
new_B_x(1) = B(1);  new_B_y(1) = B(2);
new_C_x(1) = C(1);  new_C_y(1) = C(2);
new_E_x(1) = E(1);  new_E_y(1) = E(2);
new_F_x(1) = F(1);  new_F_y(1) = F(2);
for theta=1:1:360

    %new Position of joint B
    B_new=A + [lAB*cos(inputAngle+deg2rad(theta)) lAB*sin(inputAngle+deg2rad(theta)) 0];
    %new Position of C
    [Cx,Cy]=circcirc(B_new(1),B_new(2),lBC,D(1),D(2),lCD);
    
    %Checking if there is a NAN

    circIntersect_x=any(isnan(vpa(Cx)));
    circIntersect_y=any(isnan(vpa(Cy)));

    if circIntersect_x==0 && circIntersect_y==0
        C_1=[Cx(1) Cy(1) 0];
        C_2 = [Cx(2) Cy(2) 0];

        dist1=norm(C_1-C);
        dist2=norm(C_2-C);

        if(dist1<dist2)
            C_new=vpa(C_1);
        else

            C_new=vpa(C_2);
        end

        %New Position of Joint E using B_new and C_new
        
        %new Position of Joint E
        [Ex,Ey]=circcirc(B_new(1),B_new(2),lBE,C_new(1),C_new(2),lCE);
    
        % Check for NaNs in intersection
        circIntersect_ex=any(isnan(vpa(Ex)));
        circIntersect_ey=any(isnan(vpa(Ey)));
    
        if circIntersect_ex==0 && circIntersect_ey==0
        E_1 = [Ex(1) Ey(1) 0];
        E_2 = [Ex(2) Ey(2) 0];
        
        distE1 = norm(E_1 - E);
        distE2 = norm(E_2 - E);
        
        if(distE1 < distE2)
            E_new = vpa(E_1);
        else
            E_new = vpa(E_2);
        end
        
        
        % New positions of joints F and G using E_new and fixed point G
        % Joint F from E_new and G via lengths lEF and lEG (use E_new and G(1),G(2))
        [Fx,Fy]=circcirc(E_new(1),E_new(2),lEF,G(1),G(2),lEG);
        
        %Checking if there is a NAN
        circIntersect_fx=any(isnan(vpa(Fx)));
        circIntersect_fy=any(isnan(vpa(Fy)));

        if circIntersect_fx==0 && circIntersect_fy==0
            F_1=[Fx(1) Fy(1) 0];
            F_2=[Fx(2) Fy(2) 0];

            distF1=norm(F_1-F);
            distF2=norm(F_2-F);

            if(distF1<distF2)
                F_new=vpa(F_1);
            else
                F_new=vpa(F_2);
            end
            
          %Store values for plotting
            new_B_x(theta+1)=B_new(1);
            new_B_y(theta+1)=B_new(2);
            new_C_x(theta+1)=C_new(1);
            new_C_y(theta+1)=C_new(2);
            new_E_x(theta+1)=E_new(1);
            new_E_y(theta+1)=E_new(2);
            new_F_x(theta+1)=F_new(1);
            new_F_y(theta+1)=F_new(2);
            
            % Use THIS iteration's swept positions everywhere below.
            % A, D, G are fixed ground pivots and never move.
            Bc = double(B_new);
            Cc = double(C_new);
            Ec = double(E_new);
            Fc = double(F_new);
            
            % center of mass of each link (recomputed for this position)
            S1 = (A+Bc)/2;
            S2 = (Bc+Cc+Ec)/3;
            S3 = (Cc+D)/2;
            S4 = (Ec+Fc)/2;
            S5 = (Fc+G)/2;
            
            syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin
            
            ForceA = [FAx FAy 0];
            ForceB = [FBx FBy 0];
            ForceC = [FCx FCy 0];
            ForceD = [FDx FDy 0];
            ForceE = [FEx FEy 0];
            ForceF = [FFx FFy 0];
            ForceG = [FGx FGy 0];
            InputTorque = [0 0 Tin];
            
            AppliedForce = [0 200 0];
            
            eqn1 = ForceA + ForceB + WAB == 0;
            eqn2 = cross(A-S1,ForceA) + cross(Bc-S1,ForceB) + InputTorque == 0;
            
            eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;
            eqn4 = cross(Bc - S2, -ForceB) + cross(Cc - S2, ForceC) + cross(Ec - S2, ForceE) == 0;
            
            eqn5 = -ForceC + ForceD + WCD == 0;
            eqn6 = cross(Cc - S3, -ForceC) + cross(D - S3, ForceD) == 0;
            
            eqn7 = -ForceE + ForceF + WEF == 0;
            eqn8 = cross(Ec - S4, -ForceE) + cross(Fc - S4, ForceF) == 0;
            
            eqn9  = -ForceF + ForceG + WFG + AppliedForce == 0;
            eqn10 = cross(Fc - S5, -ForceF) + cross(G - S5, ForceG) == 0;
            
            eqnMatrix = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
            StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);
            
            Force_Ax(theta+1) = double(StaticSolution.FAx);
            Force_Ay(theta+1) = double(StaticSolution.FAy);
            Force_Bx(theta+1) = double(StaticSolution.FBx);
            Force_By(theta+1) = double(StaticSolution.FBy);
            Force_Cx(theta+1) = double(StaticSolution.FCx);
            Force_Cy(theta+1) = double(StaticSolution.FCy);
            Force_Dx(theta+1) = double(StaticSolution.FDx);
            Force_Dy(theta+1) = double(StaticSolution.FDy);
            Force_Ex(theta+1) = double(StaticSolution.FEx);
            Force_Ey(theta+1) = double(StaticSolution.FEy);
            Force_Fx(theta+1) = double(StaticSolution.FFx);
            Force_Fy(theta+1) = double(StaticSolution.FFy);
            Force_Gx(theta+1) = double(StaticSolution.FGx);
            Force_Gy(theta+1) = double(StaticSolution.FGy);
            Input_Torque(theta+1) = double(StaticSolution.Tin);
            
            fprintf('Static equilibrium forces (N) at theta = %d degrees:\n', theta);
            fprintf('Force A: [%.4f, %.4f]\n', Force_Ax(theta+1), Force_Ay(theta+1));
            fprintf('Force B: [%.4f, %.4f]\n', Force_Bx(theta+1), Force_By(theta+1));
            fprintf('Force C: [%.4f, %.4f]\n', Force_Cx(theta+1), Force_Cy(theta+1));
            fprintf('Force D: [%.4f, %.4f]\n', Force_Dx(theta+1), Force_Dy(theta+1));
            fprintf('Force E: [%.4f, %.4f]\n', Force_Ex(theta+1), Force_Ey(theta+1));
            fprintf('Force F: [%.4f, %.4f]\n', Force_Fx(theta+1), Force_Fy(theta+1));
            fprintf('Force G: [%.4f, %.4f]\n', Force_Gx(theta+1), Force_Gy(theta+1));
            fprintf('Input Torque (Nm): %.4f\n', Input_Torque(theta+1));
            
            %% Angular velocity (recomputed for THIS position)
            syms wBEC wCD
            omega_AB  = [0 0 23.14];
            omega_BEC = [0 0 wBEC];
            omega_CD  = [0 0 wCD];
            
            eqn11 = cross(omega_AB,Bc-A) + cross(omega_BEC,Cc-Bc) + cross(omega_CD,D-Cc) == 0;
            loop1Solution = solve(eqn11,[wBEC wCD]);
            
            angularVelocity_BEC(theta+1) = double(loop1Solution.wBEC);
            angularVelocity_CD(theta+1)  = double(loop1Solution.wCD);
            
            omegaBEC = [0 0 angularVelocity_BEC(theta+1)];
            omegaCD  = [0 0 angularVelocity_CD(theta+1)];
            
            syms wEF wFG
            omega_EF = [0 0 wEF];
            omega_FG = [0 0 wFG];
            
            eqn12 = cross(omegaCD,D-Cc) + cross(omegaBEC,Ec-Cc) + cross(omega_EF,Fc-Ec) + cross(omega_FG,G-Fc) == 0;
            loop2solution = solve(eqn12,[wEF wFG]);
            
            angularVelocity_EF(theta+1) = double(loop2solution.wEF);
            angularVelocity_FG(theta+1) = double(loop2solution.wFG);
            
            omegaEF = [0 0 angularVelocity_EF(theta+1)];
            omegaFG = [0 0 angularVelocity_FG(theta+1)];
            
            %% Angular acceleration
            syms aBEC aCD
            alpha_AB  = [0 0 0];
            alpha_BEC = [0 0 aBEC];
            alpha_CD  = [0 0 aCD];
            
            a_B_A = cross(alpha_AB,  Bc-A)  + cross(omega_AB,  cross(omega_AB,  Bc-A));
            a_C_B = cross(alpha_BEC, Cc-Bc) + cross(omegaBEC,  cross(omegaBEC,  Cc-Bc));
            a_D_C = cross(alpha_CD,  D-Cc)  + cross(omegaCD,   cross(omegaCD,   D-Cc));
            
            eqn13 = a_B_A + a_C_B + a_D_C == 0;
            loop1AccSolution = solve(eqn13,[aBEC aCD]);
            
            alphaBEC(theta+1) = double(loop1AccSolution.aBEC);
            alphaCD(theta+1)  = double(loop1AccSolution.aCD);
            
            alphaBEC_vector(theta+1,:) = [0 0 alphaBEC(theta+1)];
            alphaCD_vector(theta+1,:)  = [0 0 alphaCD(theta+1)];
            
            syms aEF aFG
            alpha_EF = [0 0 aEF];
            alpha_FG = [0 0 aFG];
            
            a_C_D = cross(alphaCD_vector(theta+1,:),  Cc-D)  + cross(omegaCD,  cross(omegaCD,  Cc-D));
            a_E_C = cross(alphaBEC_vector(theta+1,:), Ec-Cc) + cross(omegaBEC, cross(omegaBEC, Ec-Cc));
            a_F_E = cross(alpha_EF, Fc-Ec) + cross(omegaEF, cross(omegaEF, Fc-Ec));
            a_G_F = cross(alpha_FG, G-Fc)  + cross(omegaFG, cross(omegaFG, G-Fc));
            
            eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;
            loop2AccSolution = solve(eqn14,[aEF aFG]);
            
            alphaEF(theta+1) = double(loop2AccSolution.aEF);
            alphaFG(theta+1) = double(loop2AccSolution.aFG);
            
            alphaEF_vector = [0 0 alphaEF(theta+1)];
            alphaFG_vector = [0 0 alphaFG(theta+1)];
            
            %% Joint velocities
            VB = cross(omega_AB, Bc-A);
            VC = VB + cross(omegaBEC, Cc-Bc);
            VE = VB + cross(omegaBEC, Ec-Bc);
            VF = VE + cross(omegaEF,  Fc-Ec);
            
            VD_check = VC + cross(omegaCD, D-Cc);   %#ok<NASGU> % should be ~0
            VG_check = VF + cross(omegaFG, G-Fc);   %#ok<NASGU> % should be ~0
            
            %% Joint accelerations
            AB = cross(alpha_AB, Bc-A) + cross(omega_AB, cross(omega_AB, Bc-A));
            AC = AB + cross(alphaBEC_vector(theta+1,:), Cc-Bc) + cross(omegaBEC, cross(omegaBEC, Cc-Bc));
            AE = AB + cross(alphaBEC_vector(theta+1,:), Ec-Bc) + cross(omegaBEC, cross(omegaBEC, Ec-Bc));
            AF = AE + cross(alphaEF_vector, Fc-Ec) + cross(omegaEF, cross(omegaEF, Fc-Ec));
            
            AD_check = AC + cross(alphaCD_vector(theta+1,:), D-Cc) + cross(omegaCD, cross(omegaCD, D-Cc)); %#ok<NASGU>
            AG_check = AF + cross(alphaFG_vector, G-Fc) + cross(omegaFG, cross(omegaFG, G-Fc)); %#ok<NASGU>
            
            %% Link CM velocities/accelerations
            Vs1 = cross(omega_AB, S1-A);
            Vs2 = VB + cross(omegaBEC, S2-Bc);
            Vs3 = VC + cross(omegaCD, S3-Cc);
            Vs4 = VE + cross(omegaEF, S4-Ec);
            Vs5 = VF + cross(omegaFG, S5-Fc);
            
            As1 = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));
            As2 = AB + cross(alphaBEC_vector(theta+1,:), S2-Bc) + cross(omegaBEC, cross(omegaBEC, S2-Bc));
            As3 = AC + cross(alphaCD_vector(theta+1,:), S3-Cc) + cross(omegaCD, cross(omegaCD, S3-Cc));
            As4 = AE + cross(alphaEF_vector, S4-Ec) + cross(omegaEF, cross(omegaEF, S4-Ec));
            As5 = AF + cross(alphaFG_vector, S5-Fc) + cross(omegaFG, cross(omegaFG, S5-Fc));
            
            fprintf('Kinematics at theta = %d degrees:\n', theta);
            disp('Joint velocities [VB VC VE VF]:');   disp([VB; VC; VE; VF]);
            disp('Joint accelerations [AB AC AE AF]:'); disp([AB; AC; AE; AF]);
            disp('CM velocities [Vs1..Vs5]:');          disp([Vs1; Vs2; Vs3; Vs4; Vs5]);
            disp('CM accelerations [As1..As5]:');       disp([As1; As2; As3; As4; As5]);
            
            %% Newton's second law (same pattern as statics above, this theta's position)
            syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin
            
            NForceA = [NFAx NFAy 0];
            NForceB = [NFBx NFBy 0];
            NForceC = [NFCx NFCy 0];
            NForceD = [NFDx NFDy 0];
            NForceE = [NFEx NFEy 0];
            NForceF = [NFFx NFFy 0];
            NForceG = [NFGx NFGy 0];
            NForceT = [0 0 NTin];
            
            % Link AB
            eqn15 = NForceA + NForceB + WAB == MassAB*As1;
            eqn16 = cross(A-S1,NForceA) + cross(Bc-S1,NForceB) + NForceT == J_AB*alpha_AB;
            
            % Link BEC
            eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBC*As2;
            eqn18 = cross(Bc-S2, -NForceB) + cross(Cc-S2, NForceC) + cross(Ec-S2, NForceE) == J_BC*alphaBEC_vector(theta+1,:);
            
            % Link CD
            eqn19 = -NForceC + NForceD + WCD == MassDE*As3;
            eqn20 = cross(Cc-S3, -NForceC) + cross(D-S3, NForceD) == J_DE*alphaCD_vector(theta+1,:);
            
            % Link EF (fixed: only touches E and F)
            eqn21 = -NForceE + NForceF + WEF == MassEF*As4;
            eqn22 = cross(Ec-S4, -NForceE) + cross(Fc-S4, NForceF) == J_EF*alphaEF_vector;
            
            % Link FG (fixed: ground reaction at G is NForceG, not NForceD)
            eqn23 = -NForceF + NForceG + WFG == MassFG*As5;
            eqn24 = cross(Fc-S5, -NForceF) + cross(G-S5, NForceG) == J_FG*alphaFG_vector;
            
            eqnNewton = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
            dynamicSolution = solve(eqnNewton, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);
            
            NFAx_val(theta+1) = double(dynamicSolution.NFAx);
            NFAy_val(theta+1) = double(dynamicSolution.NFAy);
            NFBx_val(theta+1) = double(dynamicSolution.NFBx);
            NFBy_val(theta+1) = double(dynamicSolution.NFBy);
            NFCx_val(theta+1) = double(dynamicSolution.NFCx);
            NFCy_val(theta+1) = double(dynamicSolution.NFCy);
            NFDx_val(theta+1) = double(dynamicSolution.NFDx);
            NFDy_val(theta+1) = double(dynamicSolution.NFDy);
            NFEx_val(theta+1) = double(dynamicSolution.NFEx);
            NFEy_val(theta+1) = double(dynamicSolution.NFEy);
            NFFx_val(theta+1) = double(dynamicSolution.NFFx);
            NFFy_val(theta+1) = double(dynamicSolution.NFFy);
            NFGx_val(theta+1) = double(dynamicSolution.NFGx);
            NFGy_val(theta+1) = double(dynamicSolution.NFGy);
            NTin_val(theta+1)  = double(dynamicSolution.NTin);
            
            fprintf('Dynamic reaction forces (N) at theta = %d degrees:\n', theta);
            fprintf('Force A: [%.4f, %.4f]\n', NFAx_val(theta+1), NFAy_val(theta+1));
            fprintf('Force B: [%.4f, %.4f]\n', NFBx_val(theta+1), NFBy_val(theta+1));
            fprintf('Force C: [%.4f, %.4f]\n', NFCx_val(theta+1), NFCy_val(theta+1));
            fprintf('Force D: [%.4f, %.4f]\n', NFDx_val(theta+1), NFDy_val(theta+1));
            fprintf('Force E: [%.4f, %.4f]\n', NFEx_val(theta+1), NFEy_val(theta+1));
            fprintf('Force F: [%.4f, %.4f]\n', NFFx_val(theta+1), NFFy_val(theta+1));
            fprintf('Force G: [%.4f, %.4f]\n', NFGx_val(theta+1), NFGy_val(theta+1));
            fprintf('Input Torque (Nm): %.4f\n', NTin_val(theta+1));

                
        else
            fprintf('new position of F cannot be determined at angle %d degree\n',theta);
        end
    else
        fprintf('new position of E cannot be determined at angle %d degree',theta);
    end
        
    else
        fprintf('new position of C cannot be determined at angle %d degree',theta);
    end

end
%% Plots across the full 360-degree input rotation
theta_deg = 0:360;   % matches indices 1:361

figure('Name','Joint Positions');
plot(new_B_x, new_B_y, new_C_x, new_C_y, new_E_x, new_E_y, new_F_x, new_F_y, 'LineWidth', 1.5);
legend('B','C','E','F'); xlabel('X (m)'); ylabel('Y (m)');
title('Joint Position Trajectories'); axis equal; grid on;

figure('Name','Angular Velocities');
plot(theta_deg, angularVelocity_BEC, theta_deg, angularVelocity_CD, ...
    theta_deg, angularVelocity_EF,  theta_deg, angularVelocity_FG, 'LineWidth', 1.5);
legend('\omega_{BEC}','\omega_{CD}','\omega_{EF}','\omega_{FG}');
xlabel('Input angle \theta (deg)'); ylabel('Angular velocity (rad/s)');
title('Link Angular Velocities vs Input Angle'); grid on;

figure('Name','Angular Accelerations');
plot(theta_deg, alphaBEC, theta_deg, alphaCD, theta_deg, alphaEF, theta_deg, alphaFG, 'LineWidth', 1.5);
legend('\alpha_{BEC}','\alpha_{CD}','\alpha_{EF}','\alpha_{FG}');
xlabel('Input angle \theta (deg)'); ylabel('Angular acceleration (rad/s^2)');
title('Link Angular Accelerations vs Input Angle'); grid on;

figure('Name','Static Reaction Forces');
plot(theta_deg, hypot(Force_Bx,Force_By), theta_deg, hypot(Force_Cx,Force_Cy), ...
    theta_deg, hypot(Force_Ex,Force_Ey), theta_deg, hypot(Force_Fx,Force_Fy), ...
    theta_deg, hypot(Force_Gx,Force_Gy), 'LineWidth', 1.5);
legend('|F_B|','|F_C|','|F_E|','|F_F|','|F_G|');
xlabel('Input angle \theta (deg)'); ylabel('Force magnitude (N)');
title('Static Joint Reaction Forces vs Input Angle'); grid on;

figure('Name','Static Input Torque');
plot(theta_deg, Input_Torque, 'LineWidth', 1.5);
xlabel('Input angle \theta (deg)'); ylabel('Torque (N\cdotm)');
title('Static Input Torque vs Input Angle'); grid on;

figure('Name','Dynamic Reaction Forces');
plot(theta_deg, hypot(NFBx_val,NFBy_val), theta_deg, hypot(NFCx_val,NFCy_val), ...
    theta_deg, hypot(NFEx_val,NFEy_val), theta_deg, hypot(NFFx_val,NFFy_val), ...
    theta_deg, hypot(NFGx_val,NFGy_val), 'LineWidth', 1.5);
legend('|F_B|','|F_C|','|F_E|','|F_F|','|F_G|');
xlabel('Input angle \theta (deg)'); ylabel('Force magnitude (N)');
title('Dynamic Joint Reaction Forces vs Input Angle'); grid on;

figure('Name','Dynamic Input Torque');
plot(theta_deg, NTin_val, 'LineWidth', 1.5);
xlabel('Input angle \theta (deg)'); ylabel('Torque (N\cdotm)');
title('Dynamic Input Torque vs Input Angle'); grid on;

