%% EEE 588: Analysis and Design of a GE-21 Engine

%Unscaled State Space Dynamics

%Operating Point 1 (OP1)

Ap1=[-3.4272,1.7566;-0.3110,-1.9281];
Bp1=[0.4909,1.1449,0.1894;0.4148,0.1405,-0.2414];
Cp1=[1 0; 0 1; 0.5211 2.6808];
Dp1=[0,0,0;0,0,0;0.5954,-0.6348,1.0206];

%Operating Point 2 (OP2)

Ap2=[-3.488, 1.169; -0.597, -1.847];
Bp2=[0.0655,4.3839; 0.0755,0.9098];
Cp2=[1,0;0,1];
Dp2=[0,0;0,0];

%Natural Modes: Poles (Eigenvalues), Eigenvectors

[evecR1,eval1,evecL1] = eig(Ap1) %evecR1 contains eigenvectors of system operating at point 1
                         %eval1 contains poles or eigenvalues of system
                         %operating at point 1
                      

[evecR2,eval2,evecL2] = eig(Ap2) %evecR2 contains eigenvectors of system operating at point 2
                         %eval2 contains poles or eigenvalues of system
                         %operating at point 2

%Transmission Zeros

z1 = tzero(ss(Ap1,Bp1,Cp1,Dp1))                % transmission zeros (No Transmission zeros)
%zdir = null([z.*eye(2)-Ap  -Bp; Cp Dp])   % transmission zero directions

z2 = tzero(ss(Ap2,Bp2,Cp2,Dp2))                % transmission zeros (No Transmission zeros)
%zdir = null([z.*eye(2)-Ap  -Bp; Cp Dp])   % transmission zero directions

%TRANSFER FUNCTIONS

sys1 = zpk(ss(Ap1,Bp1,Cp1,Dp1)) % Zeros, Poles, and Gains fron u_i to x_j for Plant at OP1
bode(sys1)
title('Frequence Responses of Transfer Functions ')
grid
xlabel('Frequency (rad/sec)')
ylabel('Bode Magnitude')
pause

sys2 = zpk(ss(Ap2,Bp2,Cp2,Dp2)) % Zeros, Poles, and Gains fron u_i to x_j for Plant at OP2
bode(sys2)
title('Frequence Responses of Transfer Functions ')
grid
xlabel('Frequency (rad/sec)')
ylabel('Bode Magnitude')
pause 

% Controllability 
cm1 = [Bp1,Ap1*Bp1]  % Controllability Matrix of Plant at OP1
rcm1= rank(cm1)              % Rank of Controllability Matrix

cm2 = [Bp2,Ap2*Bp2]  % Controllability Matrix of Plant at OP2
rcm2= rank(cm2)              % Rank of Controllability Matrix
% Observability

om1 = [Cp1; Cp1*Ap1;]          % Observability Matrix of Plant at OP1
rom1 = rank(om1)             % Rank of Observability Matrix

om2 = [Cp2; Cp2*Ap2;]          % Observability Matrix of Plant at OP2
rom2 = rank(om2)             % Rank of Observability Matrix

%% 
%Plant Frequency Response

sigma(ss(Ap1, Bp1, Cp1, Dp1));

title('Frequence Response of Plant at OP1 ')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(ss(Ap2, Bp2, Cp2, Dp2));

title('Frequency Response of Plant at OP2')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%***************************************************************************
%% SVD Analysis

% SVD Analysis at DC

DC1 =  Cp1*inv(-Ap1)*Bp1 + Dp1
[udc1,sdc1,vdc1] = svd(DC1)

DC2 =  Cp2*inv(-Ap2)*Bp2 + Dp2
[udc2,sdc2,vdc2] = svd(DC2)


%% Plant Augmentation with bank of Integrators
%Augment Plant with Integrators at Plant Input and Plot Singular Values

%Operating Point 1
[ns1 nc1] = size(Bp1);                      % ns = number of inputs(states);  nc = number of controls (u1,u2,u3);   
A1 = [ Ap1             Bp1
      0*ones(nc1,ns1)    0*ones(nc1,nc1) ];
B1 = [ 0*ones(ns1,nc1)
      eye(nc1)      ]; 
C1 = [ Cp1 , Dp1 ];
D1 = 0*ones(nc1,nc1);

sigma(ss(A1, B1, C1, D1));

title('Design Plant Singular Values for OP1')
grid on
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%**************************************************************************************
%Operating Point 2

[ns2 nc2] = size(Bp2);                      % ns = number of inputs(states);  nc = number of controls (u1,u2,u3);   
A2 = [ Ap2             Bp2
      0*ones(nc2,ns2)    0*ones(nc2,nc2) ];
B2 = [ 0*ones(ns2,nc2)
      eye(nc2)      ];
C2 = [ Cp2 , Dp2 ];
D2 = 0*ones(nc2,nc2);

sigma(ss(A2, B2, C2, D2));

title('Design Plant Singular Values for OP2')
grid on
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

%% LQR: Control Design

%PART 1: LQR Design for Operating Point 1

Q1 = diag([1, 1, 1, 0, 0]); 
rho1 = 1e-2 ;                                         
R1 = rho1*eye(nc1)


[G1, K1, clpoles1] = lqr(A1, B1, Q1, R1)


% LQ OPEN LOOP FREQUENCY RESPONSE 
%
w = logspace(-3,3,100);
sv = sigma(ss(A1, B1, G1, 0*ones(3,3)),w);
sv = 20*log10(sv); 
semilogx(w, sv)
%clear sv
title('Open Loop Singular Values: Plant Input for OP1 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%LQ CLOSED LOOP FREQUENCY RESPONSE 

sv = sigma(ss(A1-B1*G1, B1, -G1, eye(3,3)-0*ones(3,3)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('LQ  Sensitivity: Plant Input for OP1 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
t=[0:1:100];
step(ss(A1-B1*G1,B1,-G1,eye(3,3)-0*ones(3,3)),t); %Step Response of S at Input
title('Step Response from r to e at Input for OP1 (LQR)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

sv = sigma(ss(A1-B1*G1, B1, G1, 0*ones(3,3)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('LQ Complementary Sensitivity: Plant Input for OP1 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
t=[0:1:100];
step(ss(A1-B1*G1,B1,G1,0*ones(3,3)),t); %Step Response of T at Input
title('Step Response from r to y for OP1 (LQR)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

%%
%PART 2: LQR Design for Operating Point 2

Q2 = diag([1, 1, 1, 0]); 
rho2 = 1 ;                                         
R2 = rho2*eye(nc2)


[G2, K2, clpoles2] = lqr(A2, B2, Q2, R2)

% LQ OPEN LOOP FREQUENCY RESPONSE 

w = logspace(-3,3,100);
sv = sigma(ss(A2, B2, G2, 0*ones(2,2)),w);
sv = 20*log10(sv); 
semilogx(w, sv)

title('Open Loop Singular Values: Plant Input for OP2 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%LQ CLOSED LOOP FREQUENCY RESPONSE 

sv = sigma(ss(A2-B2*G2, B2, -G2, eye(2,2)-0*ones(2,2)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('LQ  Sensitivity: Plant Input for OP2 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
t=[0:1:100];
step(ss(A2-B2*G2,B2,-G2,eye(2,2)-0*ones(2,2)),t); %Step Response of S at Input
title('Step Response from r to e at Input for OP2 (LQR)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

sv = sigma(ss(A2-B2*G2, B2, G2, 0*ones(2,2)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('LQ Complementary Sensitivity: Plant Input for OP2 (LQR)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
t=[0:1:100];
step(ss(A2-B2*G2,B2,G2,0*ones(2,2)),t); %Step Response of T at Input
title('Step Response from r to y at Input for OP2 (LQR)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

%% H-Infinity: Control Design

%PART 1: Operating Point 1

s = zpk('s');
P1=ss(Ap1,Bp1,Cp1,Dp1);

W1_1 = ((.1)*(s+1)/(s+0.00001))*(eye(3))
W2_1 = ((.1)*((s+.1)/(s+1000)))*(eye(3))
W3_1 = ((.1)*(s+1))/(s+10000)*(eye(3))

[K1,CL1,GAM1] = mixsyn(P1,W1_1,W2_1,W3_1);

S1 = feedback(eye(size(P1*K1)),P1*K1);
R1 = K1*S1;
T1 = eye(size(S1))-S1;
 
w=logspace(-4,2,100);

%**************************************************************************

% Form closed loop maps
% f_CLTFM.m is a matlab function for computing OL and CL maps in a standard
% output feedback structure. This file must be in the current Matlab folder

[Lo1,Li1,So1,Si1,To1,Ti1,KS1,PS1] = f_CLTFM(P1,K1);

% PLOTS
wvec=logspace(-3,3,10000);
 
%**************************************************************************
%
% OPEN LOOP FREQUENCY RESPONSE 

sigma(Lo1,wvec);
title('Open Loop Singular Values: Error Signal for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(Li1,wvec);
title('Open Loop Singular Values: Plant Input for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%**************************************************************************
%
% CLOSED LOOP FREQUENCY RESPONSE 

sigma(So1,wvec);
title('Sensitivity: Error Signal for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
 
sigma(Si1,wvec);
title('Sensitivity: Plant Input for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(R1,w); legend('KS'); 
title('KS (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(To1,wvec);
title('Complementary Sensitivity: Plant Output for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(Ti1,wvec);
title('Complementary Sensitivity: Plant Input for OP1 (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

t=[0:1:100];
step(To1,t); %Step Response 
title('Step Response from r to y for OP1 (H-Infinity)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

%% H-Infinity: Control Design

%PART 2: Operating Point 2

s = zpk('s');
P2=ss(Ap2,Bp2,Cp2,Dp2);

W1_2 = ((0.1)*(s+1)/(s+0.000001))*(eye(2));
W2_2 = ((0.1)*((s+0.1)/(s+ 100))*(eye(2)));
W3_2 = ((0.1)*(s+1))/(s+10000)*(eye(2));

[K2,CL2,GAM2] = mixsyn(P2,W1_2,W2_2,W3_2)

S2 = feedback(eye(size(P2*K2)),P2*K2);
R2 = K2*S2;
T2 = eye(size(S2))-S2;

w=logspace(-4,2,100);

%**************************************************************************

% Form closed loop maps
% f_CLTFM.m is a matlab function for computing OL and CL maps in a standard
% output feedback structure. This file must be in the current Matlab folder

[Lo2,Li2,So2,Si2,To2,Ti2,KS2,PS2] = f_CLTFM(P2,K2);

% PLOTS
wvec=logspace(-3,3,10000);
 
%**************************************************************************
%
% OPEN LOOP FREQUENCY RESPONSE 

sigma(Lo2,wvec);
title('Open Loop Singular Values: Error Signal (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(Li2,wvec);
title('Open Loop Singular Values: Plant Input (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%**************************************************************************
%
% CLOSED LOOP FREQUENCY RESPONSE 

sigma(So2,wvec);
title('Sensitivity: Error Signal (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(Si2,wvec);
title('Sensitivity: Plant Input (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(R2,w); legend('KS   '); 
title('KS (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(To2,wvec);
title('Complementary Sensitivity: Plant Output (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

sigma(Ti2,wvec);
title('Complementary Sensitivity: Plant Input (H-Infinity)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

t=[0:1:100];
step(To2,t); %Step Response 
title('Step Response from r to y for OP2 (H-Infinity)')
grid 
xlabel('Time')
ylabel('Amplitude')
pause

%% LQG: Control Design (LQG-LTRO)

%Operating Point 2

% Use Kalman Filter to Design Target Loop (At Output) 

% DESIGN #1
%ll = inv(Cp2*inv(-Ap2)*Bp2 +Dp2);   % Match the singular values of Gfol at Low Freq
%lh = -inv(Ap2)*Bp2*ll;              % Match the singular values of Gfol High Freq
%l  = [ll;lh];    

% DESIGN #2
l2=B2;

w      = logspace(-5,5,200);            % Form vector of logarithmically spaced freq points
gfolsv = sigma(ss(A2, l2, C2, D2),w);
gfolsv = 20*log10(gfolsv);
semilogx(w, gfolsv)
%clear sv
title('G_{FOL} Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

pnint2 = eye(nc2)                                 % Process Noise Intensity Matrix
mu2 = .1;                                         % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                  % Small mu - expensive sensor   - large bandwidth
                                                  % Large mu - inexpensive sensor - small bandwidth
mnint2 = mu2*eye(nc2)                             % Measurement Noise Intensity Matrix 

[kest2, H2, sig2]= kalman(ss(A2, [B2 l2], C2, [D2 0*ones(nc2,nc2)]),pnint2, mnint2);  % Compute Filter Gain Matrix H
                         
sv = sigma(ss(A2, H2, C2, D2),w);
tsv = 20*log10(sv);
semilogx(w, tsv)
title('Target Loop (G_{KF}) Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

tolpoles = eig(A2)                              % Target Open Loop Poles
targzeros = tzero(A2,H2,C2,0*ones(nc2,nc2))     % Target Open Loop Zeros
tclpoles = eig(A2-H2*C2)                        % Target Closed Loop Poles

sv = sigma(ss(A2-H2*C2, H2, -C2, eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values for OP2 (LQG)')
grid on
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


sv = sigma(ss(A2-H2*C2, H2, C2, 0*eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv) 
%clear sv
title('Target Complementary (T_{KF}) Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% Recover Target Loop By Solving Cheap LQR Problem

Q_LQG2 = diag([1, 1, 1, 1]) ;                             % State Weighting Matrix
rho_lqg2 = 1e-9;                                          % Cheap control recovery parameter;
                                                          % The smaller the parameter, the better the recovery.
r = rho_lqg2*eye(nc2)                                     % Control Weigthing Matrix
[k2, poles, G_lqg2, rr] = care(A2,B2,Q_LQG2,r);           % Compute Control Gain Matrix G



% Compensator Analysis

Ak_lqg2 = [ A2-B2*G_lqg2-H2*C2 ,0*ones(ns2+nc2,nc2); G_lqg2,0*ones(nc2,nc2)]

Bk_lqg2 = [H2;0*ones(nc2,nc2)]

Ck_lqg2 = [0*ones(nc2, ns2+nc2), eye(nc2,nc2) ]

cpoles = eig(Ak_lqg2)                                           % Compensator Poles
czeros = tzero(A2, H2, G2, 0*ones(nc2,nc2))                     % Compensator Zeros
zerocheck = tzero(Ak_lqg2, Bk_lqg2, Ck_lqg2, 0*ones(nc2,nc2))   % Check Compensator Zeros

sv = sigma(ss(Ak_lqg2, Bk_lqg2, Ck_lqg2, 0*eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv)

title('Compensator Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% Open Loop Analysis

Aol2 = [ Ap2,Bp2*Ck_lqg2; 0*ones(ns2+nc2+nc2,ns2), Ak_lqg2]

Bol2 = [0*ones(ns2,nc2); Bk_lqg2] 
      
Col2 = [ Cp2, 0*ones(nc2,ns2+nc2+nc2) ]
    
olpoles = eig(Aol2)                                % Open Loop Poles
olzeros = tzero(Aol2,Bol2,Col2,0*ones(nc2,nc2))    % Open Loop Zeros
    
sv = sigma(ss(Aol2, Bol2, Col2, 0*eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, tsv)
%clear sv
title('Open Loop Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   


% Closed Loop Analysis

clpoles = eig(Aol2-Bol2*Col2)           % Closed Loop Poles
clpkf = eig(A2 - H2*C2)                 % Closed Loop Poles Due to Kalman Filter
clpreg = eig(A2 - B2*G_lqg2)            % Closed Loop Poles Due to Regulator


sv = sigma(ss(Aol2-Bol2*Col2, Bol2, -Col2, eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

sv = sigma(ss(Aol2-Bol2*Col2, Bol2, Col2, 0*eye(nc2)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values for OP2 (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
