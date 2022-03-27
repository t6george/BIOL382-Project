function thyroid

global aT;
global aS;
global aS2;
global a31;
global a32;

global BT;
global BS;
global BS2;
global B31;
global B32;

global GT;
global GH;
global GD1;
global GD2;
global GT3;
global GR;

global KM1;
global KM2;
global K30;
global K41;
global K42;
global K31;

global k;
global l;
global s;

global DH;
global DS;
global DT;
global DR;

global SS;
global LS;

global TRH;
global TBG;
global TBPA;
global IBS;


aT = 1.0;
aS = 1.0;
aS2 = 1.0;
a31 = 1.0;
a32 = 1.0;

BT = 1.0;
BS = 1.0;
BS2 = 1.0;
B31 = 1.0;
B32 = 1.0;

GT = 1.0;
GH = 1.0;
GD1 = 1.0;
GD2 = 1.0;
GT3 = 1.0;
GR = 1.0;

KM1 = 1.0;
KM2 = 1.0;
K30 = 1.0;
K41 = 1.0;
K42 = 1.0;
K31 = 1.0;

k = 1.0;
l = 1.0;
s = 1.0;

DH = 1.0;
DS = 1.0;
DT = 1.0;
DR = 1.0;

SS = 1.0;
LS = 1.0;

TRH = 1.0;
TBG = 1.0;
TBPA = 1.0;
IBS = 1.0;


%set simulation parameters
ODEFUN = @thyroidddt;
Tend = 20;

%set initial condition: state = [T4, T3P, T3c, TSH, TSHz]
T4_0 = 0.0;
T3P_0 = 0.0;
T3c_0 = 0.0;
TSH_0 = 0.0;
TSHz_0 = 0.0;
S0 = [T4_0, T3P_0, T3c_0, TSH_0, TSHz_0];



T0T = 1.0;
T03Z = 1.0;
T0S = 1.0;
T0S2 = 1.0;

sol=dde23(ODEFUN, [T0T, T03Z, T0S, T0S2], S0, [0, Tend]);


figure(1)
set(gca, 'fontsize', 14)
plot(sol.x, sol.y, 'LineWidth', 1.5)
axis([0 200 0 8])
%legend('mRNA (M)', 'total PER (P_T)', 'nuclear PER (P_N)')
xlabel('Time (h)')
ylabel('Concentration (\muM)')
   
end


%model dynamics
function dS = thyroidddt(t, S, Z)

    global aT;
    global aS;
    global aS2;
    global a31;
    global a32;
    
    global BT;
    global BS;
    global BS2;
    global B31;
    global B32;
    
    global GT;
    global GH;
    global GD1;
    global GD2;
    global GT3;
    global GR;
    
    global KM1;
    global KM2;
    global K30;
    global K41;
    global K42;
    global K31;
    
    global k;
    global l;
    global s;
    
    global DH;
    global DS;
    global DT;
    global DR;
    
    global SS;
    global LS;
    
    global TRH;
    global TBG;
    global TBPA;
    global IBS;

    % Current State Variables
    T4 = S(1);
    T3P = S(2);
    T3c = S(3);
    TSH = S(4);
    TSHz = S(5);

    T4th = GT * TSH / (TSH + DT) * s / l;
    FT3 = T3P / (1 + K30 * TBG);
    FT4 = T4 / (1 + K41 * TBG + K42 * TBPA);
    T3N = T3c / (1 + K31 * IBS);
    T3R = GR * T3N / (T3N + DR);


    % Delayed State Variables
    % [T0T, T03Z, T0S, T0S2]
    TSH_T0T = Z(4, 1);

    FT4_T03Z = Z(8, 2);

    TRH_T0S = TRH;
    TSHz_T0S = Z(5, 3);
    T3c_T0S = Z(3, 3);
    T3N_T0S = T3c_T0S / (1 + K31 * IBS);
    T3R_T0S = GR * T3N_T0S / (T3N_T0S + DR);


    TRH_T0S2 = TRH;
    TSHz_T0S2 = Z(5, 4);
    T3c_T0S2 = Z(3, 4);
    T3N_T0S2 = T3c_T0S2 / (1 + K31 * IBS);
    T3R_T0S2 = GR * T3N_T0S2 / (T3N_T0S2 + DR);

    
    % Differential Equations





    p0_lag = Z(1,1);
    
    phos_rate = 2 * p0_lag;
    
    dS=[vs/(1+(PN/ki)^n) - (vm*M)/(km1+M);
        ks*M - phos_rate;
        0;
        (v3*P1)/(k3+P1) - kt1*P2 + kt2*PN - (vd*P2)/(kd+P2) + phos_rate;
        kt1*P2 - kt2*PN
    ];

end
