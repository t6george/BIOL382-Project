function thyroid

    % define constants
    
    aT = 0.1;
    aS = 0.4;
    aS2 = 2.6 * 10^(-5);
    a31 = 2.6 * 10^(-2);
    a32 = 1.3 * 10^(-5);
    
    BT = 1.1 * 10^(-6);
    BS = 2.3 * 10^(-4);
    BS2 = 140.0;
    B31 = 8.0 * 10^(-6);
    B32 = 8.3 * 10^(-4);
    
    GT = 3.4;
    GH = 817.0;
    GD1 = 22.0;
    GD2 = 4.3;
    GT3 = 394.0;

    % correct value
    GR = 1.0;
    
    KM1 = 500.0;

    % correct value
    KM2 = 1.0;
    
    K30 = 2 * 10^(9);
    K41 = 2 * 10^(10);
    K42 = 2 * 10^(8);
    K31 = 2 * 10^(9);
    
    % correct value
    k = 1.0;

    % what are these??
    l = 1.0;
    s = 1.0;
    
    DH = 47.0;
    DS = 50.0;
    DT = 2.75;
    DR = 100.0;
    
    SS = 100.0;
    LS = 1.68;
    
    TRH = 6.9;
    TBG = 300.0;
    TBPA = 4.5;
    IBS = 8.0;
    
    
    % set simulation parameters
    ODEFUN = @thyroidddt;
    Tend = 3600;
    
    % set initial condition: state = [T4, T3P, T3c, TSH, TSHz]
    T4_0 = 3.0909e+05;
    T3P_0 = 1.3026e+06;
    T3c_0 = 3.4689e-09;
    TSH_0 = 1.8189e+05;
    TSHz_0 = 0.0619;
    S0 = [T4_0, T3P_0, T3c_0, TSH_0, TSHz_0];


    % set delays
    T0T = 300.0;
    T03Z = 3600.0;
    T0S = 120.0;
    T0S2 = 3240.0;

    % opts = ddeset('RelTol', 1e-2, 'AbsTol', 1e-4);
    
    % simulate system
    sol=dde23( ...
        @(t, S, Z) thyroidddt(t, S, Z, aT, aS, aS2, a31, a32, BT, BS, BS2, B31, B32, GT, GH, GD1, GD2, GT3, GR, KM1, KM2, K30, K41, K42, K31, k, l, s, DH, DS, DT, DR, SS, LS, TRH, TBG, TBPA, IBS), ...
        [T0T, T03Z, T0S, T0S2], S0, [0, Tend] ...
    );

    % sol.y(end)

    % plot
    figure(1)
    set(gca, 'fontsize', 14)
    plot(sol.x, sol.y, 'LineWidth', 1.5)
    %axis([0 20 0 1])
    legend('T4', 'T3P', 'T3c', 'TSH', 'TSHz')
    xlabel('Time (s)')
    ylabel('Concentration (\muM)')
   
end


% model dynamics
function dS = thyroidddt(t, S, Z, aT, aS, aS2, a31, a32, BT, BS, BS2, B31, B32, GT, GH, GD1, GD2, GT3, GR, KM1, KM2, K30, K41, K42, K31, k, l, s, DH, DS, DT, DR, SS, LS, TRH, TBG, TBPA, IBS)

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

    T4_T03Z = Z(1, 2);
    FT4_T03Z = T4_T03Z / (1 + K41 * TBG + K42 * TBPA);

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

    dT4 = aT * GT * TSH_T0T / (TSH_T0T + DT) - BT * T4;

    dT3P = a31 * ((GD1 * FT4) / (FT4 + KM1) + ...
        (GD2 * FT4) / (FT4 + KM2) + ...
        (GT3 * TSH) / (TSH + DT) + ...
        (GD1 * (T4th * TSH) / (TSH + k)) / (KM1 + (T4th * TSH) / (TSH + k)) + ...
        (GD2 * (T4th * TSH) / (TSH + k)) / (KM2 + (T4th * TSH) / (TSH + k))) - B31 * T3P;

    dT3c = (a32 * GD2 * FT4_T03Z) / (FT4_T03Z + KM2) - B32 * T3c;

    dTSH = (aS * GH * TRH_T0S) / ((TRH_T0S + DH) * (1 + (SS * TSHz_T0S) / (TSHz_T0S + DS)) * (1 + LS * T3R_T0S)) - BS * TSH;

    dTSHz = (aS2 * GH * TRH_T0S2) / ((TRH_T0S2 + DH) * (1 + (SS * TSHz_T0S2) / (TSHz_T0S2 + DS)) * (1 + LS * T3R_T0S2)) - BS2 * TSHz;


    dS = [dT4; dT3P; dT3c; dTSH; dTSHz];

end
