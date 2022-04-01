function steady_state
    mol = 1.0e12;
    umol = 1.0e6;
    nmol = 1.0e3;
    pmol = 1.0e0;
    fmol = 1.0e-3;

    s = 1.0;
    l = 1.0;

    mU = 100.0 / 6.0 * nmol / s;

    % Data
    aT = 0.1 / l;
    aS = 0.4 / l;
    aS2 = 2.6 * 10^(-5) / l;
    a31 = 2.6 * 10^(-2) / l;
    a32 = 1.3 * 10^(-5) / l;
    
    BT = 1.1 * 10^(-6) / s;
    BS = 2.3 * 10^(-4) / s;
    BS2 = 140.0 / s;
    B31 = 8.0 * 10^(-6) / s;
    B32 = 8.3 * 10^(-4) / s;
    
    GT = 3.4 * pmol / s;
    GH = 817.0 * mU / s;
    GD1 = 22.0 * nmol / s;
    GD2 = 4.3 * fmol / s;
    GT3 = 394.0 * fmol / s;

    % correct value
    GR = 1.0 * mol / s;
    
    KM1 = 500.0 * nmol / l;

    % correct value
    KM2 = 1.0 * nmol / l;
    
    K30 = 2 * 10^(9) * l / mol;
    K41 = 2 * 10^(10) * l / mol;
    K42 = 2 * 10^(8) * l / mol;
    K31 = 2 * 10^(9) * l / mol;
    
    % correct value
    k = 1.0 * mU / l;
    
    DH = 47.0 * nmol / l;
    DS = 50.0 * mU / l;
    DT = 2.75 * mU / l;
    DR = 100.0 * pmol / l;
    
    SS = 100.0 * l / mU;
    LS = 1.68 * l / umol;
    
    TRH = 6.9 * nmol / l;
    TBG = 300.0 * nmol / l;
    TBPA = 4.5 * umol / l;
    IBS = 8.0 * umol / l;

    % Derived values
    a1 = aT * GT / BT;
    a3 = a32 * GD2 / B32;
    b1 = 1 / (1 + K41 * TBG + K42 * TBPA);
    
    c2 = 1 + KM2 / (a1 * b1);
    c3 = DT * KM2 / (a1 * b1);
    
    a5 = aS2 * GH / BS2;
    b3 = 1 / (1 + K31 * IBS);
    c1 = TRH / (TRH + DH);
    d2 = LS * GR * b3 * a3;
    
    d3 = b3 * a3 + c2 * DR;
    d4 = c3 * DR;
    d5 = d2 + d3;
    d6 = 1 + SS;
    f1 = c1 * a5;
    
    g1 = -d4 * d6;
    g2 = f1 * d4 - d4 * DS;
    g3 = f1 * d4 * DS;
    g4 = d5 * d6;
    g5 = -f1 * d3 + d5 * DS;
    
    g6 = -f1 * d3 * DS;
    g7 = aS2 * BS / (BS2 * aS);
    
    m1 = g4 * g7 * g7;
    m2 = g5 * g7 - g1 * g7 * g7;
    m3 = g6 - g2 * g7;
    m4 = -g3;

    TSH_roots = roots([m1 m2 m3 m4]);
    TSH = TSH_roots(2)

    T4 = aT * GT * TSH / (BT * (TSH + DT))

    FT4 = T4 / (1 + K41 * TBG + K42 * TBPA);

    T4th = GT * TSH / (TSH + DT);

    T3P = a31 / B31 * (GD1 * FT4 / (FT4 + KM1) + GD2 * FT4 / (FT4 + KM2) + GT3 * TSH / (TSH + DT) + GD1 * T4th * (TSH / (TSH + k)) / (KM1 + T4th * TSH / (TSH + k)) + GD2 * T4th * (TSH / TSH + k) / (KM2 + T4th * TSH / (TSH + k)))

    T3c = a32 * GD2 / B32 * FT4 / (FT4 + KM2)

    T3N = T3c / (1 + K31 * IBS);

    T3R = GR * T3N / (T3N + DR);

    x = aS2 * GH * TRH;
    y = TRH + DH;
    z = 1 + LS * T3R;
    TSHz_roots = roots([ ...
        1 + SS ...
        DS - x / (y * z * B32) ...
        -x / (y * z * B32) ...
    ]);
    TSHz = TSHz_roots(2)

    
end
