function solve_m_t3p
    global mode
    mode = "F";
    options = optimoptions('fsolve', 'FunctionTolerance', 1e-18);

    if mode == "TSI_DT"
        fsolve(@system, [1])
    elseif mode == "TSI_k"
        fsolve(@system, [1, 1])
    elseif mode == "F"
        fsolve(@system, [1])
    elseif mode == "TSHz"
        fsolve(@system, [1])
    end
    
end

function F = system(m)
    global mode

    % Define parameters here
    % Omitted for brevity

    TSH = 0.0038;
    FT4 = 4.3e-11;
    FT3 = 2.4e-11;
    TSHz = 1.1810; % Solved numerically

    m1 = m(1);

    T4 = FT4 * (1 + K41 * TBG + K42 * TBPA);
    T3P = FT3 * (1 + K30 * TBG);
    T3c = a32 * GD2 / B32 * FT4 / (FT4 + KM2);
    T3N = T3c / (1 + K31 * IBS);
    T3R = GR * T3N / (T3N + DR);

    T4th = GT * (m1 + TSH / (TSH + DT));

    if mode == "TSI_DT"
        % dT4
        F(1) = (aT * GT * (TSH / (TSH + DT) + m1) - BT * T4) * 1e15;
    elseif mode == "TSI_k"
        m2 = m(2);

        % dT4
        F(1) = (aT * GT * (TSH / (TSH + DT) + m1) - BT * T4) * 1e15;
        % dT3P
        F(2) = (a31 * (GD1 * FT4 / (FT4 + KM1) + GD2 * FT4 / (FT4 + KM2) + GT3 * (m1 + TSH / (TSH + DT)) + GD1 * T4th * (m2 + TSH / (TSH + k)) / (KM1 + T4th * (m2 + TSH / (TSH + k))) + GD2 * T4th * (m2 + TSH / (TSH + k)) / (KM2 + T4th * (m2 + TSH / (TSH + k)))) - B31 * T3P) * 1e15;
    elseif mode == "F"
        % kft3
        F(1) = (aS * GH * TRH / ((TRH + DH) * (1 + SS * TSHz / (TSHz + DS)) * (1 + LS * T3R) * m(1)) - BS * TSH) * 1e4
    elseif mode == "TSHz"
        % Solve for TSHz
        TSHz = m(1);
        F(1) = aS2 * GH * TRH / ((TRH + DH) * (1 + SS * TSHz / (TSHz + DS)) * (1 + LS * T3R)) - BS2 * TSHz
    end
end












