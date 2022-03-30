#include "ascent/Ascent.h"

using namespace asc;

static constexpr size_t c_NumStateVariables = 10;
using FullState = std::array<double, c_NumStateVariables>;

struct Constants
{
    const double aT = 0.1;
    const double aS = 0.4;
    const double aS2 = 2.6 / pow(10.0, 5);
    const double a31 = 2.6 / pow(10.0, 2);
    const double a32 = 1.3 / pow(10.0, 5);
    
    const double BT = 1.1 / pow(10.0, 6);
    const double BS = 2.3 / pow(10.0, 4);
    const double BS2 = 140.0;
    const double B31 = 8.0 / pow(10.0, 6);
    const double B32 = 8.3 / pow(10.0, 4);
    
    const double GT = 3.4;
    const double GH = 817.0;
    const double GD1 = 22.0;
    const double GD2 = 4.3;
    const double GT3 = 394.0;
    const double GR = 1.0;
    
    const double KM1 = 500.0;
    const double KM2 = 1.0;
    
    const double K30 = 2 * pow(10.0, 9);
    const double K41 = 2 * pow(10.0, 10);
    const double K42 = 2 * pow(10.0, 8);
    const double K31 = 2 * pow(10.0, 9);
    
    const double k = 1.0;
    
    const double DH = 47.0;
    const double DS = 50.0;
    const double DT = 2.75;
    const double DR = 100.0;
    
    const double SS = 100.0;
    const double LS = 1.68;
    
    const double TRH = 6.9;
    const double TBG = 300.0;
    const double TBPA = 4.5;
    const double IBS = 8.0;
};

class Delay
{
public:

    static constexpr double dt = 0.1;

    size_t getT0TIndexOffset() const noexcept { return T0T_index_offset; }

    size_t getT03ZIndexOffset() const noexcept { return T03Z_index_offset; }

    size_t getT0SIndexOffset() const noexcept { return T0S_index_offset; }

    size_t getT0S2IndexOffset() const noexcept { return T0S2_index_offset; }

private:

    const double T0T = 300.0;
    const double T03Z = 3600.0;
    const double T0S = 120.0;
    const double T0S2 = 3240.0;

    const size_t T0T_index_offset = static_cast<size_t>(T0T / dt);
    const size_t T03Z_index_offset = static_cast<size_t>(T03Z / dt);
    const size_t T0S_index_offset = static_cast<size_t>(T0S / dt);
    const size_t T0S2_index_offset = static_cast<size_t>(T0S2 / dt);
};


int main()
{
    const auto constants = Constants();
    const auto delay = Delay();
    auto history = std::vector<FullState>();

    auto lorenz = [&cs = std::as_const(constants), &delay = std::as_const(delay), &history]
        (const state_t& x, state_t& xd, const double)
    {
        // static constexpr double sigma = 10.0;
        // static constexpr double R = 28.0;
        // static constexpr double b = 8.0 / 3.0;

        // xd[0] = sigma * (x[1] - x[0]);
        // xd[1] = R * x[0] - x[1] - x[0] * x[2];
        // xd[2] = -b * x[2] + x[0] * x[1];

        // current state variables
        const auto T4 = x[0];                                                   // history[_][0]
        const auto T3P = x[1];                                                  // history[_][1]
        const auto T3c = x[2];                                                  // history[_][2]
        const auto TSH = x[3];                                                  // history[_][3]
        const auto TSHz = x[4];                                                 // history[_][4]

        const auto T4th = cs.GT * TSH / (TSH + cs.DT);                          // history[_][5]
        const auto FT3 = T3P / (1 + cs.K30 * cs.TBG);                           // history[_][6]
        const auto FT4 = T4 / (1 + cs.K41 * cs.TBG + cs.K42 * cs.TBPA);         // history[_][7]
        const auto T3N = T3c / (1 + cs.K31 * cs.IBS);                           // history[_][8]
        const auto T3R = cs.GR * T3N / (T3N + cs.DR);                           // history[_][9]

        // delayed state variables
        const auto CURR_IDX = history.size();
        
        const auto T0T_IDX = std::max(static_cast<size_t>(0), CURR_IDX - delay.getT0TIndexOffset());
        const auto T03Z_IDX = std::max(static_cast<size_t>(0), CURR_IDX - delay.getT03ZIndexOffset());
        const auto T0S_IDX = std::max(static_cast<size_t>(0), CURR_IDX - delay.getT0SIndexOffset());
        const auto TOS2_IDX = std::max(static_cast<size_t>(0), CURR_IDX - delay.getT0S2IndexOffset());

        const auto TSH_T0T = history[T0T_IDX][0];
        const auto TSHz_T0S = history[T0S_IDX][4];
        const auto TSHz_T0S2 = history[TOS2_IDX][4];
        const auto FT4_T03Z = history[T03Z_IDX][7];
        const auto T3R_T0S = history[T0S_IDX][9];
        const auto T3R_T0S2 = history[TOS2_IDX][9];

        // save history
        history.emplace_back(FullState({T4, T3P, T3c, TSH, TSHz, T4th, FT3, FT4, T3N, T3R}));

        // differential equations
        const auto dT4 = cs.aT * cs.GT * TSH_T0T / (TSH_T0T + cs.DT) - cs.BT * T4;

        const auto dT3P = cs.a31 * ((cs.GD1 * FT4) / (FT4 + cs.KM1) +
            (cs.GD2 * FT4) / (FT4 + cs.KM2) +
            (cs.GT3 * TSH) / (TSH + cs.DT) +
            (cs.GD1 * (T4th * TSH) / (TSH + cs.k)) / (cs.KM1 + (T4th * TSH) / (TSH + cs.k)) +
            (cs.GD2 * (T4th * TSH) / (TSH + cs.k)) / (cs.KM2 + (T4th * TSH) / (TSH + cs.k))) -
            cs.B31 * T3P;

        const auto dT3c = (cs.a32 * cs.GD2 * FT4_T03Z) / (FT4_T03Z + cs.KM2) - cs.B32 * T3c;

        const auto dTSH = (cs.aS * cs.GH * cs.TRH) / ((cs.TRH + cs.DH) * (1 + (cs.SS * TSHz_T0S) / (TSHz_T0S + cs.DS)) *
            (1 + cs.LS * T3R_T0S)) - cs.BS * TSH;

        const auto dTSHz = (cs.aS2 * cs.GH * cs.TRH) / ((cs.TRH + cs.DH) * (1 + (cs.SS * TSHz_T0S2) / (TSHz_T0S2 + cs.DS)) *
            (1 + cs.LS * T3R_T0S2)) - cs.BS2 * TSHz;

        xd[0] = dT4;
        xd[1] = dT3P;
        xd[2] = dT3c;
        xd[3] = dTSH;
        xd[4] = dTSHz;
    };

    auto integrator = RK4();
    auto recorder = Recorder();

    constexpr double t_end = 10.0;
    history.reserve(static_cast<size_t>(t_end / Delay::dt) + 1);

    const auto x0 = FullState({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    history.emplace_back(x0);

    auto x = state_t(x0.begin(), x0.begin() + 5);

    double t = 0.0;

    while (t < t_end)
    {
        recorder({ t, x[0], x[1], x[2] });
        integrator(lorenz, x, t, Delay::dt);
    }

    recorder.csv("lorenz", { "t", "x0", "x1", "x2" });

    return 0;
}
