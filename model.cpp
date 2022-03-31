#include <iostream>
#include "ascent/Ascent.h"


using namespace asc;

static constexpr int c_Precision = 9;


struct Constants
{
    const double aT = 0.1;
    const double aS = 0.4;
    const double aS2 = 2.6e-5;
    const double a31 = 2.6e-2;
    const double a32 = 1.3e-5;
    
    const double BT = 1.1e-6;
    const double BS = 2.3e-4;
    const double BS2 = 140.0;
    const double B31 = 8.0e-6;
    const double B32 = 8.3e-4;
    
    const double GT = 3.4;
    const double GH = 817.0;
    const double GD1 = 22.0;
    const double GD2 = 4.3;
    const double GT3 = 394.0;
    const double GR = 1.0;
    
    const double KM1 = 500.0;
    const double KM2 = 1.0;
    
    const double K30 = 2e9;
    const double K41 = 2e10;
    const double K42 = 2e8;
    const double K31 = 2e9;
    
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


class DelayedState
{
public:
    DelayedState(const double dt, const double TSH_0, const double TSHz_0, const double FT4_0, const double T3R_0)
    {
        TSH_history.resize(static_cast<size_t>(T0T / dt) + 2);
        TSHz_history.resize(static_cast<size_t>(std::max(T0S, T0S2) / dt) + 2);
        FT4_history.resize(static_cast<size_t>(T03Z / dt) + 2);
        T3R_history.resize(static_cast<size_t>(std::max(T0S, T0S2) / dt) + 2);

        TSH_history[0] = TSH_0;
        TSHz_history[0] = TSHz_0;
        FT4_history[0] = FT4_0;
        T3R_history[0] = T3R_0;
    }

    void update_history(const double t, const double TSH, const double TSHz, const double FT4, const double T3R) noexcept
    {
        TSH_history[TSH_idx] = TSH;
        TSHz_history[TSHz_idx] = TSHz;
        FT4_history[FT4_idx] = FT4;
        T3R_history[T3R_idx] = T3R;

        TSH_idx = (TSH_idx + 1) % TSH_history.size();
        TSHz_idx = (TSHz_idx + 1) % TSHz_history.size();
        FT4_idx = (FT4_idx + 1) % FT4_history.size();
        T3R_idx = (T3R_idx + 1) % T3R_history.size();

        if (t >= T0T) ++T0T_idx;
        if (t >= T03Z) ++T03Z_idx;
        if (t >= T0S) ++T0S_idx;
        if (t >= T0S2) ++T0S2_idx;
    }

    double getTSH_T0T() const noexcept { return TSH_history[T0T_idx % TSH_history.size()]; }

    double getTSHz_T0S() const noexcept { return TSHz_history[T0S_idx % TSHz_history.size()]; }

    double getTSHz_T0S2() const noexcept { return TSHz_history[T0S2_idx % TSHz_history.size()]; }

    double getFT4_T03Z() const noexcept { return FT4_history[T03Z_idx % FT4_history.size()]; }

    double getT3R_T0S() const noexcept { return T3R_history[T0S_idx % T3R_history.size()]; }

    double getT3R_T0S2() const noexcept { return T3R_history[T0S2_idx % T3R_history.size()]; }


private:

    std::vector<double> TSH_history;
    std::vector<double> TSHz_history;
    std::vector<double> FT4_history;
    std::vector<double> T3R_history;

    const double T0T = 300.0;
    const double T03Z = 3600.0;
    const double T0S = 120.0;
    const double T0S2 = 3240.0;

    size_t T0T_idx = 0;
    size_t T03Z_idx = 0;
    size_t T0S_idx = 0;
    size_t T0S2_idx = 0;

    size_t TSH_idx = 1;
    size_t TSHz_idx = 1;
    size_t FT4_idx = 1;
    size_t T3R_idx = 1;
};


struct CurrentState
{
    double T4 = 0.0;
    double T3P = 0.0;
    double T3c = 0.0;
    double TSH = 0.0;
    double TSHz = 0.0;
    double T4th = 0.0;
    double FT3 = 0.0;
    double FT4 = 0.0;
    double T3N = 0.0;
    double T3R = 0.0;

    void populate(const double T4, const double T3P, const double T3c, const double TSH, const double TSHz, const Constants& cs) noexcept
    {
        this->T4 = T4;
        this->T3P = T3P;
        this->T3c = T3c;
        this->TSH = TSH;
        this->TSHz = TSHz;

        T4th = cs.GT * TSH / (TSH + cs.DT);
        FT3 = T3P / (1 + cs.K30 * cs.TBG);
        FT4 = T4 / (1 + cs.K41 * cs.TBG + cs.K42 * cs.TBPA);
        T3N = T3c / (1 + cs.K31 * cs.IBS);
        T3R = cs.GR * T3N / (T3N + cs.DR);
    }
};


int main()
{
    // time params
    double t = 0.0;
    constexpr double dt = 0.01;
    constexpr double t_end = 60.0 * 60.0 * 24.0 * 25.0;

    const auto cs = Constants();

    // initial conditions
    auto curr_state = CurrentState();
    const double T4_0 = 3.0909e+05 * 0.85;
    const double T3P_0 = 1.3026e+06 * 0.85;
    const double T3c_0 = 3.4689e-09 * 0.85;
    const double TSH = 1.8189e+05 * 0.85;
    const double TSHz = 0.0619 * 0.85;
    curr_state.populate(T4_0, T3P_0, T3c_0, TSH, TSHz, cs);

    // delay state
    auto delay = DelayedState(dt, curr_state.TSH, curr_state.TSHz, curr_state.FT4, curr_state.T3R);


    // model
    auto thyroid = [&cs = std::as_const(cs), &t = std::as_const(t), &state = curr_state, &delay]
        (const state_t& x, state_t& xd, const double)
    {
        // current state variables
        state.populate(x[0], x[1], x[2], x[3], x[4], cs);

        // delay state variables
        const auto TSH_T0T = delay.getTSH_T0T();
        const auto TSHz_T0S = delay.getTSHz_T0S();
        const auto TSHz_T0S2 = delay.getTSHz_T0S2();
        const auto FT4_T03Z = delay.getFT4_T03Z();
        const auto T3R_T0S = delay.getT3R_T0S();
        const auto T3R_T0S2 = delay.getT3R_T0S2();

        // update history state
        delay.update_history(t, state.TSH, state.TSHz, state.FT4, state.T3R);

        // differential equations
        const auto dT4 = cs.aT * cs.GT * TSH_T0T / (TSH_T0T + cs.DT) - cs.BT * state.T4;

        const auto dT3P = cs.a31 * ((cs.GD1 * state.FT4) / (state.FT4 + cs.KM1) +
            (cs.GD2 * state.FT4) / (state.FT4 + cs.KM2) +
            (cs.GT3 * state.TSH) / (state.TSH + cs.DT) +
            (cs.GD1 * (state.T4th * state.TSH) / (state.TSH + cs.k)) / (cs.KM1 + (state.T4th * state.TSH) / (state.TSH + cs.k)) +
            (cs.GD2 * (state.T4th * state.TSH) / (state.TSH + cs.k)) / (cs.KM2 + (state.T4th * state.TSH) / (state.TSH + cs.k))) -
            cs.B31 * state.T3P;

        const auto dT3c = (cs.a32 * cs.GD2 * FT4_T03Z) / (FT4_T03Z + cs.KM2) - cs.B32 * state.T3c;

        const auto dTSH = (cs.aS * cs.GH * cs.TRH) / ((cs.TRH + cs.DH) * (1 + (cs.SS * TSHz_T0S) / (TSHz_T0S + cs.DS)) *
            (1 + cs.LS * T3R_T0S)) - cs.BS * state.TSH;

        const auto dTSHz = (cs.aS2 * cs.GH * cs.TRH) / ((cs.TRH + cs.DH) * (1 + (cs.SS * TSHz_T0S2) / (TSHz_T0S2 + cs.DS)) *
            (1 + cs.LS * T3R_T0S2)) - cs.BS2 * state.TSHz;

        xd[0] = dT4;
        xd[1] = dT3P;
        xd[2] = dT3c;
        xd[3] = dTSH;
        xd[4] = dTSHz;
    };

    auto integrator = RK4();
    auto recorder = Recorder();
    recorder.precision = c_Precision;

    auto x = state_t({curr_state.T4, curr_state.T3P, curr_state.T3c, curr_state.TSH, curr_state.TSHz});

    unsigned num_iters = 0;
    constexpr unsigned iter_sample_size = 100;

    while (t < t_end)
    {
	    if ((num_iters++) % iter_sample_size == 0)
	    {
            recorder({
                t,
                curr_state.T4,
                curr_state.T3P, 
                curr_state.T3c, 
                curr_state.TSH, 
                curr_state.TSHz, 
                curr_state.T4th, 
                curr_state.FT3, 
                curr_state.FT4, 
                curr_state.T3N, 
                curr_state.T3R
            });
	    }

        integrator(thyroid, x, t, dt);
    }

    recorder.csv("thyroid", {"t", "T4", "T3P", "T3c", "TSH", "TSHz", "T4th", "FT3", "FT4", "T3N", "T3R"});

    std::cout << "The simulation has completed!" << std::endl;

    return 0;
}
