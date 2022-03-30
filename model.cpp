#include "ascent/Ascent.h"

using namespace asc;


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

    static constexpr double dt = 0.01;

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
    state_t x = { 10.0, 1.0, 1.0 };
    double t = 0.0;
    constexpr double t_end = 10.0;

    const auto constants = Constants();
    const auto delay = Delay();

    auto lorenz = [&constants = std::as_const(constants), &delay = std::as_const(delay)](const state_t& x, state_t& xd, const double)
    {
        static constexpr double sigma = 10.0;
        static constexpr double R = 28.0;
        static constexpr double b = 8.0 / 3.0;

        xd[0] = sigma * (x[1] - x[0]);
        xd[1] = R * x[0] - x[1] - x[0] * x[2];
        xd[2] = -b * x[2] + x[0] * x[1];
    };

   RK4 integrator;
   Recorder recorder;

   while (t < t_end)
   {
      recorder({ t, x[0], x[1], x[2] });
      integrator(lorenz, x, t, Delay::dt);
   }

   recorder.csv("lorenz", { "t", "x0", "x1", "x2" });

   return 0;
}
