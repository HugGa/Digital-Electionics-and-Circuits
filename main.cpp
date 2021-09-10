#include "circuits2.hpp"
#include "DigitalElec.hpp"
#include "conversion.hpp"
using namespace Hugh::Circuits2;
using namespace Hugh::DigitalElectronics;

int main()
{
    // Vs is ALWAYS the terminal at the lower voltage
    // Vt is threshold voltage
    // Vtp == Threshold voltage
    // Vsb = Vs - Vb
    // Vb not shown so it's tied to the lowest supply voltage for NMOS, and HIGHEST for PMOS
    // If lambdan is not given, take it as 0
    double Vg = 1;
    double Vs = 0;
    double Vb = 0;
    double Vd = .24806;
    double gamma = .3;
    double Phi = (.4);
    double Vt = .3;
    double Vsb = Vs - Vb;
    double Vgs = Vg - Vs;
    double Vtn = Vt + (Vsb == 0 ? 0 : gamma * (sqrt(fabs(2 * Phi) + Vsb) - sqrt(fabs(2 * Phi))));
    double Vds = Vd - Vs;
    std::cout << currentSCMNMOSNOVANNOCOX(Vgs, Vtn, Vds, 200, 0, 3, .6) << std::endl;
}
