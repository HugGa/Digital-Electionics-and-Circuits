#include <iostream>
#include <cassert>
#include <cmath>
namespace Hugh
{
    namespace DigitalElectronics
    {
        enum class TransistorPhase
        {
            SATURATED,
            TRIODE,
            OFF
        };
        enum class ActivationMode
        {
            NMOS,
            PMOS
        };
        template <class T>
        struct DigETuple
        {
            T state;
            double value;
            // specialize this per enum for DigETuple, may want to use MagicEnum here
            friend std::ostream &operator<<(std::ostream &os, const DigETuple<T> &tuple);
            bool operator==(DigETuple<T> b)
            {
                if ((this->value + .001) > b.value && (this->value - .001) < b.value)
                {
                    return true;
                }
                return false;
            }
        };
        std::ostream &operator<<(std::ostream &os, const DigETuple<TransistorPhase> &tuple)
        {
            os << "Phase of transistor: ";
            switch (tuple.state)
            {
            case TransistorPhase::SATURATED:
                os << "SATURATED";
                break;
            case TransistorPhase::TRIODE:
                os << "TRIODE";
                break;
            case TransistorPhase::OFF:
                os << "OFF";
                break;
            default:
                assert(false);
                break;
            }
            os << " The value of the current is: " << tuple.value << std::endl;
            return os;
        }
        // Vto == Vtn or Vtp
        DigETuple<TransistorPhase> currentShortChannelModelNMOS(double Vgs, double Vtn, double Cox, double Un, double W, double L, double Vds, double Van, double Ecn, double Ln)
        {
            return currentSCMNMOSNOVANNOCOX(Vgs, Vtn, Vds, Un * Cox, 1 / Van, W / L, Ecn * Ln);
        }
        DigETuple<TransistorPhase> currentShortChannelModelNMOS(double Vgs, double Vtn, double Cox, double Un, double W, double L, double Vds, double Van, double EcnLn)
        {
            return currentSCMNMOSNOVANNOCOX(Vgs, Vtn, Vds, Un * Cox, 1 / Van, W / L, EcnLn);
        }
        DigETuple<TransistorPhase> currentShortChannelModelNMOSNoVan(double Vgs, double Vtn, double Cox, double Un, double W, double L, double Vds, double lambdan, double EcnLn)
        {
            return currentSCMNMOSNOVANNOCOX(Vgs, Vtn, Vds, Un * Cox, lambdan, W / L, EcnLn);
        }
        // Short Channel Model using Lambdan and kprime and ecnln
        DigETuple<TransistorPhase> currentSCMNMOSNOVANNOCOX(double Vgs, double Vtn, double Vds, double kprime, double lambdan, double WL, double ecnln)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double Vgstn = Vgs - Vtn;
            double VDSsat = (Vgstn * ecnln) / (Vgstn + ecnln);
            double k = WL * kprime;
            if (VDSsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VDSsat << std::endl;
            }
            if (Vds <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vds << std::endl;
            }
            if (lambdan < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdan << std::endl;
            }
            if (Vgs < Vtn)
            {
                return returnval;
            }
            if (Vgs > Vtn && Vds <= VDSsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = WL * (kprime / (1 + (Vds / ecnln))) * ((Vgstn * Vds) - ((Vds * Vds) / 2));
                return returnval;
            }
            if (Vgs > Vtn && Vds >= VDSsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = ((k / 2) * ecnln) * ((Vgstn * Vgstn) / (Vgstn + ecnln)) * (1 + lambdan * (Vds - VDSsat));
                return returnval;
            }
            return returnval;
        }
        DigETuple<TransistorPhase> currentLongChannelNMOS(double Vgs, double Vtn, double Cox, double Un, double W, double L, double Vds, double Van)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double lambdan = 1 / (Van);
            double VDSsat = Vgs - Vtn;
            double temp = Un * Cox;
            double otemp = W / L;
            if (VDSsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VDSsat << std::endl;
            }
            if (Vds <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vds << std::endl;
            }
            if (lambdan < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdan << std::endl;
            }
            if (Vgs < Vtn)
            {
                returnval.state = TransistorPhase::OFF;
                returnval.value = 0;
                return returnval;
            }
            if (Vgs > Vtn && Vds <= VDSsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = temp * otemp * (((Vgs - Vtn) * Vds) - ((Vds * Vds) / 2));
                return returnval;
            }
            if (Vgs > Vtn && Vds >= VDSsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = (temp / 2) * otemp * (Vgs - Vtn) * (Vgs - Vtn) * (1 + lambdan * (Vds - VDSsat));
                return returnval;
            }
            return returnval;
        }

        DigETuple<TransistorPhase> currentShortChannelModelPMOS(double Vsg, double Vtp, double Cox, double Up, double W, double L, double Vsd, double Vap, double Ecp, double Lp)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double negVTP = -Vtp;
            double lambdap = 1 / Vap;
            double VsgTp = (Vsg + Vtp);
            double VSDsat = (VsgTp * Ecp * Lp) / (VsgTp + (Ecp * Lp));
            double randomVal = W * VSDsat * Cox;
            double Wl = W / L;
            if (VSDsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VSDsat << std::endl;
            }
            if (Vsd <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vsd << std::endl;
            }
            if (lambdap < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdap << std::endl;
            }
            if (Vsg < negVTP)
            {
                return returnval;
            }
            if (Vsg > negVTP && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = Wl * ((Up * Cox) / (1 + (Vsd / (Ecp * Lp))) * ((Vsg + Vtp) * Vsd - (Vsd * Vsd) / 2);
                return returnval;
            }
            if (Vsg > negVTP && Vsd >= VSDsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = randomVal * (((Vsg + Vtp) * (Vsg + Vtp)) / ((Vsg + Vtp) + Ecp * Lp)) * (1 + lambdap * (Vsd - VSDsat));
                return returnval;
            }
            return returnval;
        }
        DigETuple<TransistorPhase> currentShortChannelModelPMOS(double Vsg, double Vtp, double Cox, double Up, double W, double L, double Vsd, double Vap, double EcpLp)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double negVTP = -Vtp;
            double lambdap = 1 / Vap;
            double VsgTp = (Vsg + Vtp);
            double VSDsat = (VsgTp * EcpLp) / (VsgTp + (EcpLp));
            double randomVal = W * VSDsat * Cox;
            double Wl = W / L;
            if (VSDsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VSDsat << std::endl;
            }
            if (Vsd <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vsd << std::endl;
            }
            if (lambdap < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdap << std::endl;
            }
            if (Vsg < negVTP)
            {
                return returnval;
            }
            if (Vsg > negVTP && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = Wl * ((Up * Cox) / (1 + (Vsd / (EcpLp))) * ((Vsg + Vtp) * Vsd - (Vsd * Vsd) / 2);
                return returnval;
            }
            if (Vsg > negVTP && Vsd >= VSDsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = randomVal * (((Vsg + Vtp) * (Vsg + Vtp)) / ((Vsg + Vtp) + EcpLp)) * (1 + lambdap * (Vsd - VSDsat));
                return returnval;
            }
            return returnval;
        }
        DigETuple<TransistorPhase> currentShortChannelModelPMOSNoVan(double Vsg, double Vtp, double Cox, double Up, double W, double L, double Vsd, double lambdap, double EcpLp)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double negVTP = -Vtp;
            double VsgTp = (Vsg + Vtp);
            double VSDsat = (VsgTp * EcpLp) / (VsgTp + (EcpLp));
            double randomVal = W * VSDsat * Cox;
            double Wl = W / L;
            if (VSDsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VSDsat << std::endl;
            }
            if (Vsd <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vsd << std::endl;
            }
            if (lambdap < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdap << std::endl;
            }
            if (Vsg < negVTP)
            {
                return returnval;
            }
            if (Vsg > negVTP && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = Wl * ((Up * Cox) / (1 + (Vsd / (EcpLp))) * ((Vsg + Vtp) * Vsd - (Vsd * Vsd) / 2);
                return returnval;
            }
            if (Vsg > negVTP && Vsd >= VSDsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = randomVal * (((Vsg + Vtp) * (Vsg + Vtp)) / ((Vsg + Vtp) + EcpLp)) * (1 + lambdap * (Vsd - VSDsat));
                return returnval;
            }
            return returnval;
        }
        DigETuple<TransistorPhase> currentLongChannelModelPMOS(double Vsg, double Vtp, double Cox, double Up, double W, double L, double Vsd, double Vap)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double VSDsat = Vsg + Vtp;
            double lambdap = 1 / Vap;
            if (VSDsat <= 0)
            {
                std::cerr << "VDSsat is below 0, erroring, VDSSat: " << VSDsat << std::endl;
            }
            if (Vsd <= 0)
            {
                std::cerr << "VDs is below 0, erroring, Vds: " << Vsd << std::endl;
            }
            if (lambdap < 0)
            {
                std::cerr << "lambdan is below 0, erroring, lambdan:" << lambdap << std::endl;
            }
            if (Vsg < -Vtp)
            {
                return returnval;
            }
            if (Vsg > -Vtp && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = Up * Cox * (W / L) * (((Vsg + Vtp) * Vsd) - (Vsd * Vsd) / 2);
                return returnval;
            }
            if (Vsg > -Vtp && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = ((Up * Cox) / 2) * (W / L) * ((Vsg + Vtp) * (Vsg + Vtp)) * (1 + lambdap * (Vsd - VSDsat));
                return returnval;
            }
            return returnval;
        }
        DigETuple<TransistorPhase> currentSCMPMOSNOVANNOCOX(double Vsg, double Vtp, double Vsd, double kprime, double lambdan, double WL, double ecplp)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            double Vsgtp = Vsg + Vtp;
            double VSDsat = (Vsgtp * ecplp) / (Vsgtp + ecplp);
            double k = WL * kprime;
            if (Vsg < -Vtp)
            {
                return returnval;
            }
            if (Vsg > -Vtp && Vsd <= VSDsat)
            {
                returnval.state = TransistorPhase::TRIODE;
                returnval.value = WL * (kprime / (1 + (Vsd / ecplp))) * ((Vsgtp * Vsd) - ((Vsd * Vsd) / 2));
                return returnval;
            }
            if (Vsg > -Vtp && Vsd >= VSDsat)
            {
                returnval.state = TransistorPhase::SATURATED;
                returnval.value = ((k / 2) * ecplp) * ((Vsgtp * Vsgtp) / (Vsgtp + ecplp)) * (1 + lambdan * (Vsd - VSDsat));
                return returnval;
            }
            return returnval;
        }
        // Gamma = Body effect
        // Phif = some body effect thing
        // Vb = Voltage of Body
        // Vs = Voltage of supply
        // Vg = Voltage of the turn/off/on thing
        // Vd = Vdd being used
        DigETuple<TransistorPhase> currentShortChannel(ActivationMode mode, double Vg, double Vs, double Vd, double Vb, double gamma, double Phi, double Vt, double Cox, double Un, double W, double L, double Van, double Ecn, double Ln)
        {
            DigETuple<TransistorPhase> returnval = {TransistorPhase::OFF, 0};
            switch (mode)
            {
            case ActivationMode::NMOS:
            {
                double Vsb = Vs - Vb;
                double Vgs = Vg - Vs;
                double Vtn = Vt + (Vsb == 0 ? 0 : gamma * (sqrt(fabs(2 * Phi) + Vsb) - sqrt(fabs(2 * Phi))));
                double Vds = Vd - Vs;
                return currentShortChannelModelNMOS(Vgs, Vtn, Cox, Un, W, L, Vds, Van, Ecn, Ln);
            }
            break;
            case ActivationMode::PMOS:
            {
                double Vbs = Vb - Vs;
                double Vsg = Vs - Vg;
                double Vtp = Vt + (Vbs == 0 ? 0 : -gamma * (sqrt(fabs(2 * Phi) + Vbs) - sqrt(fabs(2 * Phi))));
                double Vsd = Vd - Vs;
                return currentShortChannelModelPMOS(Vsg, Vtp, Cox, Un, W, L, Vsd, Van, Ecn, Ln);
            }
            break;
            default:
                break;
            }
            return returnval;
        }
    }

}