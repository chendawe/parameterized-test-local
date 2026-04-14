
// /* Standard LAL */
// #include <lal/Sequence.h>
// #include <lal/Date.h>
// #include <lal/Units.h>
// #include <lal/XLALError.h>
 
// /* LAL datatypes and constants */
// #include <lal/LALDatatypes.h>
// #include <lal/LALStdlib.h>
// #include <lal/LALConstants.h>
 
// /* Time series, frequency series and spherical harmonics */
// #include <lal/TimeSeries.h>
// #include <lal/TimeFreqFFT.h>
// #include <lal/SphericalHarmonics.h>
// #include <lal/FrequencySeries.h>
 
/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>


// ==============

#include <lal/XLALError.h>
 
#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_internals.h"

static double IMRPhenomX_Inspiral_Phase_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{
    double phaseIN = 0.0;
 
    // Assemble PN phase derivative series
    phaseIN += pPhase->dphi_minus2 / powers_of_Mf->two_thirds;                   // f^{-2/3}, v =-2;  -1PN
    phaseIN += pPhase->dphi_minus1 / powers_of_Mf->one_third;                    // f^{-1/3}, v =-1;  -0.5PN
        phaseIN += pPhase->dphi0;                                                    // f^{0/3}
        phaseIN += pPhase->dphi1   * powers_of_Mf->one_third;                        // f^{1/3}
        phaseIN += pPhase->dphi2   * powers_of_Mf->two_thirds;                       // f^{2/3}
        phaseIN += pPhase->dphi3   * Mf;                                             // f^{3/3}
        phaseIN += pPhase->dphi4   * powers_of_Mf->four_thirds;                      // f^{4/3}
        phaseIN += pPhase->dphi5   * powers_of_Mf->five_thirds;                      // f^{5/3}
        phaseIN += pPhase->dphi6   * powers_of_Mf->two;                              // f^{6/3}
        phaseIN += pPhase->dphi6L  * powers_of_Mf->two * powers_of_Mf->log;          // f^{6/3}, Log[f]
        phaseIN += pPhase->dphi7   * powers_of_Mf->seven_thirds;                     // f^{7/3}
        phaseIN += pPhase->dphi8   * powers_of_Mf->eight_thirds;                     // f^{8/3}
        phaseIN += pPhase->dphi8L  * powers_of_Mf->eight_thirds * powers_of_Mf->log; // f^{8/3}
        phaseIN += pPhase->dphi9   * powers_of_Mf->three;                            // f^{9/3}
        phaseIN += pPhase->dphi9L  * powers_of_Mf->three * powers_of_Mf->log;        // f^{9/3}
 
        // Add pseudo-PN Coefficient
        phaseIN += (            pPhase->a0 * powers_of_Mf->eight_thirds
                                                                + pPhase->a1 * powers_of_Mf->three
                                                                + pPhase->a2 * powers_of_Mf->eight_thirds * powers_of_Mf->two_thirds
                                                                + pPhase->a3 * powers_of_Mf->eight_thirds * powers_of_Mf->itself
                                                                + pPhase->a4 * powers_of_Mf->eight_thirds * powers_of_Mf->four_thirds
                                                        );
 
        phaseIN  = phaseIN * powers_of_Mf->m_eight_thirds * (5.0 / (128.0 * powers_of_lalpi.five_thirds));
 
    return phaseIN;
}
 
/**
 * Ansatz for the inspiral phase.
 * The TaylorF2 coefficients are defined elsewhere.
 */
static double IMRPhenomX_Inspiral_Phase_22_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{
 
  // Assemble PN phasing series
  // Sum up phasing contributions
  double phasing    = 0.0;
 
  /* The PN Phasing series is normalised by: 3 / (128 * eta * pi^{5/3} ) */
  /* -2 */ phasing += pPhase->phi_minus2 / powers_of_Mf->two_thirds;                  // f^{-7/3}, v = -2; -1PN
  /* -1 */ phasing += pPhase->phi_minus1 / powers_of_Mf->one_third;                   // f^{-6/3}, v = -1; -0.5PN
  /* 0  */ phasing += pPhase->phi0;                                                   // f^{-5/3}, v = 0;  Newt.
  /* 1  */ phasing += pPhase->phi1  * powers_of_Mf->one_third;                        // f^{-4/3}, v = 1;  0.5PN
  /* 2  */ phasing += pPhase->phi2  * powers_of_Mf->two_thirds;                       // f^{-3/3}, v = 2;  1.0PN
  /* 3  */ phasing += pPhase->phi3  * Mf;                                             // f^{-2/3}, v = 3;  1.5PN
  /* 4  */ phasing += pPhase->phi4  * powers_of_Mf->four_thirds;                      // f^{-1/3}, v = 4;  2.0PN
  /* 5  */ phasing += pPhase->phi5  * powers_of_Mf->five_thirds;                      // f^{0},    v = 5;  2.5PN; phi_initial = - LAL_PI_4
  /* 5L */ phasing += pPhase->phi5L * powers_of_Mf->five_thirds * powers_of_Mf->log;  // f^{0},    v = 5;  2.5PN Log terms.
  /* 6  */ phasing += pPhase->phi6  * powers_of_Mf->two;                              // f^{+1/3}; v = 6;  3.0PN
  /* 6L */ phasing += pPhase->phi6L * powers_of_Mf->two * powers_of_Mf->log;          // f^{+1/3}; v = 6;  3.0PN Log terms.
  /* 7  */ phasing += pPhase->phi7  * powers_of_Mf->seven_thirds;                     // f^{+2/3}: v = 7;  3.5PN
  /* 8  */ phasing += pPhase->phi8  * powers_of_Mf->eight_thirds;                     // f^{+3/3}; v = 8;  4.0PN
  /* 8  */ phasing += pPhase->phi8L * powers_of_Mf->eight_thirds * powers_of_Mf->log; // f^{+3/3}; v = 8;  4.0PN Log terms.
  /* 9  */ phasing += pPhase->phi9  * powers_of_Mf->three;                            // f^{+4/3}; v = 9;  4.5PN
  /* 9  */ phasing += pPhase->phi9L * powers_of_Mf->three * powers_of_Mf->log;        // f^{+4/3}; v = 9;  4.5PN
 
  // Now add in the pseudo-PN Coefficients
  phasing += (  pPhase->sigma1 * powers_of_Mf->eight_thirds
              + pPhase->sigma2 * powers_of_Mf->three
              + pPhase->sigma3 * powers_of_Mf->one_third  * powers_of_Mf->three
              + pPhase->sigma4 * powers_of_Mf->two_thirds * powers_of_Mf->three
              + pPhase->sigma5 * powers_of_Mf->itself     * powers_of_Mf->three
            );
 
  // This completes the TaylorF2 PN phasing series
  phasing = phasing * pPhase->phiNorm * powers_of_Mf->m_five_thirds;
 
  /* Add initial phasing: -pi/4 */
  //phasing += pPhase->phi_initial;
 
  return phasing;
}
