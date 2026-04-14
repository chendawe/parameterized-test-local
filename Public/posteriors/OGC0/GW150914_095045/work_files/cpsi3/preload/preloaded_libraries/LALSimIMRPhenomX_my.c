/***********************************************************************************************************************************************
LALSimIMRPhenomX.c 的注释；


 **********************************************************************************************************************************************/

/* Standard LAL */
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
 
/* LAL datatypes and constants */
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
 
/* Time series, frequency series and spherical harmonics */
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SphericalHarmonics.h>
#include <lal/FrequencySeries.h>
 
/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
 
/* Standard C */
#include <math.h>
#include <complex.h>
#include <stdbool.h>
 
#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#define PHENOMXDEBUG 0
#define PHENOMXPDEBUG 0
#else
#define DEBUG 1 //print debugging info
#define PHENOMXDEBUG 1
#define PHENOMXPDEBUG 1
#endif
 
/* Link IMRPhenomX routines */
#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_ringdown.h"
#include "LALSimIMRPhenomX_intermediate.h"
#include "LALSimIMRPhenomX_inspiral.h"
#include "LALSimIMRPhenomX_internals.c"
#include "LALSimIMRPhenomX_precession.c"
#include "LALSimIMRPhenomX_PNR.c"
#include "LALSimIMRPhenomX_AntisymmetricWaveform.c"
 
/* Note: This is declared in LALSimIMRPhenomX_internals.c and avoids namespace clashes */
IMRPhenomX_UsefulPowers powers_of_lalpi;
 
#ifndef _OPENMP
#define omp ignore
#endif


/* *********************************************************************************
  *
  * The following private function generates an IMRPhenomX frequency-domain waveform
  *   - Only aligned-spin
  *   - Only the 22-mode
  *   - Physical parameters are passed via the waveform struct
  * *********************************************************************************
  */
int IMRPhenomXASGenerateFD(
  COMPLEX16FrequencySeries **htilde22, /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,       /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
  LALDict *lalParams                   /**< LAL Dictionary Structure    */
)
{
  /* Inherits debug flag from waveform struct */
  UINT4 debug = PHENOMXDEBUG;
 
  /* Set tidal version */
  NRTidal_version_type NRTidal_version;
 
  NRTidal_version=IMRPhenomX_SetTidalVersion(lalParams);
 
  if(debug)
  {
    printf("\n **** Now in IMRPhenomXASGenerateFD... **** \n");
  }
 
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}
 
  /* Initialize useful powers of LAL_PI */
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");
 
  /* Inherit minimum and maximum frequencies to generate wavefom from input frequency grid */
  double f_min = freqs_In->data[0];
  double f_max = freqs_In->data[freqs_In->length - 1];
 
  /* Size of array */
  size_t npts     = 0;
 
  /* Index shift between freqs and the frequency series */
  UINT4 offset    = 0;
 
  /* Initialize frequency sequence */
  REAL8Sequence *freqs = NULL;
 
  /* Initialize tidal corrections (following tidal code modified from LALSimIMRPhenomP.c) */
 
  REAL8Sequence *phi_tidal = NULL;
  REAL8Sequence *amp_tidal = NULL;
  REAL8Sequence *planck_taper = NULL;
 
 
  /* Set matter parameters (set to zero in pWF if NRTidal additions are not turned on) */
  REAL8 lambda1 = pWF->lambda1;
  REAL8 lambda2 = pWF->lambda2;
 
  /* New variables needed for the NRTidalv2 model */
  REAL8 X_A = pWF->m1; // Already scaled by Mtot
  REAL8 X_B = pWF->m2; // Ibid.
  REAL8 pfaN = 3./(128.*X_A*X_B);
 
  /* If deltaF is non-zero then we need to generate a uniformly sampled frequency grid of spacing deltaF. Start at f = 0. */
  if(pWF->deltaF > 0)
  {
    /* Return the closest power of 2 */
    npts = (size_t) (f_max/pWF->deltaF) + 1;
 
    /* Debug information */
    if(debug)
    {
      printf("npts     = %zu\n",npts);
      printf("fMin     = %.4f\n",f_min);
      printf("fMax     = %.4f\n",f_max);
      printf("dF       = %.4f\n",pWF->deltaF);
    }
 
    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF ), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", pWF->deltaF);
 
    /* Initialize the htilde frequency series */
    *htilde22 = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
 
    /* Check that frequency series generated okay */
    XLAL_CHECK(*htilde22,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);
 
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;
 
    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
          "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);
 
    /* Allocate memory for frequency array and terminate if this fails */
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
 
    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      freqs->data[i-iStart] = i * pWF->deltaF;
    }
    offset = iStart;
  }
  else
  {
    /* freqs is a frequency grid with non-uniform spacing, so we start at the lowest given frequency */
    npts      = freqs_In->length;
    *htilde22 = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform, 22 mode", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);
 
    XLAL_CHECK (*htilde22, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", npts);
 
    offset = 0;
    freqs  = XLALCreateREAL8Sequence(freqs_In->length);
 
    /* Allocate memory for frequency array and terminate if this fails */
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
 
    /* Populate frequency array */
    for (UINT4 i = 0; i < freqs_In->length; i++)
    {
      freqs->data[i] = freqs_In->data[i];
    }
  }
 
  memset((*htilde22)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde22)->sampleUnits), &((*htilde22)->sampleUnits), &lalSecondUnit);
 
  /* Check if LAL dictionary exists. If not, create a LAL dictionary. */
  INT4 lalParams_In = 0;
  if(lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams = XLALCreateDict();
  }
 
  if(debug)
  {
    printf("\n\n **** Initializing amplitude struct... **** \n\n");
  }
 
  /* Allocate and initialize the PhenomX 22 amplitude coefficients struct */
  IMRPhenomXAmpCoefficients *pAmp22;
  pAmp22 = XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  status = IMRPhenomXGetAmplitudeCoefficients(pWF,pAmp22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAmplitudeCoefficients failed.\n");
 
  if(debug)
  {
    printf("\n\n **** Amplitude struct initialized. **** \n\n");
    printf("\n\n **** Initializing phase struct... **** \n\n");
  }
 
  /* Allocate and initialize the PhenomX 22 phase coefficients struct */
  IMRPhenomXPhaseCoefficients *pPhase22;
  pPhase22 = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhase22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");
  // initialize coefficients of tidal phase
  if(NRTidal_version!=NoNRT_V) IMRPhenomXGetTidalPhaseCoefficients(pWF,pPhase22,NRTidal_version);
 
  if(debug)
  {
    printf("\n\n **** Phase struct initialized. **** \n\n");
  }
 
  /*
      Apply time shifts so peak amplitude is near t ~ 0.
  */
 
  /* Initialize a struct containing useful powers of Mf at fRef */
  IMRPhenomX_UsefulPowers powers_of_MfRef;
  status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");
 
  /* Linear time and phase shifts so that model peaks near t ~ 0 */
  REAL8 lina = 0;
 
  /* Get phase connection coefficients */
  IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);
  double linb=IMRPhenomX_TimeShift_22(pPhase22, pWF);
 
  //
  INT4 APPLY_PNR_DEVIATIONS = pWF->APPLY_PNR_DEVIATIONS;
  INT4 PNRForceXHMAlignment = pWF->IMRPhenomXPNRForceXHMAlignment;
  // If applying PNR deviations, then we want to be able to refer to some non-PNR waveform properties. For that, we must compute the struct for when PNR is off (and specifically, XAS is wanted).
  if ( APPLY_PNR_DEVIATIONS && PNRForceXHMAlignment ) {
 
    /*<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.
    Shift phifRef and linb so that the PNR CoPrecessing model is aligned with XHM
    ->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.*/
 
    // // Development printing
    // printf("** Enforcing XAS Phase Alignment\n");
    // printf("(1) linb = %f\n",linb);
 
 
    IMRPhenomX_PNR_EnforceXASPhaseAlignment(&linb,pWF,pPhase22);
 
    // // Development printing
    // printf("(2) linb = %f\n",linb);
 
  }
  // extra contribution to phi(fRef) due to tidal corrections
  double phiTfRef = 0.;
 
  REAL8 f_final=freqs->data[freqs->length-1];
 
  // correct for time and phase shifts due to tidal phase
  if(NRTidal_version!=NoNRT_V){
 
      REAL8 f_merger;
      REAL8 f_merger_tmp;
      switch (NRTidal_version) {
          case NRTidalv3_V:
              f_merger_tmp = XLALSimNRTunedTidesMergerFrequency_v3(pWF->Mtot, pWF->lambda1, pWF->lambda2, pWF->q, pWF->chi1L, pWF->chi2L);
              break;
          default:
              f_merger_tmp = XLALSimNRTunedTidesMergerFrequency(pWF->Mtot, pWF->kappa2T, pWF->q);
              break;
      }
      f_merger = f_merger_tmp;
 
        if(f_merger<f_final)
            f_final = f_merger;
 
        IMRPhenomX_UsefulPowers powers_of_ffinal;
        REAL8 Mf_final = f_final*pWF->M_sec;
        status = IMRPhenomX_Initialize_Powers(&powers_of_ffinal,Mf_final);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for f_final.\n");
        REAL8 dphi_fmerger=1/pWF->eta*IMRPhenomX_dPhase_22(Mf_final, &powers_of_ffinal, pPhase22, pWF)+linb-IMRPhenomX_TidalPhaseDerivative(&powers_of_ffinal, pWF, pPhase22, NRTidal_version);
        REAL8 tshift = -dphi_fmerger; //This was adapted from the PhenomPv2 implementation; the resulting BBH limit can then have a time-shift in its phase
        linb+=tshift;
        phiTfRef = -IMRPhenomX_TidalPhase(&powers_of_MfRef, pWF, pPhase22, NRTidal_version);
 
    }
 
  /* 1/eta is used to re-scale phase */
  REAL8 inveta    = (1.0 / pWF->eta);
 
  /* Calculate phase at reference frequency: phifRef = 2.0*phi0 + LAL_PI_4 + PhenomXPhase(fRef) */
  double phifRef = -(inveta * IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) + phiTfRef + linb*pWF->MfRef + lina) + 2.0*pWF->phi0 + LAL_PI_4;
 
  // Define pWF->phifRef and phifRef separately, as "phifRef" may ultimately differ superficially
  pWF->phifRef = phifRef;
 
 
  /*
      Here we declare explicit REAL8 variables for main loop in order to avoid numerous
      pointer calls.
  */
  //REAL8 MfRef     = pWF->MfRef;
  REAL8 Msec      = pWF->M_sec;
 
  REAL8 C1IM      = pPhase22->C1Int;
  REAL8 C2IM      = pPhase22->C2Int;
  REAL8 C1RD      = pPhase22->C1MRD;
  REAL8 C2RD      = pPhase22->C2MRD;
 
  REAL8 fPhaseIN  = pPhase22->fPhaseMatchIN;
  REAL8 fPhaseIM  = pPhase22->fPhaseMatchIM;
  REAL8 fAmpIN    = pAmp22->fAmpMatchIN;
  REAL8 fAmpIM    = pAmp22->fAmpRDMin;
 
  if(debug)
  {
    printf("\n\n **** Phase struct initialized. **** \n\n");
    printf("C1IM     = %.4f\n",C1IM);
    printf("C2IM     = %.4f\n",C2IM);
    printf("C1RD     = %.4f\n",C1RD);
    printf("C2RD     = %.4f\n",C2RD);
    printf("fIN      = %.4f\n",fPhaseIN);
    printf("fIM      = %.4f\n",fPhaseIM);
  }
 
  REAL8 Amp0      = pWF->amp0 * pWF->ampNorm;
 
  /* initial_status used to track  */
  UINT4 initial_status = XLAL_SUCCESS;
 
  if (NRTidal_version!=NoNRT_V) {
    int ret = 0;
    UINT4 L_fCut = freqs->length;
    phi_tidal = XLALCreateREAL8Sequence(L_fCut);
    amp_tidal = XLALCreateREAL8Sequence(L_fCut);
    planck_taper = XLALCreateREAL8Sequence(L_fCut);
    /* Get FD tidal phase correction and amplitude factor */
    ret = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(phi_tidal, amp_tidal, planck_taper, freqs, pWF->m1_SI, pWF->m2_SI, lambda1, lambda2, pWF->chi1L, pWF->chi2L, NRTidal_version);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");
  }
 
  /* Now loop over main driver to generate waveform:  h(f) = A(f) * Exp[I phi(f)] */
  #pragma omp parallel for
  for (UINT4 idx = 0; idx < freqs->length; idx++)
  {
    double Mf    = Msec * freqs->data[idx];   // Mf is declared locally inside the loop
    UINT4 jdx    = idx  + offset;             // jdx is declared locally inside the loop
 
    /* We do not want to generate the waveform at frequencies > f_max (default = 0.3 Mf) */
    if(Mf <= (pWF->f_max_prime * pWF->M_sec))
    {
 
    /* Initialize a struct containing useful powers of Mf */
    IMRPhenomX_UsefulPowers powers_of_Mf;
    initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
    if(initial_status != XLAL_SUCCESS)
    {
      status = initial_status;
      XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
    }
    else
    {
      /* Generate amplitude and phase at MfRef */
      REAL8 amp = 0.0;
      REAL8 phi = 0.0;
 
      /* The functions in this routine are inlined to help performance. */
      /* Construct phase */
      if(Mf < fPhaseIN)
      {
        phi = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pPhase22);
      }
      else if(Mf > fPhaseIM)
      {
        phi = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1RD + (C2RD * Mf);
      }
      else
      {
        phi = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1IM + (C2IM * Mf);
      }
 
        /* Scale phase by 1/eta */
        phi  *= inveta;
        phi  += linb*Mf + lina + phifRef;
 
          /* Construct amplitude */
          if(Mf < fAmpIN)
          {
                  amp = IMRPhenomX_Inspiral_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
          }
          else if(Mf > fAmpIM)
          {
                  amp = IMRPhenomX_Ringdown_Amp_22_Ansatz(Mf, pWF, pAmp22);
          }
          else
          {
        amp = IMRPhenomX_Intermediate_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }
 
        // /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
        // ((*htilde22)->data->data)[jdx] = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
 
        /* NOTE that the above lines are commented out to clarify legacy code structure.
         HERE, we allow the user to toggle output of ONLY the model's phase using lalParams.
         The intended use of this option is to enable exact output of the phase, without unwanted numerical effects that result from unwrapping the complete waveform. In particular, at low frequencies, the phase may vary so quickly that adjacent frequency points correctly differ by more than 2*pi, thus limiting unwrap routines which assume that this is not the case.
        */
        if ( pWF->PhenomXOnlyReturnPhase ) {
          //
          ((*htilde22)->data->data)[jdx] = phi;
        } else {
        /* Add NRTidal phase, if selected, code adapted from LALSimIMRPhenomP.c */
 
      if (NRTidal_version!=NoNRT_V) {
 
          REAL8 phaseTidal = phi_tidal->data[idx];
          double ampTidal = amp_tidal->data[idx];
          double window = planck_taper->data[idx];
 
          /* Add spin-induced quadrupole moment terms to tidal phasing */
 
          /* 2PN terms */
          phaseTidal += pfaN * pPhase22->c2PN_tidal* powers_of_lalpi.m_one_third * powers_of_Mf.m_one_third;
 
          /* 3PN terms */
          phaseTidal += pfaN * pPhase22->c3PN_tidal* powers_of_lalpi.one_third * powers_of_Mf.one_third;
 
          /* 3.5PN terms are only in NRTidalv2 and NRTidalv3 */
          if (NRTidal_version == NRTidalv2_V || NRTidal_version == NRTidalv3_V) {
              phaseTidal += pfaN * pPhase22->c3p5PN_tidal * powers_of_lalpi.two_thirds * powers_of_Mf.two_thirds;
          }
            /* Reconstruct waveform with NRTidal terms included: h(f) = [A(f) + A_tidal(f)] * Exp{I [phi(f) - phi_tidal(f)]} * window(f) */
          ((*htilde22)->data->data)[jdx] = pWF->amp0 * (pWF->ampNorm * powers_of_Mf.m_seven_sixths * amp + 2*sqrt(1./5.)*powers_of_lalpi.sqrt * ampTidal) * cexp(I * (phi - phaseTidal))* window;
 
      }
      else if (NRTidal_version == NoNRT_V) {
        /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
        ((*htilde22)->data->data)[jdx] = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
        }
        else {
        XLAL_PRINT_INFO("Warning: Only NRTidal, NRTidalv2, NRTidalv3, and NoNRT NRTidal_version values allowed and NRTidal is not implemented completely in IMRPhenomX*.");
      }
    }
  }
    }
  else
 
    {
        /* Mf > Mf_max, so return 0 */
        ((*htilde22)->data->data)[jdx] = 0.0 + I*0.0;
    }
 
  }
 
  // Free allocated memory
  LALFree(pAmp22);
  LALFree(pPhase22);
  XLALDestroyREAL8Sequence(freqs);
 
  // Free allocated memory for tidal extension
  XLALDestroyREAL8Sequence(phi_tidal);
  XLALDestroyREAL8Sequence(amp_tidal);
  XLALDestroyREAL8Sequence(planck_taper);
 
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }
 
  return status;
}


/* ******** PRECESSING IMR PHENOMENOLOGICAL WAVEFORM: IMRPhenomXP ********* */
 
 
/*
    Decleration for internal function to generate precessing-spin, 22 only IMRPhenomXP waveform.
    Declared here and not in header due to IMRPhenomXPrecessionStruct.
*/
int IMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde,       /**< [out] FD waveform           */
  COMPLEX16FrequencySeries **hctilde,       /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,            /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Waveform Struct */
  LALDict *lalParams                        /**< LAL Dictionary Structure    */
);


 /* *********************************************************************************
  *
  * The following private function generates an IMRPhenomXP frequency-domain waveform
  *   - Precessing spins
  *   - Only the 22 mode in the co-precessing frame
  *   - Physical parameters are passed via the waveform struct
  * *********************************************************************************
  */
int IMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde,       /**< [out] FD waveform           */
  COMPLEX16FrequencySeries **hctilde,       /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,            /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Waveform Struct */
  LALDict *lalParams                        /**< LAL Dictionary Structure    */
)
{
  #if PHENOMXPDEBUG == 1
    printf("\n **** Now in IMRPhenomXPGenerateFD... **** \n");
  #endif
 
  /* Set tidal version */
  NRTidal_version_type NRTidal_version ;
  NRTidal_version=IMRPhenomX_SetTidalVersion(lalParams);
 
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}
 
  /* Initialize useful powers of LAL_PI */
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");
 
  /* Inherit minimum and maximum frequencies to generate wavefom from input frequency grid */
  double f_min = freqs_In->data[0];
  double f_max = freqs_In->data[freqs_In->length - 1];
 
  /* Size of array */
  size_t npts     = 0;
 
  /* Index shift between freqs and the frequency series */
  UINT4 offset    = 0;
 
  /* Initialize frequency sequence */
  REAL8Sequence *freqs = NULL;
 
  /* Initialize tidal corrections (following tidal code modified from LALSimIMRPhenomP.c) */
  REAL8Sequence *phi_tidal = NULL;
  REAL8Sequence *amp_tidal = NULL;
  REAL8Sequence *planck_taper = NULL;
 
  /* Set matter parameters (set to zero in pWF if NRTidal additions are not turned on) */
  REAL8 lambda1 = pWF->lambda1;
  REAL8 lambda2 = pWF->lambda2;
 
  /* New variables needed for the NRTidalv2 model */
  REAL8 X_A = pWF->m1; // Already scaled by Mtot
  REAL8 X_B = pWF->m2; // Ibid.
  REAL8 pfaN = 3./(128.*X_A*X_B);
 
  /* If deltaF is non-zero then we need to generate a uniformly sampled frequency grid of spacing deltaF. Start at f = 0. */
  if(pWF->deltaF > 0)
  {
    /* Return the closest power of 2 */
    npts = (size_t) (f_max / pWF->deltaF) + 1;
 
    /* Debug information */
    #if PHENOMXPDEBUG == 1
      printf("npts     = %zu\n",npts);
      printf("fMin     = %.4f\n",f_min);
      printf("fMax     = %.4f\n",f_max);
      printf("dF       = %.4f\n",pWF->deltaF);
    #endif
 
    /* Coalescence time is fixed to t=0, shift by overall length in time.  */
    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.\n", pWF->deltaF);
 
    /* Initialize the htilde frequency series */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
 
    /* Check that frequency series generated okay */
    XLAL_CHECK(*hptilde,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries h_+ of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);
    XLAL_CHECK(*hctilde,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries h_x of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);
 
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;
 
    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
          "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.\n", iStart, iStop, npts);
 
    /* Allocate memory for frequency array and terminate if this fails */
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
 
    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      freqs->data[i-iStart] = i * pWF->deltaF;
    }
    offset = iStart;
  }
  else
  {
    /* freqs is a frequency grid with non-uniform spacing, so we start at the lowest given frequency */
    npts      = freqs_In->length;
    *hptilde  = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);
    *hctilde  = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);
 
    XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries h_+ of length %zu from sequence.\n", npts);
    XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries h_x of length %zu from sequence.\n", npts);
 
    offset = 0;
    freqs  = XLALCreateREAL8Sequence(freqs_In->length);
 
    /* Allocate memory for frequency array and terminate if this fails */
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
 
    /* Populate frequency array */
    for (UINT4 i = 0; i < freqs_In->length; i++)
    {
      freqs->data[i] = freqs_In->data[i];
    }
  }
 
  memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hptilde)->sampleUnits), &((*hptilde)->sampleUnits), &lalSecondUnit);
 
  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);
 
  /* Check if LAL dictionary exists. If not, create a LAL dictionary. */
  UNUSED INT4 lalParams_In = 0;
  if(lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams    = XLALCreateDict();
  }
 
  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Initializing amplitude struct... **** \n\n");
  #endif
 
    // numerical angles;
    REAL8Sequence *alpha = NULL;
    REAL8Sequence *gamma = NULL;
    REAL8Sequence *cosbeta = NULL;
 
    if(pPrec->precessing_tag==3){
 
        status = IMRPhenomXPSpinTaylorAnglesIMR(&alpha,&cosbeta,&gamma,freqs,pWF,pPrec,lalParams);
        XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPSpinTaylorAnglesIMR failed.");
 
    }
 
 
  /* Allocate and initialize the PhenomX 22 amplitude coefficients struct */
  IMRPhenomXAmpCoefficients *pAmp22;
  pAmp22 = XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  status = IMRPhenomXGetAmplitudeCoefficients(pWF,pAmp22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAmplitudeCoefficients failed.\n");
 
  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Amplitude struct initialized. **** \n\n");
    printf("\n\n **** Initializing phase struct... **** \n\n");
  #endif
 
  /* Allocate and initialize the PhenomX 22 phase coefficients struct */
  IMRPhenomXPhaseCoefficients *pPhase22;
  pPhase22 = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhase22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");
    if(NRTidal_version!=NoNRT_V) IMRPhenomXGetTidalPhaseCoefficients(pWF,pPhase22,NRTidal_version);
 
  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Phase struct initialized. **** \n\n");
  #endif
 
  /* Initialize a struct containing useful powers of Mf at fRef */
  IMRPhenomX_UsefulPowers powers_of_MfRef;
  status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");
 
  /*
      Time shift so peak amplitude is near t ~ 0.
  */
 
  /* Linear time and phase shifts so that model peaks near t ~ 0                        */
  /* lina provides an additional phase shift, which we currently don't use, but we keep */
  /* the variable in case we change our mind about phase alignement in the future       */
  REAL8 lina = 0;
 
  /* Get phase connection coefficients = phase offsets between calibration regions to achieve C^1 phase */
 
  IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);
 
  /* Apply time shift, see IMRPhenomX_TimeShift_22 for details of current implementation */
  double linb = IMRPhenomX_TimeShift_22(pPhase22, pWF);
   // extra contribution to phi(fRef) due to tidal corrections
  double phiTfRef = 0.;
 
  // correct for time and phase shifts due to tidal phase
  REAL8 f_final=freqs->data[freqs->length-1];
 
  if(NRTidal_version!=NoNRT_V){
      REAL8 f_merger;
      REAL8 f_merger_tmp;
      switch (NRTidal_version) {
          case NRTidalv3_V:
              f_merger_tmp = XLALSimNRTunedTidesMergerFrequency_v3(pWF->Mtot, pWF->lambda1, pWF->lambda2, pWF->q, pWF->chi1L, pWF->chi2L);
              break;
          default:
              f_merger_tmp = XLALSimNRTunedTidesMergerFrequency(pWF->Mtot, pWF->kappa2T, pWF->q);
              break;
      }
 
      f_merger = f_merger_tmp;
      if(f_merger<f_final)
          f_final = f_merger;
 
      IMRPhenomX_UsefulPowers powers_of_ffinal;
      REAL8 Mf_final = f_final*pWF->M_sec;
      status = IMRPhenomX_Initialize_Powers(&powers_of_ffinal,Mf_final);
      XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for f_final.\n");
      REAL8 dphi_fmerger=1/pWF->eta*IMRPhenomX_dPhase_22(Mf_final, &powers_of_ffinal, pPhase22, pWF)+linb-IMRPhenomX_TidalPhaseDerivative(&powers_of_ffinal, pWF, pPhase22, NRTidal_version);
      REAL8 tshift = -dphi_fmerger; //This was adapted from the PhenomPv2 implementation; the resulting BBH limit can then have a time-shift in its phase
      linb+=tshift;
      // tidal phase will be subtracted from the BBH phase
      phiTfRef = -IMRPhenomX_TidalPhase(&powers_of_MfRef, pWF, pPhase22, NRTidal_version);
  }
 
 
  /* Inverse of the symmetric mass ratio */
  REAL8 inveta    = (1.0 / pWF->eta);
 
  /* Construct reference phase, see discussion in Appendix A of arXiv:2001.11412 */
  pWF->phifRef    = -(inveta * IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) + linb*pWF->MfRef + lina+phiTfRef) + 2.0*pWF->phi0 + LAL_PI_4;
 
  /*
      Here we declare explicit REAL8 variables for main loop in order to avoid numerous
      pointer calls.
  */
  REAL8 C1IM              = pPhase22->C1Int;
  REAL8 C2IM              = pPhase22->C2Int;
  REAL8 C1RD              = pPhase22->C1MRD;
  REAL8 C2RD              = pPhase22->C2MRD;
 
  REAL8 fPhaseIN          = pPhase22->fPhaseMatchIN;
  REAL8 fPhaseIM          = pPhase22->fPhaseMatchIM;
  REAL8 fAmpIN            = pAmp22->fAmpMatchIN;
  REAL8 fAmpIM            = pAmp22->fAmpRDMin;
 
 /* add in PNR-specific code */
 
  int PNRUseTunedAngles     = pPrec->IMRPhenomXPNRUseTunedAngles;
 
  REAL8Sequence* alphaPNR = NULL;
  REAL8Sequence* betaPNR = NULL;
  REAL8Sequence* gammaPNR = NULL;
 
  if(PNRUseTunedAngles) /* Check for PNR angles */
  {
    alphaPNR = XLALCreateREAL8Sequence(freqs->length);
    if (!alphaPNR)
    {
      XLAL_ERROR(XLAL_EFUNC, "alphaPNR array allocation failed in LALSimIMRPhenomX.c");
    }
    betaPNR = XLALCreateREAL8Sequence(freqs->length);
    if (!betaPNR)
    {
      XLAL_ERROR(XLAL_EFUNC, "betaPNR array allocation failed in LALSimIMRPhenomX.c");
    }
    gammaPNR = XLALCreateREAL8Sequence(freqs->length);
    if (!gammaPNR)
    {
      XLAL_ERROR(XLAL_EFUNC, "gammaPNR array allocation failed in LALSimIMRPhenomX.c");
    }
 
    status = IMRPhenomX_PNR_GeneratePNRAngles(
      alphaPNR, betaPNR, gammaPNR,
      freqs, pWF, pPrec,
      lalParams);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_GeneratePNRAngles failed.\n");
 
    #if DEBUG == 1
      // Save angles into a file
      FILE *fileangle2;
      char fileSpec2[40];
      sprintf(fileSpec2, "angles_PNR.dat");
 
      printf("\nOutput angle file: %s\r\n", fileSpec2);
      fileangle2 = fopen(fileSpec2,"w");
 
      fprintf(fileangle2,"# q = %.16e m1 = %.16e m2 = %.16e chi1 = %.16e chi2 = %.16e  Mtot = %.16e distance = %.16e\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L,  pWF->Mtot, pWF->distance/LAL_PC_SI/1e6);
      fprintf(fileangle2,"#fHz   alpha   beta    gamma\n");
 
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        fprintf(fileangle2, "%.16e  %.16e  %.16e  %.16e\n", freqs->data[idx], alphaPNR->data[idx], betaPNR->data[idx], gammaPNR->data[idx]);
      }
      fclose(fileangle2);
    #endif
  }
 
  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Phase struct initialized. **** \n\n");
    printf("C1IM     = %.4f\n",C1IM);
    printf("C2IM     = %.4f\n",C2IM);
    printf("C1RD     = %.4f\n",C1RD);
    printf("C2RD     = %.4f\n",C2RD);
    printf("fIN      = %.4f\n",fPhaseIN);
    printf("fIM      = %.4f\n",fPhaseIM);
    printf("thetaJN  = %.4f\n",pPrec->thetaJN);
  #endif
 
  /* amplitude definition as in XAS */
 
  REAL8 Amp0          = pWF->amp0 * pWF->ampNorm;
 
 
 
  /* initial_status used to track  */
  UINT4 initial_status = XLAL_SUCCESS;
 
  /* Approximate frequency of the peak */
  REAL8 fmax         = pAmp22->fAmpRDMin;
 
  /* Initialize a struct containing useful powers of fmax */
  IMRPhenomX_UsefulPowers powers_of_fmax;
  initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_fmax,fmax);
  if(initial_status != XLAL_SUCCESS)
  {
    status = initial_status;
    XLALPrintError("IMRPhenomX_Initialize_Powers failed for fmax, initial_status=%d",initial_status);
  }
 
  #if DEBUG == 1
  // Save anles into a file
  FILE *fileangle;
  char fileSpec[40];
  sprintf(fileSpec, "angles_XP.dat");
 
  printf("\nOutput angle file: %s\r\n", fileSpec);
  fileangle = fopen(fileSpec,"w");
 
  fprintf(fileangle,"# q = %.16e m1 = %.16e m2 = %.16e chi1 = %.16e chi2 = %.16e  Mtot = %.16e distance = %.16e\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L,  pWF->Mtot, pWF->distance/LAL_PC_SI/1e6);
  fprintf(fileangle,"#fHz   cexp_i_alpha(re im)   cexp_i_epsilon(re im)    cexp_i_betah(re im)\n");
 
  fclose(fileangle);
  #endif
 
  if (NRTidal_version == NRTidal_V || NRTidal_version == NRTidalv2_V || NRTidal_version == NRTidalv3_V) {
    int ret = 0;
    UINT4 L_fCut = freqs->length;
    phi_tidal = XLALCreateREAL8Sequence(L_fCut);
    amp_tidal = XLALCreateREAL8Sequence(L_fCut);
    planck_taper = XLALCreateREAL8Sequence(L_fCut);
    /* Get FD tidal phase correction and amplitude factor */
    ret = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(phi_tidal, amp_tidal, planck_taper, freqs, pWF->m1_SI, pWF->m2_SI, lambda1, lambda2, pWF->chi1L, pWF->chi2L, NRTidal_version);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");
  }
 
 
  int AntisymmetricWaveform = pPrec->IMRPhenomXAntisymmetricWaveform;
 
  /** declare all variables used in anti-symmetric waveform calculation*/
  REAL8Sequence *kappa = NULL;
  REAL8 A0 = 0.0;
  REAL8 phi_A0 = 0.0;
  REAL8 phi_B0 = 0.0;
 
  REAL8 phi_antiSym = 0.0;
  REAL8 amp_antiSym = 0.0;
  double MfT = 0.85 * pWF->fRING; /* note that fRING is already in dimensionaless units  */
 
  if(AntisymmetricWaveform) /* Check for antisymmetric waveform flag */
  {
    if (!PNRUseTunedAngles)
    {
      XLAL_ERROR(XLAL_EFUNC, "Error: Antisymmetric waveform generation not supported without PNR angles, please turn on PNR angles to produce waveform with asymmetries in the (2,2) and (2,-2) modes\n");
    }
    else
    {
      kappa = XLALCreateREAL8Sequence(freqs->length);
      if (!kappa)
      {
        XLAL_ERROR(XLAL_EFUNC, "antisymmetric waveform amplitude ratio array allocation failed in LALSimIMRPhenomX.c");
      }
      status = IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(kappa, freqs, pWF, pPrec);
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNRAntisymmetricAmpRatio failed.\n");
 
      status = IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients(
      &A0, &phi_A0, &phi_B0, MfT, lina, linb, inveta, pWF,pPrec,pPhase22);
 
    }
  }
 
  /* Now loop over frequencies to generate waveform:  h(f) = A(f) * Exp[I phi(f)] */
  #pragma omp parallel for
  for (UINT4 idx = 0; idx < freqs->length; idx++)
  {
    double Mf    = pWF->M_sec * freqs->data[idx];
    UINT4 jdx    = idx  + offset;
 
    COMPLEX16 hcoprec     = 0.0;  /* Co-precessing waveform */
    COMPLEX16 hcoprec_antiSym  = 0.0;  /* Co-precessing anti-symmetric waveform */
    COMPLEX16 hplus       = 0.0;  /* h_+ */
    COMPLEX16 hcross      = 0.0;  /* h_x */
 
    /* We do not want to generate the waveform at frequencies > f_max (default = 0.3 Mf) */
    if(Mf <= (pWF->f_max_prime * pWF->M_sec))
    {
    /* Initialize a struct containing useful powers of Mf */
    IMRPhenomX_UsefulPowers powers_of_Mf;
    initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
    if(initial_status != XLAL_SUCCESS)
    {
      status = initial_status;
      XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d\n",initial_status);
    }
    else
    {
      /* Generate amplitude and phase at Mf */
 
      // initialize amplitude and phase
      REAL8 amp = 0.0;
      REAL8 phi = 0.0;
 
      /* Here we explicitly call the functions which treat the non-overlapping */
      /* inspiral, intermediate and ringdown frequency regions for the non-precessing waveform. */
 
      /* Get phase */
      if(Mf < fPhaseIN)
      {
        phi = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pPhase22);
      }
      else if(Mf > fPhaseIM)
      {
        phi = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1RD + (C2RD * Mf);
      }
      else
      {
        phi = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1IM + (C2IM * Mf);
      }
      /* Scale phase by 1/eta and apply phase and time shifts */
      phi  *= inveta;
      phi  += linb*Mf + lina + pWF->phifRef;
 
      /* Get amplitude */
      if(Mf < fAmpIN)
      {
        amp = IMRPhenomX_Inspiral_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }
      else if(Mf > fAmpIM)
      {
        amp = IMRPhenomX_Ringdown_Amp_22_Ansatz(Mf, pWF, pAmp22);
      }
      else
      {
        amp = IMRPhenomX_Intermediate_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }
 
      /* Add NRTidal phase, if selected, code adapted from LALSimIMRPhenomP.c; the following is unused in generating IMRPhenomXP_NRTidalv2/3 waveforms using GeneratorLegacy.c */
 
      if (NRTidal_version == NRTidal_V || NRTidal_version == NRTidalv2_V || NRTidal_version == NRTidalv3_V) {
 
          REAL8 phaseTidal = phi_tidal->data[idx];
          double ampTidal = amp_tidal->data[idx];
          double window = planck_taper->data[idx];
 
          /* Add spin-induced quadrupole moment terms to tidal phasing */
 
          /* 2PN terms */
          phaseTidal += pfaN * pPhase22->c2PN_tidal* powers_of_lalpi.m_one_third * powers_of_Mf.m_one_third;
          /* 3PN terms */
          phaseTidal += pfaN * pPhase22->c3PN_tidal* powers_of_lalpi.one_third* powers_of_Mf.one_third;
 
          /* 3.5PN terms are only in NRTidalv2 */
          if (NRTidal_version == NRTidalv2_V || NRTidal_version == NRTidalv3_V) {
              phaseTidal += pfaN * pPhase22->c3p5PN_tidal * powers_of_lalpi.two_thirds * powers_of_Mf.two_thirds;
          }
 
          /* Waveform in co-precessing frame with NRTidal terms included: h(f) = [A(f) + A_tidal(f)] * Exp{I [phi(f) - phi_tidal(f)]} * window(f) */
          hcoprec = pWF->amp0 * (pWF->ampNorm * powers_of_Mf.m_seven_sixths * amp + 2*sqrt(1/5.)*powers_of_lalpi.sqrt * ampTidal) * cexp(I * (phi - phaseTidal)) * window;
        }
 
      else if (NRTidal_version == NoNRT_V) {
          /* Waveform in co-precessing frame: h(f) = A(f) * Exp[I phi(f)] */
        hcoprec = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
      }
 
      else {
            XLAL_PRINT_INFO("Warning: Only NRTidal, NRTidalv2, NRTidalv3, and NoNRT NRTidal_version values allowed and NRTidal is not implemented completely in IMRPhenomX*.");
            }
 
        /* Transform modes from co-precessing frame to inertial frame */
        /* Only do this if the coprecessing model is not desired */
        if(  pWF->IMRPhenomXReturnCoPrec  )
        {
          //
          hplus  =  0.5 * (hcoprec);
          hcross = -0.5 * I * (hcoprec);
 
        }
        else
        {
 
          if(PNRUseTunedAngles) /* Look for PNR flag */
          {
            pPrec->alphaPNR = alphaPNR->data[idx];
            pPrec->betaPNR = betaPNR->data[idx];
            pPrec->gammaPNR = gammaPNR->data[idx];
          }
 
 
        if(pPrec->precessing_tag==3)
           IMRPhenomXPTwistUp22_NumericalAngles(hcoprec, alpha->data[idx], cosbeta->data[idx], gamma->data[idx], pPrec, &hplus, &hcross);
       else
           IMRPhenomXPTwistUp22(Mf,hcoprec,pWF,pPrec,&hplus,&hcross);
 
          /********************************/
          /*** anti-symmetric waveform ***/
          /*******************************/
 
          if(AntisymmetricWaveform && PNRUseTunedAngles)
          {
            /* phase of anti-symmetric waveform*/
            if(Mf < MfT)
            {
              phi_antiSym = phi/2 + alphaPNR->data[idx] + A0 *Mf + phi_A0;
            }
            else
            {
              phi_antiSym = phi + phi_B0;
            }
 
            COMPLEX16 hplus_antiSym       = 0.0;
            COMPLEX16 hcross_antiSym      = 0.0;
            amp_antiSym = cabs(hcoprec)*kappa->data[idx];
            hcoprec_antiSym = amp_antiSym* cexp(I * phi_antiSym);
            pPrec->PolarizationSymmetry = -1;
            IMRPhenomXPTwistUp22(Mf,hcoprec_antiSym,pWF,pPrec,&hplus_antiSym, &hcross_antiSym);
            pPrec->PolarizationSymmetry = 1;
            hplus += hplus_antiSym;
            hcross += hcross_antiSym;
          }
        }
 
 
      /* Populate h_+ and h_x */
      ((*hptilde)->data->data)[jdx] = hplus;
      ((*hctilde)->data->data)[jdx] = hcross;
    }
  }
    else
    {
      /* Mf > Mf_max, so return 0 */
      ((*hptilde)->data->data)[jdx] = 0.0 + I*0.0;
      ((*hctilde)->data->data)[jdx] = 0.0 + I*0.0;
    }
  }
 
XLALDestroyREAL8Sequence(alpha);
XLALDestroyREAL8Sequence(cosbeta);
XLALDestroyREAL8Sequence(gamma);
 
 
    if(pPrec->precessing_tag==3)
    {
        LALFree(pPrec->alpha_params);
        LALFree(pPrec->beta_params);
 
        gsl_spline_free(pPrec->alpha_spline);
        gsl_spline_free(pPrec->cosbeta_spline);
        gsl_spline_free(pPrec->gamma_spline);
 
        gsl_interp_accel_free(pPrec->alpha_acc);
        gsl_interp_accel_free(pPrec->gamma_acc);
        gsl_interp_accel_free(pPrec->cosbeta_acc);
 
    }
 
 
  /*
      Loop over h+ and hx and rotate waveform by 2 \zeta.
 
      Note that we do this here rather than in LALSimInpsiral.c as
      \zeta_polarization is part of the pPrec struct.
 
  */
  COMPLEX16 PhPpolp, PhPpolc;
  REAL8 cosPolFac, sinPolFac;
 
  /* If \zeta is non-zero, we need to rotate h+ and hx */
  if(fabs(pPrec->zeta_polarization) > 0.0)
  {
    cosPolFac = cos(2.0 * pPrec->zeta_polarization);
    sinPolFac = sin(2.0 * pPrec->zeta_polarization);
 
    for (UINT4 i = 0; i < (*hptilde)->data->length; i++)
    {
        PhPpolp = (*hptilde)->data->data[i];
        PhPpolc = (*hctilde)->data->data[i];
 
        (*hptilde)->data->data[i] = (cosPolFac * PhPpolp) + (sinPolFac * PhPpolc);
        (*hctilde)->data->data[i] = (cosPolFac * PhPpolc) - (sinPolFac * PhPpolp);
    }
  }
 
  /* Free allocated memory */
  LALFree(pAmp22);
  LALFree(pPhase22);
  XLALDestroyREAL8Sequence(freqs);
 
  // Free allocated memory for tidal extension
  XLALDestroyREAL8Sequence(phi_tidal);
  XLALDestroyREAL8Sequence(amp_tidal);
  XLALDestroyREAL8Sequence(planck_taper);
 
  if(PNRUseTunedAngles){
    XLALDestroyREAL8Sequence(alphaPNR);
    XLALDestroyREAL8Sequence(betaPNR);
    XLALDestroyREAL8Sequence(gammaPNR);
  }
 
  if(AntisymmetricWaveform){
    XLALDestroyREAL8Sequence(kappa);
  }
 
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }
 
  return status;
}