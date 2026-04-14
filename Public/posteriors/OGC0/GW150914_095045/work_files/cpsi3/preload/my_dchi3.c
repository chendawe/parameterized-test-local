// /***********************************************************************************************************************************************
// LALSimIMRPhenomX.c 的注释；


//  **********************************************************************************************************************************************/

/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
// #include <lal/LALDict.c>


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
// #include "LALSimIMRPhenomX_inspiral.h"
// #include "LALSimIMRPhenomX_internals.c"
#include "LALSimIMRPhenomX_precession.c"
#include "LALSimIMRPhenomX_PNR.c"
#include "LALSimIMRPhenomX_AntisymmetricWaveform.c"
 
/* Note: This is declared in LALSimIMRPhenomX_internals.c and avoids namespace clashes */
IMRPhenomX_UsefulPowers powers_of_lalpi;
 
#ifndef _OPENMP
#define omp ignore
#endif


#include "LALSimIMRPhenomXUtilities.h"

// #include "LALSimIMRPhenomX_qnm.c"
// #include "LALSimIMRPhenomX_PNR_deviations.c"
// #include <lal/LALConstants.h>
// #include <lal/Date.h>
// #include <lal/Units.h>
 
/* GSL Header Files */
#include <gsl/gsl_linalg.h>

// /* Link PN coefficients (needed for NRTidal spin-induced quadrupole terms) */
// // #include "LALSimInspiralPNCoefficients.c"



/*
 * Function to populate the IMRPhenomXPhaseCoefficients struct:
*/
int IMRPhenomXGetPhaseCoefficients(
        IMRPhenomXWaveformStruct *pWF,
        IMRPhenomXPhaseCoefficients *pPhase
)
{
        /* Get LALparams */
        LALDict *LALparams    = pWF->LALparams;
        const INT4 debug      = PHENOMXDEBUG;
 
        /* GSL objects for solving system of equations via LU decomposition */
        gsl_vector *b, *x;
        gsl_matrix *A;
        gsl_permutation *p;
        int s;
 
        REAL8 deltax;
        REAL8 xmin;
        REAL8 fi;
 
        REAL8 gpoints4[4]     = {0.0, 1.0/4.0, 3.0/4.0, 1.0};
        REAL8 gpoints5[5]     = {0.0, 1.0/2 - 1.0/(2*sqrt(2.0)), 1.0/2.0, 1.0/2 + 1.0/(2.0*sqrt(2.0)), 1.0};
 
        // Matching regions
 
        /* This is Eq. 5.11 in the paper */
        double fIMmatch = 0.6 * (0.5 * pWF->fRING + pWF->fISCO);
 
        /* This is the MECO frequency */
        double fINmatch = pWF->fMECO;
 
        /* This is Eq. 5.10 in the paper */
        double deltaf   = (fIMmatch - fINmatch) * 0.03;
 
        // Transition frequency is just below the MECO frequency and just above the RD fitting region
 
        /* These are defined in Eq. 7.7 and the text just below, f_H = fPhaseMatchIM and f_L = fPhaseMatchIN */
        pPhase->fPhaseMatchIN  = fINmatch - 1.0*deltaf;
        pPhase->fPhaseMatchIM  = fIMmatch + 0.5*deltaf;
 
        /* Defined in Eq. 7.4, this is f_L */
        pPhase->fPhaseInsMin  = 0.0026;
 
        /* Defined in Eq. 7.4, this is f_H */
        pPhase->fPhaseInsMax  = 1.020 * pWF->fMECO;
 
        /* Defined in Eq. 7.12, this is f_L */
        pPhase->fPhaseRDMin   = fIMmatch;
 
        /* Defined in Eq. 7.12, this is f_L */
        pPhase->fPhaseRDMax   = pWF->fRING + 1.25*pWF->fDAMP;
 
        pPhase->phiNorm                 = -(3.0 * powers_of_lalpi.m_five_thirds) / (128.0);
 
        /* For convenience, define some variables here */
        REAL8 chi1L           = pWF->chi1L;
        REAL8 chi2L           = pWF->chi2L;
 
        REAL8 chi1L2L         = chi1L * chi2L;
 
        REAL8 chi1L2          = pWF->chi1L * pWF->chi1L;
        REAL8 chi1L3          = pWF->chi1L * chi1L2;
 
        REAL8 chi2L2          = pWF->chi2L * pWF->chi2L;
        REAL8 chi2L3          = pWF->chi2L * chi2L2;
 
        REAL8 eta             = pWF->eta;
        REAL8 eta2            = eta*eta;
        REAL8 eta3            = eta*eta2;
 
        REAL8 delta           = pWF->delta;
 
        /* Pre-initialize all phenomenological coefficients */
        pPhase->a0 = 0.0;
        pPhase->a1 = 0.0;
        pPhase->a2 = 0.0;
        pPhase->a3 = 0.0;
        pPhase->a4 = 0.0;
 
        pPhase->b0 = 0.0;
        pPhase->b1 = 0.0;
        pPhase->b2 = 0.0;
        pPhase->b3 = 0.0;
        pPhase->b4 = 0.0;
 
        pPhase->c0 = 0.0;
        pPhase->c1 = 0.0;
        pPhase->c2 = 0.0;
        pPhase->c3 = 0.0;
        pPhase->c4 = 0.0;
        pPhase->cL = 0.0;
 
    pPhase->c2PN_tidal   = 0.;
    pPhase->c3PN_tidal   = 0.;
    pPhase->c3p5PN_tidal = 0.;
 
        pPhase->sigma0 = 0.0;
        pPhase->sigma1 = 0.0;
        pPhase->sigma2 = 0.0;
        pPhase->sigma3 = 0.0;
        pPhase->sigma4 = 0.0;
        pPhase->sigma5 = 0.0;
 
 
        /*
                The general strategy is to initialize a linear system of equations:
 
                A.x = b
 
                - A is a matrix with the coefficients of the ansatz evaluated at the collocation nodes.
                - b is a vector of the value of the collocation points
                - x is the solution vector (i.e. the coefficients) that we must solve for.
 
                We choose to do this using a standard LU decomposition.
        */
 
        /* Generate list of collocation points */
        /*
        The Gauss-Chebyshev Points are given by:
        GCPoints[n] = Table[Cos[i * pi/n] , {i, 0, n}]
 
        SamplePoints[xmin,xmax,n] = {
                pWF->delta = xmax - xmin;
                gpoints = 0.5 * (1 + GCPoints[n-1]);
 
                return {xmin + pWF->delta * gpoints}
        }
 
        gpoints4 = [1.0, 3.0/4, 1.0/4, 0.0];
        gpoints5 = [1.0, 1.0/2 + 1.0/(2.0*sqrt(2.0)), 1.0/2, 1.0/2 - 1.0/(2*sqrt(2.0)), 0.]
 
        */
 
        /*
        Ringdown phase collocation points:
        Here we only use the first N+1 points in the array where N = the
        number of pseudo PN terms.
 
        The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_RD
 
        Default is to use 5 collocation points.
        */
        deltax = pPhase->fPhaseRDMax - pPhase->fPhaseRDMin;
        xmin   = pPhase->fPhaseRDMin;
        int i;
 
        double phaseRD; // This is used in intermediate phase reconstruction.
 
        if(debug)
        {
                printf("\n");
                printf("Solving system of equations for RD phase...\n");
        }
 
        // Initialize collocation points
        for(i = 0; i < 5; i++)
        {
                pPhase->CollocationPointsPhaseRD[i] = gpoints5[i] * deltax + xmin;
        }
        // Collocation point 4 is set to the ringdown frequency ~ dip in Lorentzian
        pPhase->CollocationPointsPhaseRD[3] = pWF->fRING;
 
        if(debug)
        {
                printf("Rigndown collocation points : \n");
                printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseRD[0]);
                printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseRD[1]);
                printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseRD[2]);
                printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseRD[3]);
                printf("F5 : %.6f\n",pPhase->CollocationPointsPhaseRD[4]);
        }
 
        switch(pWF->IMRPhenomXRingdownPhaseVersion)
        {
                case 105:
                {
                        pPhase->NCollocationPointsRD = 5;
                        break;
                }
                default:
                {
                        XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXRingdownPhaseVersion is not valid.\n");
                }
        }
 
        if(debug)
        {
                printf("NCollRD = %d\n",pPhase->NCollocationPointsRD);
        }
 
        // Eq. 7.13 in arXiv:2001.11412
        double RDv4 = IMRPhenomX_Ringdown_Phase_22_v4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
 
        /* These are the calibrated collocation points, as per Eq. 7.13 */
        pPhase->CollocationValuesPhaseRD[0] = IMRPhenomX_Ringdown_Phase_22_d12(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
        pPhase->CollocationValuesPhaseRD[1] = IMRPhenomX_Ringdown_Phase_22_d24(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
        pPhase->CollocationValuesPhaseRD[2] = IMRPhenomX_Ringdown_Phase_22_d34( pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
        pPhase->CollocationValuesPhaseRD[3] = RDv4;
        pPhase->CollocationValuesPhaseRD[4] = IMRPhenomX_Ringdown_Phase_22_d54(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
 
        /* v_j = d_{j4} + v4 */
        pPhase->CollocationValuesPhaseRD[4] = pPhase->CollocationValuesPhaseRD[4] + pPhase->CollocationValuesPhaseRD[3]; // v5 = d54  + v4
        pPhase->CollocationValuesPhaseRD[2] = pPhase->CollocationValuesPhaseRD[2] + pPhase->CollocationValuesPhaseRD[3]; // v3 = d34  + v4
        pPhase->CollocationValuesPhaseRD[1] = pPhase->CollocationValuesPhaseRD[1] + pPhase->CollocationValuesPhaseRD[3]; // v2 = d24  + v4
        pPhase->CollocationValuesPhaseRD[0] = pPhase->CollocationValuesPhaseRD[0] + pPhase->CollocationValuesPhaseRD[1]; // v1 = d12  + v2
 
        // Debugging information. Leave for convenience later on.
        if(debug)
        {
                printf("\n");
                printf("Ringdown Collocation Points: \n");
                printf("v1 : %.6f\n",pPhase->CollocationValuesPhaseRD[0]);
                printf("v2 : %.6f\n",pPhase->CollocationValuesPhaseRD[1]);
                printf("v3 : %.6f\n",pPhase->CollocationValuesPhaseRD[2]);
                printf("v4 : %.6f\n",pPhase->CollocationValuesPhaseRD[3]);
                printf("v5 : %.6f\n",pPhase->CollocationValuesPhaseRD[4]);
                printf("\n");
        }
 
        phaseRD = pPhase->CollocationValuesPhaseRD[0];
 
        p = gsl_permutation_alloc(pPhase->NCollocationPointsRD);
        b = gsl_vector_alloc(pPhase->NCollocationPointsRD);
        x = gsl_vector_alloc(pPhase->NCollocationPointsRD);
        A = gsl_matrix_alloc(pPhase->NCollocationPointsRD,pPhase->NCollocationPointsRD);
 
        /*
        Populate the b vector
        */
        gsl_vector_set(b,0,pPhase->CollocationValuesPhaseRD[0]);
        gsl_vector_set(b,1,pPhase->CollocationValuesPhaseRD[1]);
        gsl_vector_set(b,2,pPhase->CollocationValuesPhaseRD[2]);
        gsl_vector_set(b,3,pPhase->CollocationValuesPhaseRD[3]);
        gsl_vector_set(b,4,pPhase->CollocationValuesPhaseRD[4]);
 
        /*
                        Eq. 7.12 in arXiv:2001.11412
 
                        ansatzRD(f) = a_0 + a_1 f^(-1/3) + a_2 f^(-2) + a_3 f^(-3) + a_4 f^(-4) + ( aRD ) / ( (f_damp^2 + (f - f_ring)^2 ) )
 
                        Canonical ansatz sets a_3 to 0.
        */
 
        /*
                We now set up and solve a linear system of equations.
                First we populate the matrix A_{ij}
        */
 
        /* Note that ff0 is always 1 */
        REAL8 ff, invff, ff0, ff1, ff2, ff3, ff4;
 
        /* A_{0,i} */
        ff    = pPhase->CollocationPointsPhaseRD[0];
        invff = 1.0 / ff;
        ff1   = cbrt(invff);   // f^{-1/3}
        ff2   = invff * invff; // f^{-2}
        ff3   = ff2 * ff2;               // f^{-4}
        ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
        gsl_matrix_set(A,0,0,1.0); // Constant
        gsl_matrix_set(A,0,1,ff1); // f^(-1/3) term
        gsl_matrix_set(A,0,2,ff2); // f^(-2) term
        gsl_matrix_set(A,0,3,ff3); // f^(-4) term
        gsl_matrix_set(A,0,4,ff4); // Lorentzian term
 
        if(debug)
        {
                printf("For row 0: a0 + a1 %.6f + a2 %.6f + a4 %.6f + aRD %.6f\n",ff1,ff2,ff3,ff4);
        }
 
        /* A_{1,i} */
        ff    = pPhase->CollocationPointsPhaseRD[1];
        invff = 1.0 / ff;
        ff1   = cbrt(invff);
        ff2   = invff * invff;
        ff3   = ff2 * ff2;
        ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
        gsl_matrix_set(A,1,0,1.0);
        gsl_matrix_set(A,1,1,ff1);
        gsl_matrix_set(A,1,2,ff2);
        gsl_matrix_set(A,1,3,ff3);
        gsl_matrix_set(A,1,4,ff4);
 
        /* A_{2,i} */
        ff    =  pPhase->CollocationPointsPhaseRD[2];
        invff = 1.0 / ff;
        ff1   = cbrt(invff);
        ff2   = invff * invff;
        ff3   = ff2 * ff2;
        ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
        gsl_matrix_set(A,2,0,1.0);
        gsl_matrix_set(A,2,1,ff1);
        gsl_matrix_set(A,2,2,ff2);
        gsl_matrix_set(A,2,3,ff3);
        gsl_matrix_set(A,2,4,ff4);
 
        /* A_{3,i} */
        ff    = pPhase->CollocationPointsPhaseRD[3];
        invff = 1.0 / ff;
        ff1   = cbrt(invff);
        ff2   = invff * invff;
        ff3   = ff2 * ff2;
        ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
        gsl_matrix_set(A,3,0,1.0);
        gsl_matrix_set(A,3,1,ff1);
        gsl_matrix_set(A,3,2,ff2);
        gsl_matrix_set(A,3,3,ff3);
        gsl_matrix_set(A,3,4,ff4);
 
        /* A_{4,i} */
        ff    = pPhase->CollocationPointsPhaseRD[4];
        invff = 1.0 / ff;
        ff1   = cbrt(invff);
        ff2   = invff * invff;
        ff3   = ff2 * ff2;
        ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
        gsl_matrix_set(A,4,0,1.0);
        gsl_matrix_set(A,4,1,ff1);
        gsl_matrix_set(A,4,2,ff2);
        gsl_matrix_set(A,4,3,ff3);
        gsl_matrix_set(A,4,4,ff4);
 
        /* We now solve the system A x = b via an LU decomposition */
        gsl_linalg_LU_decomp(A,p,&s);
        gsl_linalg_LU_solve(A,p,b,x);
 
        pPhase->c0  = gsl_vector_get(x,0); // x[0];     // a0
        pPhase->c1  = gsl_vector_get(x,1); // x[1];             // a1
        pPhase->c2  = gsl_vector_get(x,2); // x[2];     // a2
        pPhase->c4  = gsl_vector_get(x,3); // x[3];     // a4
        pPhase->cRD = gsl_vector_get(x,4);
        pPhase->cL  = -(pWF->dphase0 * pPhase->cRD); // ~ x[4] // cL = - a_{RD} * dphase0
 
        /* Apply NR tuning for precessing cases (500) */
        pPhase->cL = pPhase->cL + ( pWF->PNR_DEV_PARAMETER * pWF->NU4 );
        // pPhase->c0 = pPhase->c0 + ( pWF->PNR_DEV_PARAMETER * pWF->NU0 );
 
        if(debug)
        {
                printf("\n");
                printf("Ringdown Coefficients: \n");
                printf("c0  : %.6f\n",pPhase->c0);
                printf("c1  : %.6f\n",pPhase->c1);
                printf("c2  : %.6f\n",pPhase->c2);
                printf("c4  : %e\n",pPhase->c4);
                printf("cRD : %.6f\n",gsl_vector_get(x,4));
                printf("d0  : %.6f\n",pWF->dphase0);
                printf("cL  : %e\n",pPhase->cL);
                printf("\n");
 
                printf("Freeing arrays...\n");
        }
 
        /* Tidy up in preparation for next GSL solve ... */
        gsl_vector_free(b);
        gsl_vector_free(x);
        gsl_matrix_free(A);
        gsl_permutation_free(p);
 
        /*
        Inspiral phase collocation points:
        Here we only use the first N+1 points in the array where N = the
        number of pseudo PN terms. E.g. for 4 pseudo-PN terms, we will
        need 4 collocation points. An ansatz with n free coefficients
        needs n pieces of information in order to constrain the ansatz.
 
        The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_INS
 
        Default is to use 4 pseudo-PN coefficients and hence 4 collocation points.
 
        GC points as per Eq. 7.4 and 7.5, where f_L = pPhase->fPhaseInsMin and f_H = pPhase->fPhaseInsMax
        */
        deltax      = pPhase->fPhaseInsMax - pPhase->fPhaseInsMin;
        xmin        = pPhase->fPhaseInsMin;
 
        /*
                        Set number of pseudo-PN coefficients:
                                - If you add a new PN inspiral approximant, update with new version here.
        */
        switch(pWF->IMRPhenomXInspiralPhaseVersion)
        {
                case 104:
                {
                        pPhase->NPseudoPN = 4;
                        pPhase->NCollocationPointsPhaseIns = 4;
                        break;
                }
                case 105:
                {
                        pPhase->NPseudoPN = 5;
                        pPhase->NCollocationPointsPhaseIns = 5;
                        break;
                }
                case 114:
                {
                        pPhase->NPseudoPN = 4;
                        pPhase->NCollocationPointsPhaseIns = 4;
                        break;
                }
                case 115:
                {
                        pPhase->NPseudoPN = 5;
                        pPhase->NCollocationPointsPhaseIns = 5;
                        break;
                }
                default:
                {
                        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error: IMRPhenomXInspiralPhaseVersion is not valid.\n");
                }
        }
 
        if(debug)
        {
                printf("\n");
                printf("NPseudoPN : %d\n",pPhase->NPseudoPN);
                printf("NColl : %d\n",pPhase->NCollocationPointsPhaseIns);
                printf("\n");
        }
 
        p = gsl_permutation_alloc(pPhase->NCollocationPointsPhaseIns);
        b = gsl_vector_alloc(pPhase->NCollocationPointsPhaseIns);
        x = gsl_vector_alloc(pPhase->NCollocationPointsPhaseIns);
        A = gsl_matrix_alloc(pPhase->NCollocationPointsPhaseIns,pPhase->NCollocationPointsPhaseIns);
 
        /*
        If we are using 4 pseudo-PN coefficients, call the routines below.
        The inspiral phase version is still passed to the individual functions.
        */
        if(pPhase->NPseudoPN == 4)
        {
                // By default all models implemented use the following GC points.
                // If a new model is calibrated with different choice of collocation points, edit this.
                for(i = 0; i < pPhase->NCollocationPointsPhaseIns; i++)
                {
                        fi = gpoints4[i] * deltax + xmin;
                        pPhase->CollocationPointsPhaseIns[i] = fi;
                }
 
                // Calculate the value of the differences between the ith and 3rd collocation points at the GC nodes
                pPhase->CollocationValuesPhaseIns[0] = IMRPhenomX_Inspiral_Phase_22_d13(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[1] = IMRPhenomX_Inspiral_Phase_22_d23(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[2] = IMRPhenomX_Inspiral_Phase_22_v3( pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[3] = IMRPhenomX_Inspiral_Phase_22_d43(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
 
                // Calculate the value of the collocation points at GC nodes via: v_i = d_i3 + v3
                pPhase->CollocationValuesPhaseIns[0] = pPhase->CollocationValuesPhaseIns[0] + pPhase->CollocationValuesPhaseIns[2];
                pPhase->CollocationValuesPhaseIns[1] = pPhase->CollocationValuesPhaseIns[1] + pPhase->CollocationValuesPhaseIns[2];
                pPhase->CollocationValuesPhaseIns[3] = pPhase->CollocationValuesPhaseIns[3] + pPhase->CollocationValuesPhaseIns[2];
 
                if(debug)
                {
                        printf("\n");
                        printf("Inspiral Collocation Points and Values:\n");
                        printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseIns[0]);
                        printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseIns[1]);
                        printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseIns[2]);
                        printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseIns[3]);
                        printf("\n");
                        printf("V1 : %.6f\n",pPhase->CollocationValuesPhaseIns[0]);
                        printf("V2 : %.6f\n",pPhase->CollocationValuesPhaseIns[1]);
                        printf("V3 : %.6f\n",pPhase->CollocationValuesPhaseIns[2]);
                        printf("V4 : %.6f\n",pPhase->CollocationValuesPhaseIns[3]);
                        printf("\n");
                }
 
                gsl_vector_set(b,0,pPhase->CollocationValuesPhaseIns[0]);
                gsl_vector_set(b,1,pPhase->CollocationValuesPhaseIns[1]);
                gsl_vector_set(b,2,pPhase->CollocationValuesPhaseIns[2]);
                gsl_vector_set(b,3,pPhase->CollocationValuesPhaseIns[3]);
 
                /* A_{0,i} */
                ff  = pPhase->CollocationPointsPhaseIns[0];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff0 = 1.0;
                gsl_matrix_set(A,0,0,1.0);
                gsl_matrix_set(A,0,1,ff1);
                gsl_matrix_set(A,0,2,ff2);
                gsl_matrix_set(A,0,3,ff3);
 
                /* A_{1,i} */
                ff  = pPhase->CollocationPointsPhaseIns[1];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff0 = 1.0;
                gsl_matrix_set(A,1,0,1.0);
                gsl_matrix_set(A,1,1,ff1);
                gsl_matrix_set(A,1,2,ff2);
                gsl_matrix_set(A,1,3,ff3);
 
                /* A_{2,i} */
                ff  = pPhase->CollocationPointsPhaseIns[2];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff0 = 1.0;
                gsl_matrix_set(A,2,0,1.0);
                gsl_matrix_set(A,2,1,ff1);
                gsl_matrix_set(A,2,2,ff2);
                gsl_matrix_set(A,2,3,ff3);
 
                /* A_{3,i} */
                ff  = pPhase->CollocationPointsPhaseIns[3];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff0 = 1.0;
                gsl_matrix_set(A,3,0,1.0);
                gsl_matrix_set(A,3,1,ff1);
                gsl_matrix_set(A,3,2,ff2);
                gsl_matrix_set(A,3,3,ff3);
 
                /* We now solve the system A x = b via an LU decomposition */
                gsl_linalg_LU_decomp(A,p,&s);
                gsl_linalg_LU_solve(A,p,b,x);
 
                /* Set inspiral phenomenological coefficients from solution to A x = b */
                pPhase->a0 = gsl_vector_get(x,0); // x[0]; // alpha_0
                pPhase->a1 = gsl_vector_get(x,1); // x[1]; // alpha_1
                pPhase->a2 = gsl_vector_get(x,2); // x[2]; // alpha_2
                pPhase->a3 = gsl_vector_get(x,3); // x[3]; // alpha_3
                pPhase->a4 = 0.0;
 
                /*
                                PSEUDO PN TERMS WORK:
                                        - 104 works.
                                        - 105 not tested.
                                        - 114 not tested.
                                        - 115 not tested.
                */
                if(debug)
                {
                        printf("\n");
                        printf("3pPN\n");
                        printf("Inspiral Pseudo-PN Coefficients:\n");
                        printf("a0 : %.6f\n",pPhase->a0);
                        printf("a1 : %.6f\n",pPhase->a1);
                        printf("a2 : %.6f\n",pPhase->a2);
                        printf("a3 : %.6f\n",pPhase->a3);
                        printf("a4 : %.6f\n",pPhase->a4);
                        printf("\n");
                }
 
                /* Tidy up in preparation for next GSL solve ... */
                gsl_vector_free(b);
                gsl_vector_free(x);
                gsl_matrix_free(A);
                gsl_permutation_free(p);
 
        }
        else if(pPhase->NPseudoPN == 5)
        {
                // Using 5 pseudo-PN coefficients so set 5 collocation points
                for(i = 0; i < 5; i++)
                {
                        fi = gpoints5[i] * deltax + xmin;
                        pPhase->CollocationPointsPhaseIns[i] = fi;
                }
                pPhase->CollocationValuesPhaseIns[0] = IMRPhenomX_Inspiral_Phase_22_d13(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[1] = IMRPhenomX_Inspiral_Phase_22_d23(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[2] = IMRPhenomX_Inspiral_Phase_22_v3( pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[3] = IMRPhenomX_Inspiral_Phase_22_d43(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
                pPhase->CollocationValuesPhaseIns[4] = IMRPhenomX_Inspiral_Phase_22_d53(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
 
                /* v_j = d_j3 + v_3 */
                pPhase->CollocationValuesPhaseIns[0] = pPhase->CollocationValuesPhaseIns[0] + pPhase->CollocationValuesPhaseIns[2];
                pPhase->CollocationValuesPhaseIns[1] = pPhase->CollocationValuesPhaseIns[1] + pPhase->CollocationValuesPhaseIns[2];
                pPhase->CollocationValuesPhaseIns[3] = pPhase->CollocationValuesPhaseIns[3] + pPhase->CollocationValuesPhaseIns[2];
                pPhase->CollocationValuesPhaseIns[4] = pPhase->CollocationValuesPhaseIns[4] + pPhase->CollocationValuesPhaseIns[2];
 
                if(debug)
                {
                        printf("\n");
                        printf("Inspiral Collocation Points and Values:\n");
                        printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseIns[0]);
                        printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseIns[1]);
                        printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseIns[2]);
                        printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseIns[3]);
                        printf("F5 : %.6f\n",pPhase->CollocationPointsPhaseIns[4]);
                        printf("\n");
                        printf("V1 : %.6f\n",pPhase->CollocationValuesPhaseIns[0]);
                        printf("V2 : %.6f\n",pPhase->CollocationValuesPhaseIns[1]);
                        printf("V3 : %.6f\n",pPhase->CollocationValuesPhaseIns[2]);
                        printf("V4 : %.6f\n",pPhase->CollocationValuesPhaseIns[3]);
                        printf("V5 : %.6f\n",pPhase->CollocationValuesPhaseIns[4]);
                        printf("\n");
                }
 
                gsl_vector_set(b,0,pPhase->CollocationValuesPhaseIns[0]);
                gsl_vector_set(b,1,pPhase->CollocationValuesPhaseIns[1]);
                gsl_vector_set(b,2,pPhase->CollocationValuesPhaseIns[2]);
                gsl_vector_set(b,3,pPhase->CollocationValuesPhaseIns[3]);
                gsl_vector_set(b,4,pPhase->CollocationValuesPhaseIns[4]);
 
                /* A_{0,i} */
                ff  = pPhase->CollocationPointsPhaseIns[0];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff4 = ff * ff1;
                gsl_matrix_set(A,0,0,1.0);
                gsl_matrix_set(A,0,1,ff1);
                gsl_matrix_set(A,0,2,ff2);
                gsl_matrix_set(A,0,3,ff3);
                gsl_matrix_set(A,0,4,ff4);
 
                /* A_{1,i} */
                ff  = pPhase->CollocationPointsPhaseIns[1];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff4 = ff * ff1;
                gsl_matrix_set(A,1,0,1.0);
                gsl_matrix_set(A,1,1,ff1);
                gsl_matrix_set(A,1,2,ff2);
                gsl_matrix_set(A,1,3,ff3);
                gsl_matrix_set(A,1,4,ff4);
 
                /* A_{2,i} */
                ff  = pPhase->CollocationPointsPhaseIns[2];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff4 = ff * ff1;
                gsl_matrix_set(A,2,0,1.0);
                gsl_matrix_set(A,2,1,ff1);
                gsl_matrix_set(A,2,2,ff2);
                gsl_matrix_set(A,2,3,ff3);
                gsl_matrix_set(A,2,4,ff4);
 
                /* A_{3,i} */
                ff  = pPhase->CollocationPointsPhaseIns[3];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff4 = ff * ff1;
                gsl_matrix_set(A,3,0,1.0);
                gsl_matrix_set(A,3,1,ff1);
                gsl_matrix_set(A,3,2,ff2);
                gsl_matrix_set(A,3,3,ff3);
                gsl_matrix_set(A,3,4,ff4);
 
                /* A_{4,i} */
                ff  = pPhase->CollocationPointsPhaseIns[4];
                ff1 = cbrt(ff);
                ff2 = ff1 * ff1;
                ff3 = ff;
                ff4 = ff * ff1;
                gsl_matrix_set(A,4,0,1.0);
                gsl_matrix_set(A,4,1,ff1);
                gsl_matrix_set(A,4,2,ff2);
                gsl_matrix_set(A,4,3,ff3);
                gsl_matrix_set(A,4,4,ff4);
 
                /* We now solve the system A x = b via an LU decomposition */
                gsl_linalg_LU_decomp(A,p,&s);
                gsl_linalg_LU_solve(A,p,b,x);
 
                /* Set inspiral phenomenological coefficients from solution to A x = b */
                pPhase->a0 = gsl_vector_get(x,0); // x[0];
                pPhase->a1 = gsl_vector_get(x,1); // x[1];
                pPhase->a2 = gsl_vector_get(x,2); // x[2];
                pPhase->a3 = gsl_vector_get(x,3); // x[3];
                pPhase->a4 = gsl_vector_get(x,4); // x[4];
 
                if(debug)
                {
                        printf("\n");
                        printf("4pPN\n");
                        printf("Inspiral Pseudo-PN Coefficients:\n");
                        printf("a0 : %.6f\n",pPhase->a0);
                        printf("a1 : %.6f\n",pPhase->a1);
                        printf("a2 : %.6f\n",pPhase->a2);
                        printf("a3 : %.6f\n",pPhase->a3);
                        printf("a4 : %.6f\n",pPhase->a4);
                        printf("\n");
                }
 
                /* Tidy up in preparation for next GSL solve ... */
                gsl_vector_free(b);
                gsl_vector_free(x);
                gsl_matrix_free(A);
                gsl_permutation_free(p);
        }
        else
        {
                XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: NPseudoPN requested is not valid.\n");
        }
 
        /* The pseudo-PN coefficients are normalized such that: (dphase0 / eta) * f^{8/3} * a_j */
        /* So we must re-scale these terms by an extra factor of f^{-8/3} in the PN phasing */
        pPhase->sigma1 = (-5.0/3.0) * pPhase->a0;
        pPhase->sigma2 = (-5.0/4.0) * pPhase->a1;
        pPhase->sigma3 = (-5.0/5.0) * pPhase->a2;
        pPhase->sigma4 = (-5.0/6.0) * pPhase->a3;
        pPhase->sigma5 = (-5.0/7.0) * pPhase->a4;
 
        /* Initialize TaylorF2 PN coefficients  */
        pPhase->dphi0  = 0.0;
        pPhase->dphi1  = 0.0;
        pPhase->dphi2  = 0.0;
        pPhase->dphi3  = 0.0;
        pPhase->dphi4  = 0.0;
        pPhase->dphi5  = 0.0;
        pPhase->dphi6  = 0.0;
        pPhase->dphi7  = 0.0;
        pPhase->dphi8  = 0.0;
        pPhase->dphi9  = 0.0;
        pPhase->dphi10 = 0.0;
        pPhase->dphi11 = 0.0;
        pPhase->dphi12 = 0.0;
 
        pPhase->dphi5L = 0.0;
        pPhase->dphi6L = 0.0;
        pPhase->dphi8L = 0.0;
        pPhase->dphi9L = 0.0;
 
        pPhase->phi0   = 0.0;
        pPhase->phi1   = 0.0;
        pPhase->phi2   = 0.0;
        pPhase->phi3   = 0.0;
        pPhase->phi4   = 0.0;
        pPhase->phi5   = 0.0;
        pPhase->phi6   = 0.0;
        pPhase->phi7   = 0.0;
        pPhase->phi8   = 0.0;
        pPhase->phi9   = 0.0;
        pPhase->phi10  = 0.0;
        pPhase->phi11  = 0.0;
        pPhase->phi12  = 0.0;
 
        pPhase->phi5L  = 0.0;
        pPhase->phi6L  = 0.0;
        pPhase->phi8L  = 0.0;
        pPhase->phi9L  = 0.0;
 
        /* **** TaylorF2 PN Coefficients: Phase **** */
 
        /*
                        - These are the PN coefficients normalised by: 3 / (128 * eta * [pi M f]^{5/3} ).
                        - We add in powers of (M f)^{N/3} later but add powers of pi^{N/3} here
                        - The log terms are *always* in terms of log(v), so we multiply by log(v) when summing PN phasing series.
                        - We *do not* overwrite the PN phasing series with pseudo-PN terms. These are added separately.
 
                        PN terms can be found in:
                                - Marsat et al, CQG, 32, 085008, (2015)
                                - Bohe et al, CQG, 32, 195010, (2015)
                                - Bernard et al, PRD, 95, 044026, (2017)
                                - Bernard et al, PRD, 93, 084037, (2016)
                                - Damour et al, PRD, 89, 064058, (2014)
                                - Damour et al, PRD, 95, 084005, (2017)
                                - Bernard et al, PRD, 96, 104043, (2017)
                                - Marchand et al, PRD, 97, 044023, (2018)
                                - Marchand et al, CQG, 33, 244003, (2016)
                                - Nagar et al, PRD, 99, 044007, (2019)
                                - Messina et al, PRD, 97, 084016, (2018)
        */
 
        /* Split into non-spinning and spin-dependent coefficients */
        UNUSED REAL8 phi0S = 0.0, phi1S = 0.0, phi2S = 0.0;
        REAL8 phi0NS  = 0.0,  phi1NS = 0.0,  phi2NS = 0.0;
        REAL8 phi3NS  = 0.0,  phi3S  = 0.0,  phi4NS = 0.0,  phi4S   = 0.0,  phi5NS  = 0.0,  phi5S  = 0.0;
        REAL8 phi5LNS = 0.0,  phi5LS = 0.0,  phi6NS = 0.0,  phi6S   = 0.0,  phi6LNS = 0.0,  phi6LS = 0.0;
        REAL8 phi7NS  = 0.0,  phi7S  = 0.0,  phi8NS = 0.0,  phi8S   = 0.0,  phi8LNS = 0.0;
        REAL8 phi8LS  = 0.0,  phi9NS = 0.0,  phi9S  = 0.0,  phi9LNS = 0.0,  phi9LS  = 0.0;
 
 
        /* Analytically known PN coefficients */
        /* Newtonian */
        phi0NS         = 1.0;
 
        /* ~~ 0.5 PN ~~ */
        phi1NS         = 0.0;
 
        /* ~~ 1.0 PN ~~ */
        /* 1.0PN, Non-Spinning */
        phi2NS         = (3715/756. + (55*eta)/9.) * powers_of_lalpi.two_thirds;
 
        /* ~~ 1.5 PN ~~ */
        /* 1.5PN, Non-Spinning */
        phi3NS         = -16.0 * powers_of_lalpi.two;
        /* 1.5PN, Spin-Orbit */
        phi3S          = ( (113*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 76*(chi1L + chi2L)*eta)/6. ) * powers_of_lalpi.itself;
 
        /* ~~ 2.0 PN ~~ */
        /* 2.0PN, Non-Spinning */
        phi4NS         = ( 15293365/508032. + (27145*eta)/504. + (3085*eta2)/72. ) * powers_of_lalpi.four_thirds;
        /* 2.0PN, Spin-Spin */
        phi4S          = ( (-5*(81*chi1L2*(1 + delta - 2*eta) + 316*chi1L2L*eta - 81*chi2L2*(-1 + delta + 2*eta)))/16. ) * powers_of_lalpi.four_thirds;
 
        /* ~~ 2.5 PN ~~ */
        phi5NS         = 0.0;
        phi5S          = 0.0;
 
        /* ~~ 2.5 PN, Log Term ~~ */
        /* 2.5PN, Non-Spinning */
        phi5LNS        = ( (5*(46374 - 6552*eta)*LAL_PI)/4536. ) * powers_of_lalpi.five_thirds;
        /* 2.5PN, Spin-Orbit */
        phi5LS         = ( (-732985*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 560*(-1213*(chi1L + chi2L)
                                                                                                + 63*(chi1L - chi2L)*delta)*eta + 85680*(chi1L + chi2L)*eta2)/4536. ) * powers_of_lalpi.five_thirds;
 
        /* ~~ 3.0 PN ~~ */
        /* 3.0 PN, Non-Spinning */
        phi6NS         = ( 11583231236531/4.69421568e9 - (5*eta*(3147553127 + 588*eta*(-45633 + 102260*eta)))/3.048192e6 - (6848*LAL_GAMMA)/21.
                                                - (640*powers_of_lalpi.two)/3. + (2255*eta*powers_of_lalpi.two)/12. - (13696*log(2))/21. - (6848*powers_of_lalpi.log)/63. ) * powers_of_lalpi.two;
        /* 3.0 PN, Spin-Orbit */
        phi6S          = ( (5*(227*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 156*(chi1L + chi2L)*eta)*LAL_PI)/3. ) * powers_of_lalpi.two;
        /* 3.0 PN, Spin-Spin */
        phi6S         += ( (5*(20*chi1L2L*eta*(11763 + 12488*eta) + 7*chi2L2*(-15103*(-1 + delta) + 2*(-21683 + 6580*delta)*eta - 9808*eta2) -
                                                        7*chi1L2*(-15103*(1 + delta) + 2*(21683 + 6580*delta)*eta + 9808*eta2)))/4032. ) * powers_of_lalpi.two;
 
        /* ~~ 3.0 PN, Log Term ~~ */
        phi6LNS        = (-6848/63.) * powers_of_lalpi.two;
        phi6LS         = 0.0;
 
        /* ~~ 3.5 PN ~~ */
        /* 3.5 PN, Non-Spinning */
        phi7NS         = ( (5*(15419335 + 168*(75703 - 29618*eta)*eta)*LAL_PI)/254016. ) * powers_of_lalpi.seven_thirds;
        /* 3.5 PN, Spin-Orbit */
        phi7S          = ( (5*(-5030016755*(chi1L + chi2L + chi1L*delta - chi2L*delta) + 4*(2113331119*(chi1L + chi2L) + 675484362*(chi1L - chi2L)*delta)*eta - 1008*(208433*(chi1L + chi2L) + 25011*(chi1L - chi2L)*delta)*eta2 + 90514368*(chi1L + chi2L)*eta3))/6.096384e6 ) * powers_of_lalpi.seven_thirds;
        /* 3.5 PN, Spin-Spin */
        phi7S         += ( -5*(57*chi1L2*(1 + delta - 2*eta) + 220*chi1L2L*eta - 57*chi2L2*(-1 + delta + 2*eta))*LAL_PI ) * powers_of_lalpi.seven_thirds;
        /* 3.5 PN, Cubic-in-Spin */
        phi7S         += ( (14585*(-(chi2L3*(-1 + delta)) + chi1L3*(1 + delta)) - 5*(chi2L3*(8819 - 2985*delta) + 8439*chi1L*chi2L2*(-1 + delta) - 8439*chi1L2*chi2L*(1 + delta) + chi1L3*(8819 + 2985*delta))*eta + 40*(chi1L + chi2L)*(17*chi1L2 - 14*chi1L2L + 17*chi2L2)*eta2)/48. ) * powers_of_lalpi.seven_thirds;
 
                /* ~~ 4.0 PN ~~ */
        /* 4.0 PN, Non-Spinning */
        phi8NS         = 0.0;
        /* 4.0 PN, Spin-Orbit */
        phi8S          = ( (-5*(1263141*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 2*(794075*(chi1L + chi2L) + 178533*(chi1L - chi2L)*delta)*eta + 94344*(chi1L + chi2L)*eta2)*LAL_PI*(-1 + powers_of_lalpi.log))/9072. ) * powers_of_lalpi.eight_thirds;
 
        /* ~~ 4.0 PN, Log Term ~~ */
        /* 4.0 PN, log term, Non-Spinning */
        phi8LNS        = 0.0;
        /* 4.0 PN, log term, Spin-Orbit */
        phi8LS         = ((-5*(1263141*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 2*(794075*(chi1L + chi2L) + 178533*(chi1L - chi2L)*delta)*eta
                                                        + 94344*(chi1L + chi2L)*eta2)*LAL_PI)/9072.) * powers_of_lalpi.eight_thirds;
 
        /* ~~ 4.5 PN ~~ */
        phi9NS         = 0.0;
        phi9S          = 0.0;
 
                /* ~~ 4.5 PN, Log Term ~~ */
        phi9LNS        = 0.0;
        phi9LS         = 0.0;
 
        /* This version of TaylorF2 contains an additional 4.5PN tail term and a LO-SS tail term at 3.5PN */
        if(pWF->IMRPhenomXInspiralPhaseVersion == 114 || pWF->IMRPhenomXInspiralPhaseVersion == 115)
        {
                        /* 3.5PN, Leading Order Spin-Spin Tail Term */
                        phi7S         += ( (5*(65*chi1L2*(1 + delta - 2*eta) + 252*chi1L2L*eta - 65*chi2L2*(-1 + delta + 2*eta))*LAL_PI)/4. ) * powers_of_lalpi.seven_thirds;
 
                        /* 4.5PN, Tail Term */
                        phi9NS        += ( (5*(-256 + 451*eta)*powers_of_lalpi.three)/6. + (LAL_PI*(105344279473163 + 700*eta*(-298583452147 + 96*eta*(99645337 + 14453257*eta)) -
                                                                                                                                                                                12246091038720*LAL_GAMMA - 24492182077440*log(2.0)))/1.877686272e10 - (13696*LAL_PI*powers_of_lalpi.log)/63. ) * powers_of_lalpi.three;
 
                        /* 4.5PN, Log Term */
                        phi9LNS       += (  (-13696*LAL_PI)/63.0  ) * powers_of_lalpi.three;
        }
 
        /* 0.0 PN */
        pPhase->phi0   = phi0NS;
 
        /* 0.5 PN */
        pPhase->phi1   = phi1NS;
 
        /* 1.0 PN */
        pPhase->phi2   = phi2NS;
 
        /* 1.5 PN */
        pPhase->phi3   = phi3NS + phi3S;
 
        /* 2.0 PN */
        pPhase->phi4   = phi4NS + phi4S;
 
        /* 2.5 PN */
        pPhase->phi5   = phi5NS + phi5S;
 
        /* 2.5 PN, Log Terms */
        pPhase->phi5L  = phi5LNS + phi5LS;
 
        /* 3.0 PN */
        pPhase->phi6   = phi6NS + phi6S;
 
        /* 3.0 PN, Log Term */
        pPhase->phi6L  = phi6LNS + phi6LS;
 
        /* 3.5PN */
        pPhase->phi7   = phi7NS + phi7S;
 
        /* 4.0PN */
        pPhase->phi8   = phi8NS + phi8S;
 
        /* 4.0 PN, Log Terms */
        pPhase->phi8L  = phi8LNS + phi8LS;
 
        /* 4.5 PN */
        pPhase->phi9   = phi9NS + phi9S;
 
        /* 4.5 PN, Log Terms */
        pPhase->phi9L  = phi9LNS + phi9LS;
 
 
        if(debug)
        {
                printf("TaylorF2 PN Coefficients: \n");
                printf("phi0   : %.6f\n",pPhase->phi0);
                printf("phi1   : %.6f\n",pPhase->phi1);
                printf("phi2   : %.6f\n",pPhase->phi2);
                printf("phi3   : %.6f\n",pPhase->phi3);
                printf("phi4   : %.6f\n",pPhase->phi4);
                printf("phi5   : %.6f\n",pPhase->phi5);
                printf("phi6   : %.6f\n",pPhase->phi6);
                printf("phi7   : %.6f\n",pPhase->phi7);
                printf("phi8   : %.6f\n",pPhase->phi8);
 
                printf("phi5L  : %.6f\n",pPhase->phi5L);
                printf("phi6L  : %.6f\n",pPhase->phi6L);
                printf("phi8L  : %.6f\n",pPhase->phi8L);
 
                printf("phi8P  : %.6f\n",pPhase->sigma1);
                printf("phi9P  : %.6f\n",pPhase->sigma2);
                printf("phi10P : %.6f\n",pPhase->sigma3);
                printf("phi11P : %.6f\n",pPhase->sigma4);
                printf("phi12P : %.6f\n",pPhase->sigma5);
        }
 
        pPhase->phi_initial = - LAL_PI_4;
 
        /* **** TaylorF2 PN Coefficients: Normalized Phase Derivative **** */
        pPhase->dphi0  = pPhase->phi0;
        pPhase->dphi1  = 4.0 / 5.0 * pPhase->phi1;
        pPhase->dphi2  = 3.0 / 5.0 * pPhase->phi2;
        pPhase->dphi3  = 2.0 / 5.0 * pPhase->phi3;
        pPhase->dphi4  = 1.0 / 5.0 * pPhase->phi4;
        pPhase->dphi5  = -3.0 / 5.0 * pPhase->phi5L;
        pPhase->dphi6  = -1.0 / 5.0 * pPhase->phi6 - 3.0 / 5.0 * pPhase->phi6L;
        pPhase->dphi6L = -1.0 / 5.0 * pPhase->phi6L;
        pPhase->dphi7  = -2.0 / 5.0 * pPhase->phi7;
        pPhase->dphi8  = -3.0 / 5.0 * pPhase->phi8 - 3.0 / 5.0 * pPhase->phi8L;
        pPhase->dphi8L = -3.0 / 5.0 * pPhase->phi8L;
        pPhase->dphi9  = -4.0 / 5.0 * pPhase->phi9 - 3.0 / 5.0 * pPhase->phi9L;
        pPhase->dphi9L = -3.0 / 5.0 * pPhase->phi9L;
 
        if(debug)
        {
                printf("\nTaylorF2 PN Derivative Coefficients\n");
                printf("dphi0  : %.6f\n",pPhase->dphi0);
                printf("dphi1  : %.6f\n",pPhase->dphi1);
                printf("dphi2  : %.6f\n",pPhase->dphi2);
                printf("dphi3  : %.6f\n",pPhase->dphi3);
                printf("dphi4  : %.6f\n",pPhase->dphi4);
                printf("dphi5  : %.6f\n",pPhase->dphi5);
                printf("dphi6  : %.6f\n",pPhase->dphi6);
                printf("dphi7  : %.6f\n",pPhase->dphi7);
                printf("dphi8  : %.6f\n",pPhase->dphi8);
                printf("dphi9  : %.6f\n",pPhase->dphi9);
                printf("\n");
                printf("dphi6L : %.6f\n",pPhase->dphi6L);
                printf("dphi8L : %.6f\n",pPhase->dphi8L);
                printf("dphi9L : %.6f\n",pPhase->dphi9L);
        }
 
        /*
                        Calculate phase at fmatchIN. This will be used as the collocation point for the intermediate fit.
                        In practice, the transition point is just below the MECO frequency.
        */
        if(debug)
        {
                printf("\nTransition frequency for ins to int : %.6f\n",pPhase->fPhaseMatchIN);
        }
 
        IMRPhenomX_UsefulPowers powers_of_fmatchIN;
        IMRPhenomX_Initialize_Powers(&powers_of_fmatchIN,pPhase->fPhaseMatchIN);
 
        double phaseIN;
        phaseIN  = pPhase->dphi0;                                                                                                                                       // f^{0/3}
        phaseIN += pPhase->dphi1        * powers_of_fmatchIN.one_third;                                                                 // f^{1/3}
        phaseIN += pPhase->dphi2        * powers_of_fmatchIN.two_thirds;                                                                // f^{2/3}
        phaseIN += pPhase->dphi3        * powers_of_fmatchIN.itself;                                                                    // f^{3/3}
        phaseIN += pPhase->dphi4        * powers_of_fmatchIN.four_thirds;                                                               // f^{4/3}
        phaseIN += pPhase->dphi5        * powers_of_fmatchIN.five_thirds;                                                               // f^{5/3}
        phaseIN += pPhase->dphi6        * powers_of_fmatchIN.two;                                                                               // f^{6/3}
        phaseIN += pPhase->dphi6L       * powers_of_fmatchIN.two * powers_of_fmatchIN.log;                              // f^{6/3}, Log[f]
        phaseIN += pPhase->dphi7        * powers_of_fmatchIN.seven_thirds;                                                              // f^{7/3}
        phaseIN += pPhase->dphi8        * powers_of_fmatchIN.eight_thirds;                                                              // f^{8/3}
        phaseIN += pPhase->dphi8L       * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.log;             // f^{8/3}
        phaseIN += pPhase->dphi9        * powers_of_fmatchIN.three;                                                                             // f^{9/3}
        phaseIN += pPhase->dphi9L       * powers_of_fmatchIN.three * powers_of_fmatchIN.log;                    // f^{9/3}
 
        // Add pseudo-PN Coefficient
        phaseIN += (            pPhase->a0 * powers_of_fmatchIN.eight_thirds
                                                                + pPhase->a1 * powers_of_fmatchIN.three
                                                                + pPhase->a2 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.two_thirds
                                                                + pPhase->a3 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.itself
                                                                + pPhase->a4 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.four_thirds
                                                        );
 
        phaseIN  = phaseIN * powers_of_fmatchIN.m_eight_thirds * pWF->dphase0;
 
        /*
        Intermediate phase collocation points:
        Here we only use the first N points in the array where N = the
        number of intermediate collocation points.
 
        The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_INT
 
        Default is to use 5 collocation points.
 
        See. Eq. 7.7 and 7.8 where f_H = pPhase->fPhaseMatchIM and f_L = pPhase->fPhaseMatchIN
        */
        deltax      = pPhase->fPhaseMatchIM - pPhase->fPhaseMatchIN;
        xmin        = pPhase->fPhaseMatchIN;
 
        switch(pWF->IMRPhenomXIntermediatePhaseVersion)
        {
                case 104:
                {
                        // Fourth order polynomial ansatz
                        pPhase->NCollocationPointsInt = 4;
                        break;
                }
                case 105:
                {
                        // Fifth order polynomial ansatz
                        pPhase->NCollocationPointsInt = 5;
                        break;
                }
                default:
                {
                        XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
                }
        }
 
        if(debug)
        {
        printf("\nNColPointsInt : %d\n",pPhase->NCollocationPointsInt);
        }
 
        p = gsl_permutation_alloc(pPhase->NCollocationPointsInt);
        b = gsl_vector_alloc(pPhase->NCollocationPointsInt);
        x = gsl_vector_alloc(pPhase->NCollocationPointsInt);
        A = gsl_matrix_alloc(pPhase->NCollocationPointsInt,pPhase->NCollocationPointsInt);
 
        // Canonical intermediate model using 4 collocation points
        if(pWF->IMRPhenomXIntermediatePhaseVersion == 104)
        {
                // Using 4 collocation points in intermediate region
                for(i = 0; i < 4; i++)
                {
                        fi = gpoints4[i] * deltax + xmin;
 
                        pPhase->CollocationPointsPhaseInt[i] = fi;
                }
 
                // v2IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
                double v2IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v2mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                // v3IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
                double v3IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v3mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                // Direct fit to the collocation point at F2. We will take a weighted average of the direct and conditioned fit.
                double v2IM      = IMRPhenomX_Intermediate_Phase_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                /* Evaluate collocation points */
                pPhase->CollocationValuesPhaseInt[0] = phaseIN;
 
                // Take a weighted average for these points? Can help condition the fit.
                pPhase->CollocationValuesPhaseInt[1] = 0.75*(v2IMmRDv4 + RDv4) + 0.25*v2IM;
 
                // Use just v2 - v4RD to reconstruct the fit?
                //pPhase->CollocationValuesPhaseInt[1] = v2IMmRDv4 + RDv4);
 
                pPhase->CollocationValuesPhaseInt[2] = v3IMmRDv4 + RDv4;
 
                pPhase->CollocationValuesPhaseInt[3] = phaseRD;
 
                /* A_{0,i} */
                ff  = pPhase->CollocationPointsPhaseInt[0];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,0,0,1.0);
                gsl_matrix_set(A,0,1,ff1);
                gsl_matrix_set(A,0,2,ff2);
                gsl_matrix_set(A,0,3,ff3);
                gsl_vector_set(b,0,pPhase->CollocationValuesPhaseInt[0] - ff0);
 
                /* A_{1,i} */
                ff  = pPhase->CollocationPointsPhaseInt[1];
                ff1 = 1.0 / (ff / pWF->fRING);
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,1,0,1);
                gsl_matrix_set(A,1,1,ff1);
                gsl_matrix_set(A,1,2,ff2);
                gsl_matrix_set(A,1,3,ff3);
                gsl_vector_set(b,1,pPhase->CollocationValuesPhaseInt[1] - ff0);
 
                /* A_{2,i} */
                ff  = pPhase->CollocationPointsPhaseInt[2];
                ff1 = 1.0 / (ff / pWF->fRING);
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,2,0,1);
                gsl_matrix_set(A,2,1,ff1);
                gsl_matrix_set(A,2,2,ff2);
                gsl_matrix_set(A,2,3,ff3);
                gsl_vector_set(b,2,pPhase->CollocationValuesPhaseInt[2] - ff0);
 
                /* A_{3,i} */
                ff  = pPhase->CollocationPointsPhaseInt[3];
                ff1 = 1.0 / (ff / pWF->fRING);
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,3,0,1);
                gsl_matrix_set(A,3,1,ff1);
                gsl_matrix_set(A,3,2,ff2);
                gsl_matrix_set(A,3,3,ff3);
                gsl_vector_set(b,3,pPhase->CollocationValuesPhaseInt[3] - ff0);
 
                /* We now solve the system A x = b via an LU decomposition */
                gsl_linalg_LU_decomp(A,p,&s);
                gsl_linalg_LU_solve(A,p,b,x);
 
                /* Set intermediate phenomenological coefficients from solution to A x = b */
                pPhase->b0 = gsl_vector_get(x,0);                                                        // x[0] // Constant
                pPhase->b1 = gsl_vector_get(x,1) * pWF->fRING;                                           // x[1] // f^{-1}
                pPhase->b2 = gsl_vector_get(x,2) * pWF->fRING * pWF->fRING;                              // x[2] // f^{-2}
                //pPhase->b3 = 0.0;
                pPhase->b4 = gsl_vector_get(x,3) * pWF->fRING * pWF->fRING * pWF->fRING * pWF->fRING;    // x[3] // f^{-4}
 
                /* Tidy up in preparation for next GSL solve ... */
                gsl_vector_free(b);
                gsl_vector_free(x);
                gsl_matrix_free(A);
                gsl_permutation_free(p);
        }
        // Canonical intermediate model using 5 collocation points
        else if(pWF->IMRPhenomXIntermediatePhaseVersion == 105)
        {
                // Using 5 collocation points in intermediate region
                for(i = 0; i < 5; i++)
                {
                        fi = gpoints5[i] * deltax + xmin;
 
                        pPhase->CollocationPointsPhaseInt[i] = fi;
                }
 
                /* Evaluate collocation points */
 
                /* The first and last collocation points for the intermediate region are set from the inspiral fit and ringdown respectively */
                pPhase->CollocationValuesPhaseInt[0] = phaseIN;
                pPhase->CollocationValuesPhaseInt[4] = phaseRD;
 
                // v2IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
                double v2IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v2mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                // v3IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
                double v3IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v3mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                // Direct fit to the collocation point at F2. We will take a weighted average of the direct and conditioned fit.
                double v2IM      = IMRPhenomX_Intermediate_Phase_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
                // Take a weighted average for these points. Helps condition the fit.
                pPhase->CollocationValuesPhaseInt[1] = 0.75*(v2IMmRDv4 + RDv4) + 0.25*v2IM;
 
                pPhase->CollocationValuesPhaseInt[2] = v3IMmRDv4 + RDv4;
                pPhase->CollocationValuesPhaseInt[3] = IMRPhenomX_Intermediate_Phase_22_d43(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);
 
 
                // Collocation points: v4 = d43 + v3
                pPhase->CollocationValuesPhaseInt[3] = pPhase->CollocationValuesPhaseInt[3] + pPhase->CollocationValuesPhaseInt[2];
 
                /* A_{0,i} */
                ff  = pPhase->CollocationPointsPhaseInt[0];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff4 = ff2 * ff2;
                ff0 = (4.0 * pPhase->cL) / ((2.0*pWF->fDAMP)*(2.0*pWF->fDAMP) + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,0,0,1.0);
                gsl_matrix_set(A,0,1,ff1);
                gsl_matrix_set(A,0,2,ff2);
                gsl_matrix_set(A,0,3,ff3);
                gsl_matrix_set(A,0,4,ff4);
                gsl_vector_set(b,0,pPhase->CollocationValuesPhaseInt[0] - ff0);
 
                if(debug)
                {
                printf("For row 0: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f , ff0 = %.6f, ff = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[0] - ff0,ff0,ff);
                }
 
                /* A_{1,i} */
                ff  = pPhase->CollocationPointsPhaseInt[1];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff4 = ff2 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,1,0,1.0);
                gsl_matrix_set(A,1,1,ff1);
                gsl_matrix_set(A,1,2,ff2);
                gsl_matrix_set(A,1,3,ff3);
                gsl_matrix_set(A,1,4,ff4);
                gsl_vector_set(b,1,pPhase->CollocationValuesPhaseInt[1] - ff0);
 
                if(debug)
                {
                printf("For row 1: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[1] - ff0);
                }
 
                /* A_{2,i} */
                ff  = pPhase->CollocationPointsPhaseInt[2];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff4 = ff2 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,2,0,1.0);
                gsl_matrix_set(A,2,1,ff1);
                gsl_matrix_set(A,2,2,ff2);
                gsl_matrix_set(A,2,3,ff3);
                gsl_matrix_set(A,2,4,ff4);
                gsl_vector_set(b,2,pPhase->CollocationValuesPhaseInt[2] - ff0);
 
                if(debug)
                {
                printf("For row 2: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[2] - ff0);
                }
 
                /* A_{3,i} */
                ff  = pPhase->CollocationPointsPhaseInt[3];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff4 = ff2 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,3,0,1.0);
                gsl_matrix_set(A,3,1,ff1);
                gsl_matrix_set(A,3,2,ff2);
                gsl_matrix_set(A,3,3,ff3);
                gsl_matrix_set(A,3,4,ff4);
                gsl_vector_set(b,3,pPhase->CollocationValuesPhaseInt[3] - ff0);
 
                if(debug)
                {
                printf("For row 3: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[3] - ff0);
                }
 
                /* A_{4,i} */
                ff  = pPhase->CollocationPointsPhaseInt[4];
                ff1 = pWF->fRING / ff;
                ff2 = ff1 * ff1;
                ff3 = ff1 * ff2;
                ff4 = ff2 * ff2;
                ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
                gsl_matrix_set(A,4,0,1.0);
                gsl_matrix_set(A,4,1,ff1);
                gsl_matrix_set(A,4,2,ff2);
                gsl_matrix_set(A,4,3,ff3);
                gsl_matrix_set(A,4,4,ff4);
                gsl_vector_set(b,4,pPhase->CollocationValuesPhaseInt[4] - ff0);
 
                if(debug)
                {
                printf("For row 4: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[4] - ff0);
                }
 
                /* We now solve the system A x = b via an LU decomposition */
                gsl_linalg_LU_decomp(A,p,&s);
                gsl_linalg_LU_solve(A,p,b,x);
 
                /* Set intermediate phenomenological coefficients from solution to A x = b */
                pPhase->b0 = gsl_vector_get(x,0);                                                      // x[0] // Const.
                pPhase->b1 = gsl_vector_get(x,1) * pWF->fRING;                                         // x[1] // f^{-1}
                pPhase->b2 = gsl_vector_get(x,2) * pWF->fRING * pWF->fRING;                            // x[2] // f^{-2}
                pPhase->b3 = gsl_vector_get(x,3) * pWF->fRING * pWF->fRING * pWF->fRING;               // x[3] // f^{-3}
                pPhase->b4 = gsl_vector_get(x,4) * pWF->fRING * pWF->fRING * pWF->fRING * pWF->fRING;  // x[4] // f^{-4}
 
                if(debug)
                {
                printf("\n");
                printf("Intermediate Collocation Points and Values:\n");
                printf("F1 : %.7f\n",pPhase->CollocationPointsPhaseInt[0]);
                printf("F2 : %.7f\n",pPhase->CollocationPointsPhaseInt[1]);
                printf("F3 : %.7f\n",pPhase->CollocationPointsPhaseInt[2]);
                printf("F4 : %.7f\n",pPhase->CollocationPointsPhaseInt[3]);
                printf("F5 : %.7f\n",pPhase->CollocationPointsPhaseInt[4]);
                printf("\n");
                printf("V's agree with Mathematica...\n");
                printf("V1 : %.7f\n",pPhase->CollocationValuesPhaseInt[0]);
                printf("V2 : %.7f\n",pPhase->CollocationValuesPhaseInt[1]);
                printf("V3 : %.7f\n",pPhase->CollocationValuesPhaseInt[2]);
                printf("V4 : %.7f\n",pPhase->CollocationValuesPhaseInt[3]);
                printf("V5 : %.7f\n",pPhase->CollocationValuesPhaseInt[4]);
                printf("\n");
                printf("g0 : %.7f\n",gsl_vector_get(x,0));
                printf("g1 : %.7f\n",gsl_vector_get(x,1));
                printf("g2 : %.7f\n",gsl_vector_get(x,2));
                printf("g3 : %.7f\n",gsl_vector_get(x,3));
                printf("g4 : %.7f\n",gsl_vector_get(x,4));
                printf("\n");
                printf("b0 : %.7f\n",pPhase->b0);
                printf("b1 : %.7f\n",pPhase->b1);
                printf("b2 : %.7f\n",pPhase->b2);
                printf("b3 : %.7f\n",pPhase->b3);
                printf("b4 : %.7f\n",pPhase->b4);
                printf("\n");
                }
 
                /* Tidy up */
                gsl_vector_free(b);
                gsl_vector_free(x);
                gsl_matrix_free(A);
                gsl_permutation_free(p);
        }
        else
        {
                XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
        }
 
        /* Ringdown coefficients */
        REAL8 nonGR_dc1   = XLALSimInspiralWaveformParamsLookupNonGRDC1(LALparams);
        REAL8 nonGR_dc2   = XLALSimInspiralWaveformParamsLookupNonGRDC2(LALparams);
        REAL8 nonGR_dc4   = XLALSimInspiralWaveformParamsLookupNonGRDC4(LALparams);
        REAL8 nonGR_dcl   = XLALSimInspiralWaveformParamsLookupNonGRDCL(LALparams);
 
        /* Intermediate coefficients */
        REAL8 nonGR_db1   = XLALSimInspiralWaveformParamsLookupNonGRDB1(LALparams);
        REAL8 nonGR_db2   = XLALSimInspiralWaveformParamsLookupNonGRDB2(LALparams);
        REAL8 nonGR_db3   = XLALSimInspiralWaveformParamsLookupNonGRDB3(LALparams);
        REAL8 nonGR_db4   = XLALSimInspiralWaveformParamsLookupNonGRDB4(LALparams);
 
        /* Inspiral coefficients */
        REAL8 dchi_minus2 = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus2(LALparams);
        REAL8 dchi_minus1 = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus1(LALparams);
        REAL8 dchi0       = XLALSimInspiralWaveformParamsLookupNonGRDChi0(LALparams);
        REAL8 dchi1       = XLALSimInspiralWaveformParamsLookupNonGRDChi1(LALparams);
        REAL8 dchi2       = XLALSimInspiralWaveformParamsLookupNonGRDChi2(LALparams);
        REAL8 dchi3       = XLALSimInspiralWaveformParamsLookupNonGRDChi3(LALparams);
        REAL8 dchi4       = XLALSimInspiralWaveformParamsLookupNonGRDChi4(LALparams);
        REAL8 dchi5       = XLALSimInspiralWaveformParamsLookupNonGRDChi5(LALparams);
        REAL8 dchi5L      = XLALSimInspiralWaveformParamsLookupNonGRDChi5L(LALparams);
        REAL8 dchi6       = XLALSimInspiralWaveformParamsLookupNonGRDChi6(LALparams);
        REAL8 dchi6L      = XLALSimInspiralWaveformParamsLookupNonGRDChi6L(LALparams);
        REAL8 dchi7       = XLALSimInspiralWaveformParamsLookupNonGRDChi7(LALparams);
 
        /* Can include these terms in the future as desired... */
        REAL8 dchi8       = 0.0;
        REAL8 dchi8L      = 0.0;
        REAL8 dchi9       = 0.0;
        REAL8 dchi9L      = 0.0;
 
        /* ~~~~ RINGDOWN ~~~~ */
        pPhase->cLGR  = pPhase->cL; // Store GR value for reference
        pPhase->c1   *= (1.0 + nonGR_dc1);
        pPhase->c2   *= (1.0 + nonGR_dc2);
        pPhase->c4   *= (1.0 + nonGR_dc4);
        pPhase->cL   *= (1.0 + nonGR_dcl);
 
        /* Set pre-cached variables */
        pPhase->c4ov3   = pPhase->c4 / 3.0;
        pPhase->cLovfda = pPhase->cL / pWF->fDAMP;
 
        /* Apply NR tuning for precessing cases (500) */
        pPhase->b1 = pPhase->b1  +  ( pWF->PNR_DEV_PARAMETER * pWF->ZETA2 );
        pPhase->b4 = pPhase->b4  +  ( pWF->PNR_DEV_PARAMETER * pWF->ZETA1 );
 
        /* ~~~~ INTERMEDIATE ~~~~ */
        if(pWF->IMRPhenomXIntermediatePhaseVersion == 104)
        {
                pPhase->b1 *= (1.0 + nonGR_db1);
                pPhase->b2 *= (1.0 + nonGR_db2);
                pPhase->b4 *= (1.0 + nonGR_db4);
        }
        else if(pWF->IMRPhenomXIntermediatePhaseVersion == 105)
        {
                pPhase->b1 *= (1.0 + nonGR_db1);
                pPhase->b2 *= (1.0 + nonGR_db2);
                pPhase->b3 *= (1.0 + nonGR_db3);
                pPhase->b4 *= (1.0 + nonGR_db4);
        }
        else
        {
                XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
        }
 
        /* ~~~~ INSPIRAL ~~~~ */
        /* Initialize -1PN coefficient*/
        pPhase->phi_minus2   = 0.0;
        pPhase->dphi_minus2  = 0.0;
 
        pPhase->phi_minus1   = 0.0;
        pPhase->dphi_minus1  = 0.0;
 
        /*
                If tgr_parameterization = 1, deform complete PN coefficient. This is an FTA-like parameterization.
                If tgr_parameterization = 0, only deform non-spinning coefficient. This is the original TIGER-like implementation.
        */
        int tgr_parameterization = 0;
        tgr_parameterization     = XLALSimInspiralWaveformParamsLookupNonGRParameterization(LALparams);
 
        if(tgr_parameterization == 1)
        {
                        /* -1.0 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi_minus2 = dchi_minus2 / powers_of_lalpi.two_thirds;
 
                        /* -0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi_minus1 = dchi_minus1 / powers_of_lalpi.one_third;
 
                        /* 0.0 PN */
                        pPhase->phi0       = (phi0NS + phi0S)*(1.0 + dchi0);
 
                        /* 0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi1       = dchi1 * powers_of_lalpi.one_third;
 
                        /* 1.0 PN */
                        pPhase->phi2       = (phi2NS + phi2S)*(1.0 + dchi2);
 
                        // /* 1.5 PN */
                        // pPhase->phi3       = (phi3NS + phi3S)*(1.0 + dchi3);
                        /* 1.5 PN */
                        pPhase->phi3       = dchi3;

                        /* 2.0 PN */
                        pPhase->phi4       = (phi4NS + phi4S)*(1.0 + dchi4);
 
                        /* 2.5 PN */
                        pPhase->phi5       = (phi5NS + phi5S)*(1.0 + dchi5);
 
                        /* 2.5 PN, Log Terms */
                        pPhase->phi5L      = (phi5LNS + phi5LS)*(1.0 + dchi5L);
 
                        /* 3.0 PN */
                        pPhase->phi6       = (phi6NS + phi6S)*(1.0 + dchi6);
 
                        /* 3.0 PN, Log Term */
                        pPhase->phi6L      = (phi6LNS + phi6LS)*(1.0 + dchi6L);
 
                        /* 3.5PN */
                        pPhase->phi7       = (phi7NS + phi7S)*(1.0 + dchi7);
 
                        /* 4.0PN */
                        pPhase->phi8       = (phi8NS + phi8S)*(1.0 + dchi8);
 
                        /* 4.0 PN, Log Terms */
                        pPhase->phi8L      = (phi8LNS + phi8LS)*(1.0 + dchi8L);
 
                        /* 4.0 PN */
                        pPhase->phi9       = (phi9NS + phi9S)*(1.0 + dchi9);
 
                        /* 4.0 PN, Log Terms */
                        pPhase->phi9L      = (phi9LNS + phi9LS)*(1.0 + dchi9L);
        }
        else if(tgr_parameterization == 0)
        {
                        /* -1.0 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi_minus2 = dchi_minus2 / powers_of_lalpi.two_thirds;
 
                        /* -0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi_minus1 = dchi_minus1 / powers_of_lalpi.one_third;
 
                        /* 0.0 PN */
                        pPhase->phi0       = phi0NS*(1.0 + dchi0) + phi0S;
 
                        /* 0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
                        pPhase->phi1       = dchi1 * powers_of_lalpi.one_third;
 
                        /* 1.0 PN */
                        pPhase->phi2       = phi2NS*(1.0 + dchi2) + phi2S;
 
                        // /* 1.5 PN */
                        // pPhase->phi3       = phi3NS*(1.0 + dchi3)+ phi3S;
                        /* 1.5 PN */
                        pPhase->phi3       = dchi3;

                        /* 2.0 PN */
                        pPhase->phi4       = phi4NS*(1.0 + dchi4) + phi4S;
 
                        /* 2.5 PN */
                        pPhase->phi5       = phi5NS*(1.0 + dchi5) + phi5S;
 
                        /* 2.5 PN, Log Terms */
                        pPhase->phi5L      = phi5LNS*(1.0 + dchi5L) + phi5LS;
 
                        /* 3.0 PN */
                        pPhase->phi6       = phi6NS*(1.0 + dchi6) + phi6S;
 
                        /* 3.0 PN, Log Term */
                        pPhase->phi6L      = phi6LNS*(1.0 + dchi6L) + phi6LS;
 
                        /* 3.5PN */
                        pPhase->phi7       = phi7NS*(1.0 + dchi7) + phi7S;
 
                        /* 4.0PN */
                        pPhase->phi8       = phi8NS*(1.0 + dchi8) + phi8S;
 
                        /* 4.0 PN, Log Terms */
                        pPhase->phi8L      = phi8LNS*(1.0 + dchi8L) + phi8LS;
 
                        /* 4.0 PN */
                        pPhase->phi9       = phi9NS*(1.0 + dchi9) + phi9S;
 
                        /* 4.0 PN, Log Terms */
                        pPhase->phi9L      = phi9LNS*(1.0 + dchi9L) + phi9LS;
        }
        else
        {
                        XLALPrintError("Error in IMRPhenomXGetPhaseCoefficients: TGR Parameterizataion is not valid.\n");
        }
 
        /* Recalculate phase derivatives including TGR corrections */
        pPhase->dphi_minus2 = +(7.0 / 5.0) * pPhase->phi_minus2;
        pPhase->dphi_minus1 = +(6.0 / 5.0) * pPhase->phi_minus1;
        pPhase->dphi0       = +(5.0 / 5.0) * pPhase->phi0;
        pPhase->dphi1       = +(4.0 / 5.0) * pPhase->phi1;
        pPhase->dphi2       = +(3.0 / 5.0) * pPhase->phi2;
        pPhase->dphi3       = +(2.0 / 5.0) * pPhase->phi3;
        pPhase->dphi4       = +(1.0 / 5.0) * pPhase->phi4;
        pPhase->dphi5       = -(3.0 / 5.0) * pPhase->phi5L;
        pPhase->dphi6       = -(1.0 / 5.0) * pPhase->phi6 - (3.0 / 5.0) * pPhase->phi6L;
        pPhase->dphi6L      = -(1.0 / 5.0) * pPhase->phi6L;
        pPhase->dphi7       = -(2.0 / 5.0) * pPhase->phi7;
        pPhase->dphi8       = -(3.0 / 5.0) * pPhase->phi8 - (3.0 / 5.0) * pPhase->phi8L;
        pPhase->dphi8L      = -(3.0 / 5.0) * pPhase->phi8L;
        pPhase->dphi9       = -(4.0 / 5.0) * pPhase->phi9 - (3.0 / 5.0) * pPhase->phi9L;
        pPhase->dphi9L      = -(3.0 / 5.0) * pPhase->phi9L;
 
        /* Initialize connection coefficients */
        pPhase->C1Int = 0;
        pPhase->C2Int = 0;
        pPhase->C1MRD = 0;
        pPhase->C2MRD = 0;
 
        /* END OF ROUTINE; RETURN STRUCT */
        return XLAL_SUCCESS;
}
 




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


/***********************************************************************************************************************************************
 * 模块名称：调试辅助
 * 功能描述：提供调试信息输出和错误跟踪功能
 * 备注信息：发布版本中会被编译宏禁用
 **********************************************************************************************************************************************/
void IMRPhenomX_Phase_22_ConnectionCoefficients(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase)
{
        int debug = pWF->debug;
 
        double fIns = pPhase->fPhaseMatchIN;
        double fInt = pPhase->fPhaseMatchIM;
 
        /*
                        Assume an ansatz of the form:
 
                        phi_Inspiral (f) = phi_Intermediate (f) + C1 + C2 * f
 
                        where transition frequency is fIns
 
                        phi_Inspiral (fIns) = phi_Intermediate (fIns) + C1 + C2 * fIns
                        phi_Inspiral'(fIns) = phi_Intermediate'(fIns) + C2
 
                        Solving for C1 and C2
 
                        C2 = phi_Inspiral'(fIns) - phi_Intermediate'(fIns)
                        C1 = phi_Inspiral (fIns) - phi_Intermediate (fIns) - C2 * fIns
 
        */
        IMRPhenomX_UsefulPowers powers_of_fIns;
        IMRPhenomX_Initialize_Powers(&powers_of_fIns,fIns);
 
        double DPhiIns = IMRPhenomX_Inspiral_Phase_22_Ansatz(fIns,&powers_of_fIns,pPhase);
        double DPhiInt = IMRPhenomX_Intermediate_Phase_22_Ansatz(fIns,&powers_of_fIns,pWF,pPhase);
 
        pPhase->C2Int  = DPhiIns - DPhiInt;
 
        double phiIN = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(fIns,&powers_of_fIns,pPhase);
        double phiIM = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(fIns,&powers_of_fIns,pWF,pPhase);
 
        if(debug)
        {
        printf("\n");
        printf("dphiIM = %.6f and dphiIN = %.6f\n",DPhiInt,DPhiIns);
        printf("phiIN(fIns)  : %.7f\n",phiIN);
        printf("phiIM(fIns)  : %.7f\n",phiIM);
        printf("fIns         : %.7f\n",fIns);
        printf("C2           : %.7f\n",pPhase->C2Int);
        printf("\n");
        }
 
        pPhase->C1Int = phiIN - phiIM - (pPhase->C2Int * fIns);
 
        /*
                        Assume an ansatz of the form:
 
                        phi_Intermediate (f)    = phi_Ringdown (f) + C1 + C2 * f
 
                        where transition frequency is fIM
 
                        phi_Intermediate (fIM) = phi_Ringdown (fRD) + C1 + C2 * fIM
                        phi_Intermediate'(fIM) = phi_Ringdown'(fRD) + C2
 
                        Solving for C1 and C2
 
                        C2 = phi_Inspiral'(fIM) - phi_Intermediate'(fIM)
                        C1 = phi_Inspiral (fIM) - phi_Intermediate (fIM) - C2 * fIM
 
        */
        IMRPhenomX_UsefulPowers powers_of_fInt;
        IMRPhenomX_Initialize_Powers(&powers_of_fInt,fInt);
 
        double phiIMC         = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(fInt,&powers_of_fInt,pWF,pPhase) + pPhase->C1Int + pPhase->C2Int*fInt;
        double phiRD          = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(fInt,&powers_of_fInt,pWF,pPhase);
        double DPhiIntC       = IMRPhenomX_Intermediate_Phase_22_Ansatz(fInt,&powers_of_fInt,pWF,pPhase) + pPhase->C2Int;
        double DPhiRD         = IMRPhenomX_Ringdown_Phase_22_Ansatz(fInt,&powers_of_fInt,pWF,pPhase);
 
        pPhase->C2MRD = DPhiIntC - DPhiRD;
        pPhase->C1MRD = phiIMC - phiRD - pPhase->C2MRD*fInt;
 
        if(debug)
        {
        printf("\n");
        printf("phiIMC(fInt) : %.7f\n",phiIMC);
        printf("phiRD(fInt)  : %.7f\n",phiRD);
        printf("fInt         : %.7f\n",fInt);
        printf("C2           : %.7f\n",pPhase->C2Int);
        printf("\n");
        }
 
        if(debug)
        {
        printf("dphiIM = %.6f and dphiRD = %.6f\n",DPhiIntC,DPhiRD);
        printf("\nContinuity Coefficients\n");
        printf("C1Int : %.6f\n",pPhase->C1Int);
        printf("C2Int : %.6f\n",pPhase->C2Int);
        printf("C1MRD : %.6f\n",pPhase->C1MRD);
        printf("C2MRD : %.6f\n",pPhase->C2MRD);
        }
 
  return;
}


/*
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 * This function computes the IMRPhenomX phase given a phase coefficients struct pPhase.
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 */
double IMRPhenomX_Phase_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF)
{
        // Inspiral region, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIN))
  {
          double PhiIns = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(ff, powers_of_f, pPhase);
          return PhiIns;
  }
 
        // Ringdown region, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIM))
  {
          double PhiMRD = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(ff, powers_of_f, pWF, pPhase)
                            + pPhase->C1MRD + (pPhase->C2MRD * ff);
          return PhiMRD;
  }
 
  //    Intermediate region, fPhaseInsMax < f < fPhaseIntMax
  double PhiInt = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(ff, powers_of_f, pWF, pPhase)
                            + pPhase->C1Int + (pPhase->C2Int * ff);
 
  return PhiInt;
}
 
/*
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 * This function computes the IMRPhenomX phase derivative given a phase coefficients struct pPhase.
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 */
double IMRPhenomX_dPhase_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF)
{
        // Inspiral region, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIN))
  {
          double dPhiIns = IMRPhenomX_Inspiral_Phase_22_Ansatz(ff, powers_of_f, pPhase);
          return dPhiIns;
  }
 
        // Ringdown region, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIM))
  {
          double dPhiMRD = IMRPhenomX_Ringdown_Phase_22_Ansatz(ff, powers_of_f, pWF, pPhase)
                            + (pPhase->C2MRD);
          return dPhiMRD;
  }
 
  //    Intermediate region, fPhaseInsMax < f < fPhaseIntMax
  double dPhiInt = IMRPhenomX_Intermediate_Phase_22_Ansatz(ff, powers_of_f, pWF, pPhase)
                            + (pPhase->C2Int);
 
  return dPhiInt;
}
 


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