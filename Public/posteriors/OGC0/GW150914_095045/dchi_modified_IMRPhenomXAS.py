from WF4Py import waveforms as WF;
from WF4Py import WFutils as utils;

import numpy as np
from astropy.cosmology import Planck18

import pycbc.conversions as cvs;
import pycbc.waveform as wf;
import pycbc.types as tps;

import time
import sys
import os

import matplotlib.pyplot as pp
from matplotlib import gridspec
from matplotlib import rc
rc('text', usetex=False)

ourcodenameString = 'WF4Py'


def get_dchi_modified_IMRPhenomXAS(**kwargs) :
    # // IMPORTANT CHECK: verify that inclination_azimuth is the longitude of ascending nodes

    # print(f"kwargs: {kwargs}")
    mass1 = kwargs["mass1"]
    mass2 = kwargs["mass2"]
    chi1z = kwargs["spin1z"]
    chi2z = kwargs["spin2z"]
    dL = kwargs["distance"]

    phi0 = kwargs["coa_phase"]
    # kwargs["long_asc_nodes"]
    # kwargs["mean_per_ano"]

    events = dict(
        Mc =  np.array([cvs.mchirp_from_mass1_mass2(mass1, mass2)]),
        eta = np.array([cvs.eta_from_mass1_mass2(mass1, mass2)   ]),
        chi1z = np.array([ chi1z ]),
        chi2z = np.array([ chi2z ]),
        dL    = np.array([ dL/1000. ]),
        f_lower = np.array([ kwargs["f_lower"] ]) if "f_lower" in kwargs else 15
    )

    if "dchi0" in kwargs and kwargs.get("if_dchi0", False):
        events["dchi0"] = kwargs["dchi0"]
    if "dchi1" in kwargs and kwargs.get("if_dchi1", False):
        events["dchi1"] = kwargs["dchi1"]
    if "dchi2" in kwargs and kwargs.get("if_dchi2", False):
        events["dchi2"] = kwargs["dchi2"]
    if "dchi3" in kwargs and kwargs.get("if_dchi3", False):
        events["dchi3"] = kwargs["dchi3"]
    if "dchi4" in kwargs and kwargs.get("if_dchi4", False):
        events["dchi4"] = kwargs["dchi4"]
    if "dchi5l" in kwargs and kwargs.get("if_dchi5l", False):
        events["dchi5l"] = kwargs["dchi5l"]
    if "dchi5" in kwargs and kwargs.get("if_dchi5", False):
        events["dchi5"] = kwargs["dchi5"]
    if "dchi6" in kwargs and kwargs.get("if_dchi6", False):
        events["dchi6"] = kwargs["dchi6"]
    if "dchi6l" in kwargs and kwargs.get("if_dchi6l", False):
        events["dchi6l"] = kwargs["dchi6l"]
    if "dchi7" in kwargs and kwargs.get("if_dchi7", False):
        events["dchi7"] = kwargs["dchi7"]
    if "dchi8" in kwargs and kwargs.get("if_dchi8", False):
        events["dchi8"] = kwargs["dchi8"]
    if "dchi8l" in kwargs and kwargs.get("if_dchi8l", False):
        events["dchi8l"] = kwargs["dchi8l"]
    


    tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=104, IntPhaseVersion=105, IntAmpVersion=104)

    if "f_final" in kwargs :
        if kwargs["f_final"] < kwargs["f_lower"] :
            events["f_final"] = tmpWF.fcut(**events)
        else :
            events["f_final"] = np.array([kwargs["f_final"]])
    else :
        events["f_final"] = tmpWF.fcut(**events)

    # fcut = tmpWF.fcut(**events)
    # fminarr = np.full(events["f_final"], kwargs["f_lower"])
    delta_f = kwargs["delta_f"]
    fgrids = np.arange(0, events["f_final"], delta_f)
    fgrids = fgrids.reshape(-1, 1)

    xas_ampl = tmpWF.Ampl(fgrids, **events)
    xas_phase = tmpWF.Phi(fgrids, **events) + 2*phi0
    # xas_phase0 = tmpWF.Phi0(fgrids, **events) + 2*phi0

    def htilde22(phase, ampl) :
        return ampl * np.exp(-1j*phase)

    iota = kwargs["inclination"]

    def hphc(phase, ampl, iota) :
        cfac = np.cos(iota)
        Ylmfactor = 2.0 * np.sqrt(5.0 / (64.0 * np.pi)) * np.exp(-1j * 2 * (np.pi/2))
        pfac = 0.5 * (1. + cfac*cfac)
        hp = htilde22(phase, ampl)
        hc = -1j * cfac * hp * Ylmfactor
        hp *= pfac*Ylmfactor

        return hp, hc
    
    Y_NORM = 2. * np.sqrt(5./(64.*np.pi))
    
    hp, hc = hphc(xas_phase, xas_ampl, iota)

    return tps.FrequencySeries(-np.conj(hp.reshape(1, -1)[0])/Y_NORM, delta_f=delta_f), tps.FrequencySeries(np.conj(hc.reshape(1, -1)[0])/Y_NORM, delta_f=delta_f)


# def get_dchi_modified_IMRPhenomXAS(**kwargs) :
#     # // IMPORTANT CHECK: verify that inclination_azimuth is the longitude of ascending nodes

#     mass1 = kwargs["mass1"]
#     mass2 = kwargs["mass2"]
#     chi1z = kwargs["spin1z"]
#     chi2z = kwargs["spin2z"]
#     dL = kwargs["distance"]

#     phi0 = kwargs["coa_phase"]
#     # kwargs["long_asc_nodes"]
#     # kwargs["mean_per_ano"]

#     events = dict(
#         Mc =  np.array([cvs.mchirp_from_mass1_mass2(mass1, mass2)]),
#         eta = np.array([cvs.eta_from_mass1_mass2(mass1, mass2)   ]),
#         chi1z = np.array([ chi1z ]),
#         chi2z = np.array([ chi2z ]),
#         dL    = np.array([ dL/1000. ]),
#         f_lower = np.array([ kwargs["f_lower"] ]) if "f_lower" in kwargs else 15
#     )

#     print(events)

#     if "dchi0" in kwargs:
#         events["dchi0"] = kwargs["dchi0"]
#     elif "dchi1" in kwargs:
#         events["dchi1"] = kwargs["dchi1"]
#     elif "dchi2" in kwargs:
#         events["dchi2"] = kwargs["dchi2"]
#     elif "dchi3" in kwargs:
#         events["dchi3"] = kwargs["dchi3"]
#     elif "dchi4" in kwargs:
#         events["dchi4"] = kwargs["dchi4"]
#     elif "dchi5L" in kwargs:
#         events["dchi5L"] = kwargs["dchi5L"]
#     elif "dchi5" in kwargs:
#         events["dchi5"] = kwargs["dchi5"]
#     elif "dchi6" in kwargs:
#         events["dchi6"] = kwargs["dchi6"]
#     elif "dchi6L" in kwargs:
#         events["dchi6L"] = kwargs["dchi6L"]
#     elif "dchi7" in kwargs:
#         events["dchi7"] = kwargs["dchi7"]
#     elif "dchi8" in kwargs:
#         events["dchi8"] = kwargs["dchi8"]
#     elif "dchi8L" in kwargs:
#         events["dchi8L"] = kwargs["dchi8L"]
#     else:
#         pass

#     tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=104, IntPhaseVersion=105, IntAmpVersion=104)

#     events["f_final"] = np.array([kwargs["f_final"]]) if "f_final" in kwargs else tmpWF.fcut(**events)
#     # fcut = tmpWF.fcut(**events)
#     fminarr = np.full(events["f_final"], kwargs["f_lower"])
#     delta_f = kwargs["delta_f"]
#     fgrids = np.arange(0, events["f_final"], delta_f)
#     fgrids = fgrids.reshape(-1, 1)

#     print(fgrids)

#     xas_ampl = tmpWF.Ampl(fgrids, **events) + 2*phi0
#     xas_phase22 = tmpWF.Phi(fgrids, **events)

#     iota = kwargs["inclination"]
        
#     Y_NORM = 2. * np.sqrt(5./(64.*np.pi))
    
#     hp, hc = utils.Add_Higher_Modes(Ampl={"22": xas_ampl}, Phi={"22": xas_phase22}, iota=iota)
#     # hp, hc = utils.Add_Higher_Modes(Ampl={"22": xas_ampl, "2-2": xas_ampl}, Phi={"22": xas_phase22, "2-2": -xas_phase22}, iota=iota)
#     return tps.FrequencySeries(-np.conj(hp.reshape(1, -1)[0])/Y_NORM, delta_f=delta_f), tps.FrequencySeries(np.conj(hc.reshape(1, -1)[0])/Y_NORM, delta_f=delta_f)
