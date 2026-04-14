from WF4Py import waveforms as WF;
from WF4Py import WFutils as utils;

import numpy as np
from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
from astropy import units as u

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

import numpy as np
# 保留原代码的依赖导入（根据你的环境补充）
# import cvs, WF, tps


import lal as _lal

from pycbc.filter import match
from pycbc.psd import aLIGOZeroDetHighPower



# def get_dchi_modified_IMRPhenomXAS(**params) :
#     params.pop('approximant')
#     return wf.get_fd_waveform(approximant='IMRPhenomXAS', **params)


# def get_dchi_modified_IMRPhenomXAS(**kwargs) :
#     # // IMPORTANT CHECK: verify that inclination_azimuth is the longitude of ascending nodes
#     # print(f"kwargs: {kwargs}")
    
    
    
#     hp_gr, hc_gr = wf.get_fd_waveform(approximant='IMRPhenomXAS',
#                                                     mass1=kwargs['mass1'],
#                                                     mass2=kwargs['mass2'],
#                                                     spin1z=kwargs['spin1z'],
#                                                     spin2z=kwargs['spin2z'],
#                                                     delta_f=kwargs["delta_f"],
#                                                     distance=kwargs['distance'],
#                                                     inclination=kwargs['inclination'],
#                                                     f_lower=kwargs["f_lower"],
#                                                     )
    
#     # 1. 新增极化角参数校验，避免KeyError（兼容旧代码调用）
#     psi = 0
#     # psi = kwargs.get("polarization", 0.0)  # 无polarization参数时默认0.0

#     mass1 = kwargs["mass1"]
#     mass2 = kwargs["mass2"]
#     chi1z = kwargs["spin1z"]
#     chi2z = kwargs["spin2z"]
#     dL = kwargs["distance"]
    
#     iota = kwargs["inclination"]
#     phi0 = kwargs["coa_phase"]
    
#     f_ref = kwargs["f_ref"]
    
#     long_asc_nodes = kwargs.get("long_asc_nodes", 0.0)  # 升交点经度Ω，默认0
#     mean_per_ano = kwargs.get("mean_per_ano", 0.0)      # 近心点幅角ω，默认0
    
#     # print(f"long_asc_nodes: {long_asc_nodes}")
#     # print(f"mean_per_ano: {mean_per_ano}")

#     if mass1 < mass2 :
#         temp = mass1
#         mass1 = mass2
#         mass2 = temp
        
#         temp = chi1z
#         chi1z = chi2z
#         chi2z = temp
    
#     # zs = z_at_value(Planck18.luminosity_distance, dL * u.Mpc).value
#     # print(f"zs: {zs}")
#     events = dict(
#         Mc =  np.array([cvs.mchirp_from_mass1_mass2(mass1, mass2)])*(1),
#         eta = np.array([cvs.eta_from_mass1_mass2(mass1, mass2)   ]),
#         chi1z = np.array([ chi1z ]),
#         chi2z = np.array([ chi2z ]),
#         dL    = np.array([ dL/1000. ]),
#         f_lower = np.array([ kwargs["f_lower"] ]) if "f_lower" in kwargs else np.array([15]),
#         f_ref = np.array([ f_ref ])
#     )
        
#     # 双精度浮点数（Python float 类型）的机器epsilon
#     eps = np.finfo(np.float64).eps
#     if events["eta"] > 0.25 and np.isclose(events["eta"], 0.25, atol=eps) :
#         events["eta"] = 0.25
#     if events["eta"] < 0 and np.isclose(events["eta"], 0, atol=eps) :
#         events["eta"] = 0

#     if events["chi1z"] > 1 and np.isclose(events["chi1z"], 1, atol=eps) :
#         events["chi1z"] = 1
#     if events["chi1z"] < -1 and np.isclose(events["chi1z"], -1, atol=eps) :
#         events["chi1z"] = -1

#     if events["chi2z"] > 1 and np.isclose(events["chi2z"], 1, atol=eps) :
#         events["chi2z"] = 1
#     if events["chi2z"] < -1 and np.isclose(events["chi2z"], -1, atol=eps) :
#         events["chi2z"] = -1

#     if events["f_lower"] < 0 and np.isclose(events["f_lower"], 0, atol=eps) :
#         events["f_lower"] = 0

#     dchi_keys = ["dchi0", "dchi1", "dchi2", "dchi3", "dchi4", "dchi5l", "dchi5", "dchi6", "dchi6l", "dchi7", "dchi8", "dchi8l"]
#     for key in dchi_keys:
#         if key in kwargs and kwargs.get(f"if_{key}", False):
#             events[key] = kwargs[key]
#             events[f"if_{key}"] = kwargs[f"if_{key}"]

#     # tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=104, IntPhaseVersion=105, IntAmpVersion=104)
#     tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=115, IntPhaseVersion=105, IntAmpVersion=104)

#     if "f_final" in kwargs :
#         if kwargs["f_final"] < kwargs["f_lower"] :
#             # events["f_final"] = tmpWF.fcut(**events)
#             events["f_final"] = np.array([2**(np.floor(np.log2(kwargs["f_lower"]))+1)])
#         else :
#             events["f_final"] = np.array([2**(np.floor(np.log2(kwargs["f_final"]))+1)])
#     else :
        
#         events["f_final"] = 2**(np.floor(np.log2(tmpWF.fcut(**events)))+1)

#     events["f_final"] = np.array([2048])
    
#     delta_f = kwargs["delta_f"]
#     # fgrids = np.arange(0, events["f_final"]+delta_f/10., delta_f).reshape(-1, 1)  # 保留频率网格

#     fgrids = hp_gr.sample_frequencies.data.reshape(-1, 1)

#     # 获取现象学振幅和相位
#     xas_ampl = tmpWF.Ampl(fgrids, **events)
#     # 2. 修正相位计算：还原旧代码规范，移除多余的-2psi和错误外层负号
#     xas_phase = tmpWF.Phi(fgrids, **events) + 2*phi0 + 2*(long_asc_nodes + mean_per_ano)
    
#     # Y22 = 2.0 * np.sqrt(5.0 / (64.0 * np.pi))，wf4py里面amp乘过了
#     htilde22 = xas_ampl * np.exp(-1j * xas_phase)
    
#     # 4. 修正极化因子计算：匹配旧代码规范，修复hc虚数符号
#     cfac = np.cos(iota)
#     pfac_hp = 0.5 * (1.0 + cfac ** 2)  # hp极化因子（不变）
#     pfac_hc = -cfac                    # hc极化因子：移除多余的-1j（后续共轭会处理虚部）

#     # hp = htilde22
#     # hc = htilde22
    
#     hp = np.conj(htilde22)
#     hc = 1j*hp

#     # # 5. 恢复Y22归一化：乘入hp/hc，修正振幅量级
#     hp *= pfac_hp
#     hc *= pfac_hc
    
#     cos2psi = np.cos(2 * psi)
#     sin2psi = np.sin(2 * psi)
#     hp, hc = hp * cos2psi - hc * sin2psi, hp * sin2psi + hc * cos2psi

#     # # 6. 还原旧代码核心规范：复共轭+符号修正（最关键，匹配观测数据极化相位）
#     # hp = np.conj(hp)
#     # hc = hp * np.exp(-1j*np.pi/2)

#     # 7. 统一维度处理：与旧代码一致，显式转为一维，避免解析异常
#     hp_flat = hp.reshape(1, -1)[0]
#     hc_flat = hc.reshape(1, -1)[0]

#     # 转换为FrequencySeries
#     # hp_series = tps.FrequencySeries(hp_flat, delta_f=delta_f)
#     # hc_series = tps.FrequencySeries(hc_flat, delta_f=delta_f)
    

    
#     # psd = aLIGOZeroDetHighPower(len(hp_gr.sample_frequencies.data), kwargs["delta_f"], 0)

#     # hp_flat = np.nan_to_num(hp_flat, 0)
#     # hc_flat = np.nan_to_num(hc_flat, 0)
    
#     hp_series = tps.FrequencySeries(hp_flat, delta_f=delta_f, epoch=_lal.LIGOTimeGPS(-16.0))
#     hc_series = tps.FrequencySeries(hc_flat, delta_f=delta_f, epoch=_lal.LIGOTimeGPS(-16.0))

#     # m, _ = match(hp_gr, hp_series, psd=psd, low_frequency_cutoff=kwargs["f_lower"])
    
#     # if 1 - m > 1e-4:
#     #     hp_series, hc_series = hp_gr, hc_gr
#     # else:
#     #     pass
    


#     return hp_series, hc_series



def get_dchi_modified_IMRPhenomXAS(**kwargs) :
    # // IMPORTANT CHECK: verify that inclination_azimuth is the longitude of ascending nodes
    # print(f"kwargs: {kwargs}")
    
    # 1. 新增极化角参数校验，避免KeyError（兼容旧代码调用）
    psi = 0
    # psi = kwargs.get("polarization", 0.0)  # 无polarization参数时默认0.0

    mass1 = kwargs["mass1"]
    mass2 = kwargs["mass2"]
    chi1z = kwargs["spin1z"]
    chi2z = kwargs["spin2z"]
    dL = kwargs["distance"]
    
    iota = kwargs["inclination"]
    phi0 = kwargs["coa_phase"]
    
    f_ref = kwargs["f_ref"]
    
    long_asc_nodes = kwargs.get("long_asc_nodes", 0.0)  # 升交点经度Ω，默认0
    mean_per_ano = kwargs.get("mean_per_ano", 0.0)      # 近心点幅角ω，默认0
    
    # print(f"long_asc_nodes: {long_asc_nodes}")
    # print(f"mean_per_ano: {mean_per_ano}")

    if mass1 < mass2 :
        temp = mass1
        mass1 = mass2
        mass2 = temp
        
        temp = chi1z
        chi1z = chi2z
        chi2z = temp
    
    # zs = z_at_value(Planck18.luminosity_distance, dL * u.Mpc).value
    # print(f"zs: {zs}")
    events = dict(
        Mc =  np.array([cvs.mchirp_from_mass1_mass2(mass1, mass2)])*(1),
        eta = np.array([cvs.eta_from_mass1_mass2(mass1, mass2)   ]),
        chi1z = np.array([ chi1z ]),
        chi2z = np.array([ chi2z ]),
        dL    = np.array([ dL/1000. ]),
        f_lower = np.array([ kwargs["f_lower"] ]) if "f_lower" in kwargs else np.array([15]),
        f_ref = np.array([ f_ref ])
    )
        
    # 双精度浮点数（Python float 类型）的机器epsilon
    eps = np.finfo(np.float64).eps
    if events["eta"] > 0.25 and np.isclose(events["eta"], 0.25, atol=eps) :
        events["eta"] = 0.25
    if events["eta"] < 0 and np.isclose(events["eta"], 0, atol=eps) :
        events["eta"] = 0

    if events["chi1z"] > 1 and np.isclose(events["chi1z"], 1, atol=eps) :
        events["chi1z"] = 1
    if events["chi1z"] < -1 and np.isclose(events["chi1z"], -1, atol=eps) :
        events["chi1z"] = -1

    if events["chi2z"] > 1 and np.isclose(events["chi2z"], 1, atol=eps) :
        events["chi2z"] = 1
    if events["chi2z"] < -1 and np.isclose(events["chi2z"], -1, atol=eps) :
        events["chi2z"] = -1

    if events["f_lower"] < 0 and np.isclose(events["f_lower"], 0, atol=eps) :
        events["f_lower"] = 0

    dchi_keys = ["dchi0", "dchi1", "dchi2", "dchi3", "dchi4", "dchi5l", "dchi5", "dchi6", "dchi6l", "dchi7", "dchi8", "dchi8l"]
    for key in dchi_keys:
        if key in kwargs and kwargs.get(f"if_{key}", False):
            events[key] = kwargs[key]
            events[f"if_{key}"] = kwargs[f"if_{key}"]

    # choose_origin=False, InsPhaseVersion=104, IntPhaseVersion=105, IntAmpVersion=105, 
    InsPhaseVersion = kwargs.get("InsPhaseVersion", 115)
    IntPhaseVersion = kwargs.get("IntPhaseVersion", 105)
    IntAmpVersion = kwargs.get("IntAmpVersion", 105)
    # tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=115, IntPhaseVersion=105, IntAmpVersion=105)
    tmpWF = WF.IMRPhenomXAS(InsPhaseVersion=InsPhaseVersion, IntPhaseVersion=IntPhaseVersion, IntAmpVersion=IntAmpVersion)

    # if "f_final" in kwargs :
    #     if kwargs["f_final"] < kwargs["f_lower"] :
    #         # events["f_final"] = tmpWF.fcut(**events)
    #         events["f_final"] = np.array([2**(np.floor(np.log2(kwargs["f_lower"]))+1)])
    #     else :
    #         events["f_final"] = np.array([2**(np.floor(np.log2(kwargs["f_final"]))+1)])
    # else :
    #     events["f_final"] = 2**(np.floor(np.log2(tmpWF.fcut(**events)))+1)
    
    
    if "f_final" in kwargs :
        if kwargs["f_final"] < kwargs["f_lower"] :
            # events["f_final"] = tmpWF.fcut(**events)
            events["f_final"] = tmpWF.fcut(**events)
        else :
            events["f_final"] = np.array([kwargs["f_final"]])
    else :
        events["f_final"] = tmpWF.fcut(**events)

    delta_f = kwargs["delta_f"]
    fgrids = np.arange(0, events["f_final"]+delta_f/10., delta_f).reshape(-1, 1)  # 保留频率网格

    
    # 获取现象学振幅和相位
    xas_ampl = tmpWF.Ampl(fgrids, **events)
    # 2. 修正相位计算：还原旧代码规范，移除多余的-2psi和错误外层负号

    
    # if kwargs.get("choose_origin", False) :
    #     xas_phase = tmpWF.Phi_origin(fgrids, **events) + 2*phi0 + 2*(long_asc_nodes + mean_per_ano)
    # else :
    #     xas_phase = tmpWF.Phi(fgrids, **events) + 2*phi0 + 2*(long_asc_nodes + mean_per_ano)
    xas_phase = tmpWF.Phi(fgrids, **events) + 2*phi0 + 2*(long_asc_nodes + mean_per_ano)
    # Y22 = 2.0 * np.sqrt(5.0 / (64.0 * np.pi))，wf4py里面amp乘过了
    xas_phase = xas_phase.astype(np.float128)
    htilde22 = xas_ampl * np.exp(-1j * xas_phase)



    # 4. 修正极化因子计算：匹配旧代码规范，修复hc虚数符号
    cfac = np.cos(iota)
    pfac_hp = 0.5 * (1.0 + cfac ** 2)  # hp极化因子（不变）
    pfac_hc = -cfac                    # hc极化因子：移除多余的-1j（后续共轭会处理虚部）

    
    # hp = htilde22
    # hc = htilde22
    
    hp = np.conj(htilde22)
    hc = 1j*hp

    # # 5. 恢复Y22归一化：乘入hp/hc，修正振幅量级
    hp *= pfac_hp
    hc *= pfac_hc

    cos2psi = np.cos(2 * psi)
    sin2psi = np.sin(2 * psi)
    hp, hc = hp * cos2psi - hc * sin2psi, hp * sin2psi + hc * cos2psi

    x_phase = np.unwrap(np.angle(hp))
    pp.plot(np.sin(np.unwrap(x_phase - xas_phase))[1000:])

    # # 6. 还原旧代码核心规范：复共轭+符号修正（最关键，匹配观测数据极化相位）
    # hp = np.conj(hp)
    # hc = hp * np.exp(-1j*np.pi/2)

    # 7. 统一维度处理：与旧代码一致，显式转为一维，避免解析异常
    hp_flat = hp.reshape(1, -1)[0]
    hc_flat = hc.reshape(1, -1)[0]

    # print(hp_flat)

    hp_series = tps.FrequencySeries(hp_flat, delta_f=delta_f, epoch=_lal.LIGOTimeGPS(-16.0))
    hc_series = tps.FrequencySeries(hc_flat, delta_f=delta_f, epoch=_lal.LIGOTimeGPS(-16.0))
    # import matplotlib.pyplot as pp
    # pp.plot(hp_series.sample_frequencies, hp_series)
    return hp_series, hc_series
