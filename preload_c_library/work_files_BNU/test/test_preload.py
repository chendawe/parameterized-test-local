from pycbc import waveform as wf

import numpy as np

import matplotlib.pyplot as pp


f_lower = 15
delta_f = 1.0/64
f_final = 2048
f_len = int(np.ceil((f_final) / delta_f)) + 1

my_nongr_args = dict(
    approximant = "IMRPhenomXPHM",
    mass1 = 30.0,
    mass2 = 40.0,
    
    distance    = 100.0,
    inclination = 0.0,

    f_lower = f_lower,
    delta_f = delta_f,
    f_final = f_final,
    dchi3 = 0.3,

    mode_array = [[2,-2],[2,2]]
)

gr_args = dict(
    approximant = "IMRPhenomXPHM",
    mass1 = 30.0,
    mass2 = 40.0,
    
    distance    = 100.0,
    inclination = 0.0,

    f_lower = f_lower,
    delta_f = delta_f,
    f_final = f_final,

    mode_array = [[2,-2],[2,2]]
)


my_nongr_wf = wf.get_fd_waveform(**my_nongr_args)
gr_wf = wf.get_fd_waveform(**gr_args)

lal_nongr_wf = np.load('lal-nongr_waveform.npy')

pp.figure()
pp.plot(my_nongr_wf[0].sample_frequencies, np.abs(my_nongr_wf[0].data), label="my-non-GR")
pp.plot(my_nongr_wf[0].sample_frequencies, np.abs(lal_nongr_wf), label="LAL-non-GR")
pp.plot(my_nongr_wf[0].sample_frequencies, np.abs(gr_wf[0].data), label="GR")
pp.legend()
pp.savefig("./2-2_mode_Amplitude.png")

pp.figure()
pp.plot(my_nongr_wf[0].sample_frequencies, np.unwrap(np.angle(my_nongr_wf[0].data)), label="my-non-GR")
pp.plot(my_nongr_wf[0].sample_frequencies, np.unwrap(np.angle(lal_nongr_wf)), label="LAL-non-GR")
pp.plot(my_nongr_wf[0].sample_frequencies, np.unwrap(np.angle(gr_wf[0].data)), label="GR")
pp.legend()
pp.savefig("./2-2_mode_phase.png")

pp.figure()
pp.plot(my_nongr_wf[0].sample_frequencies, np.unwrap(np.angle(my_nongr_wf[0].data))-np.unwrap(np.angle(lal_nongr_wf)), label="my-non-GR minus LAL-non-GR")
pp.legend()
pp.savefig("./difference1.png")

pp.figure()
pp.plot(my_nongr_wf[0].sample_frequencies, np.unwrap(np.angle(my_nongr_wf[0].data))-np.unwrap(np.angle(gr_wf[0].data)), label="my-non-GR minus GR")
pp.legend()
pp.savefig("./difference2.png")