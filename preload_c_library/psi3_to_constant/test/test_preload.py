from pycbc import waveform as wf

import numpy as np

import matplotlib.pyplot as pp


f_lower = 15
delta_f = 1.0/64
f_final = 2048
f_len = int(np.ceil((f_final) / delta_f)) + 1

test_args = dict(
    approximant = "IMRPhenomXPHM",
    mass1 = 30.0,
    mass2 = 40.0,
    
    distance    = 100.0,
    inclination = 0.0,

    f_lower = f_lower,
    delta_f = delta_f,
    f_final = f_final,
    
)


test_wf = wf.get_fd_waveform(**test_args)

pp.plot(test_wf[0].sample_frequencies, np.abs(test_wf[0]))

pp.plot(test_wf[0].sample_frequencies, np.unwrap(np.angle(test_wf[0])))