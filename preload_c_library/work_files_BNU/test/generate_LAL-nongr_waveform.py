import numpy as np
from pycbc import waveform as wf

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
    dchi3 = 0.3,
    mode_array = [[2,-2],[2,2]]
)


test_wf = wf.get_fd_waveform(**test_args)


# 创建示例复数数组
original_array = test_wf[0].data

# 保存到文件（自动保留复数类型和形状）
np.save('lal-nongr_waveform.npy', original_array)

