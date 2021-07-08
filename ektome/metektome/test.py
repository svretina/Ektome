#!/usr/bin/env python3


import ektome.metektome.error as e

folder = "vanilla_q1000_b6040_px10_py1-0.41_pz10_sx10.9_sy10.1_sz10_px20_py20.41_pz20_sx20_sy20_sz2-0.1_ex_r500"
err = e.Error(folder, dim=3)
info = err.error_report()
print(info)
