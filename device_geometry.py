"""Holds the device gemoetry parameters (Table 5), taken from Wu et al.,
>> A Predictive 3-D Source/Drain Resistance Compact Model and the Impact on 7 nm and Scaled FinFets<<, 2020, with interpolation for 4nm. 16nm is taken from PTM HP.
"""

node_names = [16, 7, 5, 4, 3]
GP = [64, 56, 48, 44, 41]
FP = [40, 30, 28, 24, 22]
GL = [20, 18, 16, 15, 14]
FH = [26, 35, 45, 50, 55]
FW = [12, 6.5, 6, 5.5, 5.5]
vdd = [0.85, 0.75, 0.7, 0.65, 0.65]
