.LIB 4NM_FINFET_HP

** ASAP FF models v1.0 8/3/16

** Hspice modelcard
.model nmos nmos level = 72 

*NOTE: Paste here the slvt NMOS model from ASAP7_PDKandLIB_v1p6/pdk/asap7PDK_r1p6/models/hspice/7nm_FF.pm
*      then change hfin to 5e-008 and tfin to 5.5e-009 (Table 5 in the ISFPGA'21 paper). 

** Hspice modelcard
.model pmos pmos level = 72 

*NOTE: Paste the PMOS model and apply the same changes as for the NMOS.

.ENDL
