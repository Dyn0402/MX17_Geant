# Draft: questions for Irakli (DREAM bench data, FRESH_LOOK §19)

For Dylan to send. Nothing sent yet.

---

Hi Irakli,

Thanks again for the NSW and TPOT bench measurements. They already helped: the TPOT DC gain
(18.05 mV/fC at 0 pF, 180 ns) suggests our simulated chain is ~10 % high in ADC/fC.

A few questions so we can use the NSW AC numbers as well:

1. In the NSW "AC" measurement, was the coupling a series capacitor added on the test board,
   or the spark-protection network of the protected FEU variant? If it was a discrete capacitor,
   what value? Our fit to your gain-vs-C points gives an effective C_c ≈ 510 pF.
2. For the protected FEU (the version built for non-resistive Micromegas): do you have the
   protection schematic or BOM, in particular the series capacitor/resistor values and the
   diodes? Our cosmic bench ran on protected FEUs, and the X and Y strip capacitances differ,
   so a series C there could change the X/Y charge balance.
3. Has the 200 fC range been measured on the DC (TPOT) setup? The TPOT slides only show
   100 fC (the NSW sheets have both, ~x1.95 apart), and our data use 200 fC.
4. Is the 18.05 mV/fC with common mode removed, or raw? And at which DREAM peaking-time
   register code (we assume code 2, 180 ns nominal)?

Thanks,
Dylan
