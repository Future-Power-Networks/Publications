# Supplementary Files for the Large-Signal Stability of Power Systems with Mixtures of GFL, GFM and GSP Inverters
The provided three folders correspond to three different inverter combination scenarios:

- `GFL-GFL`: Both inverters are Grid-Following (GFL).
- `GFL-GFMorGSP`: One inverter is Grid-Following (GFL), and the other is either Grid-Forming (GFM) or Grid-Supporting (GSP).
- `GFM-GSPorGFM`: One inverter is Grid-Forming (GFM), and the other inverter is either Grid-Following (GFL) or Grid-Supporting (GSP).

## Instructions for Running the Programs

1. **Run** `**main.m**` **first:**
   All inverter parameters are configured within this file. The file performs the following functions:

   - Plotting Domain of Attraction (DOA)
   - Analyzing manifold structures
   - Calculating the Stability Radius
   - Executing reduced-order model simulations

   Key parameters within `main.m`:

   - The variable `system`:
     - `"GFL"`: Both inverters are GFL.
     - `"GFM"`: Both inverters are GFM.
     - `"STATCOM"`: IBR1 is GFL, IBR2 is GSP.
     - `"GFMGFL"`: IBR1 is GFL, IBR2 is GFM.
     - `"GFMGFLV"`: IBR1 is GFM, IBR2 is GSP.

2. **Next, run** `**para.m**` **from the Simulink folder** to initialize and assign parameters to the Simulink model.

3. **Finally, run** `**plot_sim.m**` to compare simulation results with theoretical analysis.

