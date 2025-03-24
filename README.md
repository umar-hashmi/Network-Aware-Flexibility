# Network aware energy market participation of industrial flexibility

The key contributions of this work are as follows:

•	Two-Stage Network-Aware (NA) Market Participation: Introduces a novel day-ahead (DA) and real-time (RT) framework that leverages industrial flexibility for energy markets while preventing network incidents (e.g., voltage violations, congestion).

•	Dynamic Operating Envelopes (DOEs) for Safe Flexibility: Proposes DA DOEs that account for load fluctuations, ensuring flexibility commitments do not compromise power quality or grid stability.

•	Myopic Real-Time Control for Grid Stability: Implements a fast RT control loop (minute-level updates) using active and reactive power adjustments to minimize deviations from DA commitments while maintaining voltage stability.

•	Trade-off Analysis Between Profit and Grid Health: Demonstrates through case studies that NA participation reduces profits by up to 40% but eliminates local voltage violations, while also analyzing converter sizing impacts on market gains and voltage correction.

## Code Description

• Case 0: No flexibility considered 

• Case 1: Industrial flexibility energy market participation with no network consideration, 

• Case 2: Industrial flexibility energy market participation with DA network consideration, 

– DA uncertainty is modelled using 100 Monte Carlo scenarios generated using 30% forecast error based on the scenario generation method.

• Case 3: Industrial flexibility energy market participation with RT network consideration, and 

• Case 4: Industrial flexibility energy market participation with DA and RT network consideration.

## KPI considered
KPI1: Profit from market participation

KPI2: Battery number of cycles of operation

KPI3: Cumulative voltage correction (CVC)

KPI4: Voltage correction index (VCI)

KPI5: CO2 emission reduction

KPI6: Power factor variation

KPI7: Statistical attributes of aggregated load

KPI8: Mean Converter Usage (MCU)

## Citation
@article{hashmi2025NetworkAware,

  title={Network aware energy market participation of industrial flexibility},
  
  author={Md Umar Hashmi, Manna Carlo and Dirk Van Hertem},
  
  journal={Under review},
  
  year={2025}
  
}


