# functional_controllability
 
This toolkit is for computing the controllability measurements, minimum control energy, and kernel ridge prediction for behavior tasks.
## Test environment
1. Matlab 2018a
2. Python version 3.6.7
3. sklearn version 0.19.1
## minimum control energy

Use Energy Efficiency/pipelineCompareStaticDynamicEnergy.m to compute the minimum static and dynamic control energy.
`[Es, Ed, delta_s] = pipelineCompareStaticDynamicEnergy(time_series, control_node, window_size)`

Efficiency/pipelineCompareDynamicEnergyAfterShuffle.m to compute the dynamic control energy before and after shuffle.
`Er = pipelineCompareDynamicEnergyAfterShuffle(time_series, control_node, window_size, shuffle_times)`
