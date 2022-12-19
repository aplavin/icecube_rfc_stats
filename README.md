How to reproduce p-values:
- Install Julia, download from https://julialang.org/downloads/
- Set the required parameters in `compute.jl`: additional error magnitudes and number of MC simulations
- Run `julia -tauto compute.jl`, takes on the order of 1 minute for 10^4 simulations
- The output would look like this:
```julia
┌ Info: 
│   c.add° = 0.0
└   p = 0.016979660406791865
┌ Info: 
│   c.add° = 0.45
└   p = 5.999880002399952e-5
```
The obtained p-values should be in agreement with the black line in Figure 2 of https://arxiv.org/abs/2211.09631.