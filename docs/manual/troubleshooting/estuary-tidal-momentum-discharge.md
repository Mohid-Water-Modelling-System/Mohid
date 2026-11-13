# Estuary tidal momentum discharge

## The problem

Regional scale models and large scale models can't reproduce accurately estuary dynamics and their influence in the larger scale circulation. However, relevant processes are present when a river plume is discharged into a larger sea or ocean. Especially the density driven effects such as salinity fronts combined with an extra momentum discharge. The extra momentum discharge is mainly due to the M2 tidal wave component on entering and leaving the estuary and creating a "toilet-flush" effect. A probably critical parameter that would regulate the type of estuary tidal momentum discharge would be, besides the estuary tidal amplitude volume, the estuary cross-section. If it's narrower, then the jet-like discharge ejects the freshwater plume at the surface farther in the ocean than if it is broader. The question arise: How can we reproduce the tidal wave momentum discharge of an Estuary when no Estuary is present in the model?

## The solution

Let 

$$
Q=Q(t)
$$

 be the water flux exchanged between the estuary and the coastal region. Let A be the surface of the estuary. Let H be the mean tidal deviation from neap-tide to spring-tide inside the estuary. Let 

$$
\\tau
$$

 be the M2 tidal period. Then there exists a phase of the M2 tide such that the water volume that circulates during half an M2 tidal period is solely inward, or solely outward, and is worth 

$$
V = A \\cdot H
$$

, the whole estuary volume change from neap-tide to spring-tide. Conversely, we can equate the instant water flux exchange



$$
V=\\int\_0^{\\tau/2} \\, Q(t) \\, dt
$$



We can further argument that the water flux is sinusoidal for the M2 tidal period. Hence the solution to the above equation would be something like



$$
Q(t) = Q\_0 \\, \\sin\\left(\\frac{2\\,\\pi}{\\tau}\\,\\left(t - t\_0\\right)\\right)
$$



Solving the system we obtain



$$
Q(t) = \\frac{V\\,\\pi}{\\tau}\\,\\sin\\left( \\frac{2\\,\\pi}{\\tau}\\,\\left(t - t\_0\\right) \\right)
$$
