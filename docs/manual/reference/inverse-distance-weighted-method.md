# Inverse-distance weighted method

In this method of interpolation the value of a property at a target location is determined as the weighted sum of the property value at the sampled locations, being the weight of each sample value a power of the inverse of the distance between the target location and the sample location.

The value of the inverse-distance weighting power controls the region of influence of each sampled location: as its value increase the region of influence decreases until it becomes the area which is closer to the target location than to any other location, a null value causes the method of interpolation to be equivalent to averaging the sample values.

References:

Lam, N. S., 1983, "Spatial interpolation methods review", *The American* Cartographer*, 10, pp. 129-149.*
