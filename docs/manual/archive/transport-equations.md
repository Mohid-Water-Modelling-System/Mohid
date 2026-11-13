# Transport equations

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

The variation rate of a quantity concentration **c** inside a material volume **vol** is equal to the flow that comes in minus the flow that comes out, transported by the velocity field **v**, through the volume's boundary **A**, plus the turbulent flow that comes in minus the turbulent flow that comes out (turbulent diffusion), through the same boundary **A**, plus the source terms **Sc** minus the sink terms **Sk**.

Below's the correspondent mathematical formulation:

$lim_{\Delta t \rightarrow 0} \frac{ \left( \int_{vol} c\, dvol \right)^{t+\Delta t} - \left( \int_{vol} c\, dvol \right)^t }{\Delta t} = - \int_{A} c \, \left( \overrightarrow{v} \cdot \overrightarrow{n} \right) dA \; - \; \int_{A} \left( -\nu \, \overrightarrow{\nabla c} \right)\cdot \overrightarrow{n} dA \; + \; \{ Sc - Sk \}$

You'll notice that on the left hand side of the equation, quantities are evaluated either at time **t**, either at time **t + dt**. However, at what time should the right hand side quantities be evaluated?

- **t**: Explicit numerical method. Low computational cost, but low precision and unstable.
- **t + dt**: Implicit numerical method. High computational cost and low precision, but stable.
- **1/2t + 1/2(t + dt)**: Semi-implicit numerical method. Highest computational cost but higher precision. More stable than explicit method, but less stable than implicit method.
