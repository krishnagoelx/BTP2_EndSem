
## Analysis of Turbulent Velocity Wavenumber Spectrum

### 1. 2D Contour of Φ(kₓ,kᵧ)
The 2D contour plot shows that most of the turbulent energy is concentrated at low wavenumbers 
(large eddies), forming a bell-shaped distribution centered at the origin. This is characteristic 
of homogeneous, isotropic turbulence. The steep flanks in the plot demonstrate the inertial-range 
decay with the -7/3 roll-off in the denominator of the von Karman spectrum. This steep decline in 
energy as wavenumber increases aligns with Kolmogorov's theory of energy cascade in turbulence.

### 2. Radial Cut and Slope Verification
The log-log plot of spectrum magnitude versus wavenumber confirms that the von Karman model follows 
a k^(-7/3) slope in the inertial range. This corresponds to the -5/3 power law in physical space 
predicted by Kolmogorov's theory, but appears as -7/3 in the von Karman spectrum due to the 
formulation. At low wavenumbers (energy-containing range), the spectrum deviates from this slope 
as these large scales contain the most energy. At high wavenumbers, the spectrum tends to fall off 
even more rapidly, representing the dissipation range where viscous effects become important.

### 3. Effect of Integral Length Scale
Changing the integral length scale Λ shifts the peak of the spectrum along the wavenumber axis.
A larger Λ shifts the peak toward smaller wavenumbers, indicating that more energy is contained in 
larger turbulent structures. This parameter is crucial for aeroacoustic modeling as it defines the 
scale of the most energetic eddies interacting with the airfoil. In practical terms, increasing the 
integral length scale typically leads to higher noise predictions at low frequencies and lower noise 
at high frequencies, as the energy is redistributed across the spectrum.

### 4. Comparison of von Kármán and Liepmann Models
The two turbulence models show similar behavior at low wavenumbers but diverge significantly at high 
wavenumbers. The von Kármán model has the theoretically correct -7/3 roll-off in the inertial range, 
making it more physically accurate according to Kolmogorov theory. The Liepmann model, while 
mathematically simpler, exhibits a different high-wavenumber behavior. For aeroacoustic predictions, 
the choice between these models can impact high-frequency noise estimates, particularly when small-scale 
turbulent structures interact with airfoil leading edges.
