# NBodySimulation

Minimum Requirements:

    Python 3.9.4

    --> numpy==1.20.2

    --> matplotlib==3.4.1

    --> scipy==1.7.1
    
    --> astropy==4.2.1 (optional)



**Synopsis:**

The codes in this repository are used to create a simulation of `N` bodies, each of which is mutually attracted to the other `N-1` bodies. Various methods can be used to solve the system of ODEs and time-step through the simulation; as of now, only the Velocity-Verlet and Euler-Cromer methods are deubgged and implemented. These methods employ symplectic integration, a consequence of which is the conservation of energy - this is especially important for orbits over large time-scales. While any bodies with mass can be simulated, I thought it would be neat to use the major bodies of our solar system as a test-case; one has the option of using `astropy` to get the position and velocity coordinates of each major body. The figure below depicts the barycentric reference frame of the 100-year simulation (using time-steps of 0.2 days each) of the eight planets and the Sun, starting from 1/1/1900.

![Solar System: 2-D position vectors (COM frame)](https://github.com/mikeysflix/NBodySimulation/blob/master/Figures/Solar_System-Velocity_Verlet_(Leap_Frog)-vector_position-dims_0_1-ref_Center_of_Mass.png?raw=true) 

One has the option of changing reference frames. The figure below depicts the same simulation as the image above, but from the reference frame of the Earth.

![Solar System: 2-D position vectors (Earth frame)](https://github.com/mikeysflix/NBodySimulation/blob/master/Figures/Solar_System-Velocity_Verlet_(Leap_Frog)-vector_position-dims_0_1-ref_Earth.png?raw=true) 


It is known that the mutual gravitation between multiple bodies results in the perturbation of their orbits. These perturbations can accumulate in time to affect orbits in observable ways. An example of this can be seen by observing the tug-of-war between Neptune and Uranus, which was used to predict the existence of Neptune before it was observed. The figure below depicts the phase angle of Uranus from the reference frame of Earth by comparing the orbits with and without Neptune.

![Solar System: Uranus Phase X-Neptune (Earth frame)](https://github.com/mikeysflix/NBodySimulation/blob/master/Figures/Solar_System_(with_and_without_Neptune)-Velocity_Verlet_(Leap_Frog)_AND_Velocity_Verlet_(Leap_Frog)-offset_phase_position-ref_Earth.png?raw=true) 
  

One can use the vector positions of each body to find the scalar distance of each body from the barycenter. The figure below depicts these scalar distances.

![Solar System: 2-D position scalars (COM frame)](https://github.com/mikeysflix/NBodySimulation/blob/master/Figures/Solar_System-Velocity_Verlet_(Leap_Frog)-scalar_velocity-ref_Center_of_Mass.png?raw=true) 

One can use these scalar distances to compute the approximate period of each orbit. There are many methods that can be used for this purpose; as of now, the methods employed include boundary-crossings, auto-correlation, Fourier transform, and the power spectral density (Fourier transform of the auto-correlation). It is worth mentioning that auto-correlation tends to do worse for a body when multiple periods are present, and that the fourier methods are more robust (especially for larger `N`). The table below shows the estimated and true values for orbital periodicity, along with the relative error.

![Solar System: Period Table](https://github.com/mikeysflix/NBodySimulation/blob/master/Figures/Solar_System-Velocity_Verlet_(Leap_Frog)-Table-Orbital_Period.png?raw=true) 


**Things To Add:**

1) Compare simulation trajectory to the actual trajectory for each body.

2) Play with exotic orbits (figure-eight, horse-shoe, etc).

3) Implement non-symplectic methods with adaptive time-stepping.

4) Apply rotation to reference frames and find Larange points.

5) Launch rocket from one body to another.



#
