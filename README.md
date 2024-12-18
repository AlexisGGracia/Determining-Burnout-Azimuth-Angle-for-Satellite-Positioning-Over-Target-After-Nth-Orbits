# Determination of the Burnout Azimuth Angle to Position a Satellite Over a Selected Location on Earth After Nth Orbital Periods in MATLAB

# Introduction

This project utilizes Ms. Katherine Johnson's NASA technical notes from 1960 to determine the burnout azimuth angle required for a spacecraft to position itself over a desired location on Earth after a specific number of orbital periods accounting for J2 pertubations. The simulation was implemented in MATLAB and is designed for satellite mission planning. It calculates the orbital elements necessary to position the spacecraft for applications such as reconnaissance, Earth observation, or reentry positioning.

Who is it for? This project was developed at The University of Texas at Austin as part of my Final Project for my Spacecraft Dynamics course taken during Fall 2023

# Code Breakdown and Workflow

This section provides a detailed explanation of the logic and steps involved in the code.

1. **Defines parameters such as semi-major axis, eccentricity, Greenwich Sidereal Meridian Time (GMST) at burnout, etc.**

   - Orbital elements such as semi-major axis, eccentricity are determined based on mission constraints such as altitude of orbiting spacecraft and desired orbit type.
   - Parameters defined are as follows:

   ```MATLAB
   % Defining parameters
   format longG
   a = 6500;                       % Semi-major axis (km)
   e = 0.001;                      % Orbit's eccentricity
   theta1 = 20*pi/180;             % Initial true anomaly at burnout
   GSMT = 0;                       % Greenwich Sidereal Meridian time at burnout
   mu_earth = 398600.4415;         % Gravitational parameter of Earth (km^3/s^2)
   R_earth = 6378.1363;            % Radius of Earth (km)
   J2 = 0.0010826267;              % J2 perturbations due to Earth's obliqueness
   w_earth = 2*pi / 86164;         % Earth's rotation rate
   n = 10;
   iterations = 4;
   i = 1;                          % Counting variable
   f = 3.353*10^-3;               % Earth's obliqueness constant
   height = 0;

     
         
2. **Using Katherine Johnson and Skopinski NASA technical notes, determine the burnout azimuth of the spacecraft using your preferred iterative method.**

3. **Determine the orbital elements of the satellite at burnout** to achieve the goal of passing over a desired position after a precise number of revolutions.**

4. **Propagate the orbit using ODE45 based on two-body dynamics, including J2 perturbations.**

5. **Determine the accuracy of the solution by:**
   - (a) **Generating groundtracks for the satellite**, including coastlines and the desired location on the map.
   - (b) **Plotting the satellite elevation relative to the desired location** over the time period when it should pass overhead.


### Important Notes

- **J2 Perturbations**: In this project, J2 perturbations (due to Earth's oblateness) are only considered during the **first iteration** of the burnout azimuth calculation.
                        This is because the primary goal of this calculation is to determine the **initial azimuth angle at burnout** required to position the spacecraft
                        over a target location immediately after burnout.
                        
  
- **First Iteration**: During the first iteration, J2 perturbations are included to refine the orbital parameters (such as inclination and node) at the moment of burnout.
                       These perturbations slightly affect the spacecraft’s trajectory and need to be accounted for to ensure accurate initial positioning.

- **Subsequent Iterations**: For later iterations, J2 perturbations are **ignored**. This simplification is made because the subsequent calculations focus primarily on
                             placing the satellite in the desired position **immediately after burnout**. The long-term effects of J2 on the orbit (such as orbital precession)
                             are considered minimal for this immediate mission goal.

- **Why This Approach?**: The neglect of J2 in later iterations is a practical choice, as the impact of J2 on short-term positioning is minimal compared to the complexity
                          it would add to the calculations. For long-term mission planning, J2 and other perturbations should be modeled, but for the purposes of this burnout
                          azimuth determination, they do not significantly affect the accuracy of the placement.




  
## Requirements

To run this project, you will need the following software and tools:

- **MATLAB**: Ensure that MATLAB is installed. You can download it from [here](https://www.mathworks.com/products/matlab.html).
- **MATLAB Toolboxes**: The following toolboxes are recommended for running the simulations:
  - Aerospace Toolbox
  - Mapping Toolbox
  - MATLAB functions which are included as separate files under this project

### Additional Information:
The project uses MATLAB's built-in functions and toolboxes to perform orbital simulations, including two-body propagation, J2 perturbations, and groundtrack generation.
If you cannot acess the aerospace and mapping toolbox for free, I have included the workspace to upload a 2D map to plot groundtracks. You just need to download it into your workspace

You can install necessary MATLAB packages directly through the MATLAB interface.


# Visual Results
In this section, I have included the groundtracks and an elevation plot relative to my desire location that demonstrates the accuracy of this model. Additionally, I have included the groundtracks and elevation plot replicating the results from the technical case presented by Ms. Katherine Johnson

## Groundtracks (Austin, Texas Case): 

- **Groundtracks**: Starting location is along the coast of california and final destination is Austin, Texas. My mission was designed such that it completes this after 10th orbital periods.
![Austin_Grountrack_Pic](https://github.com/user-attachments/assets/95396d28-e579-4f98-964d-ea55ac6b8ee7)
    
## Elevation vs Time plot (Austin, Texas case): 

- **Elevation Plot**: The elevation plot versus time plot helps track how the altitude of a spacecraft changes over time, providing valuable insight into its orbital path and visibility from a ground station or observation point. 
   - Since the final destination is Austin, Texas (right above the ground station) it shows that the final elevation is 90 degrees

![Austin_Texas_ELEVATION](https://github.com/user-attachments/assets/dbe41be2-75e8-48ee-8ffe-a421fe4186d0)

## Groundtracks (Technical Case): 
![Technical_Case_Groundtrack](https://github.com/user-attachments/assets/7b3b37a2-e049-4298-bf66-f68f1830c0ad)

## Elevation vs Time plot (Techncial Case): 
![TechnicalCase_elevation](https://github.com/user-attachments/assets/b12f5b10-09d6-4d5a-a502-db3e17b99f3c)

## Citations

H. Skopinski and K. G. Johnson, *"Determination of Azimuth Angle at Burnout for Placing a Satellite Over a Selected Earth Position"*, NASA Tech Note D-233, September 1960. [NASA Technical Report](https://ntrs.nasa.gov/citations/19980227091)


## Contributor

This project was developed and maintained by:

- **Alexis Gracia**  
  [GitHub Profile](https://github.com/AlexisGGracia)  
  [Email](mailto:agg3455@my.utexas.edu)
