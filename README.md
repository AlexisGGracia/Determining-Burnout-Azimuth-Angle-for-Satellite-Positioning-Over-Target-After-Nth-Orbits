# Determination of the Burnout Azimuth Angle to Position a Satellite Over a Selected Location on Earth After Nth Orbital Periods in MATLAB

# Introduction

This project utilizes Ms. Katherine Johnson's NASA technical notes to determine the burnout azimuth angle required for a spacecraft to position itself over a desired location on Earth after a specific number of orbital periods accounting for J2 pertubations. The simulation was implemented in MATLAB and is designed for satellite mission planning. It calculates the orbital elements necessary to position the spacecraft for applications such as reconnaissance, Earth observation, or reentry positioning.

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

     
         
2. **Using Katherine Johnson and Skopinski NASA technical notes, determine the burnout azimuth of the spacecraft.**

3. **Determine the orbital elements of the satellite at burnout** to achieve the goal of passing over a desired position after a precise number of revolutions.

4. **Propagate the orbit based on two-body dynamics**, including J2 perturbations.

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

- **Python 3.x**: Ensure Python 3.x is installed. You can download it from [here](https://www.python.org/downloads/).
- **ArcGIS Pro**: You will need ArcGIS Pro to run the script and use its mapping features. You can get it from [ArcGIS Pro](https://www.esri.com/en-us/arcgis/products/arcgis-pro).
- **ArcPy**: This Python library is used for automating GIS workflows in ArcGIS Pro. It comes pre-installed with ArcGIS Pro.
- **Additional Python Libraries**:
  - `numpy`
  - `pandas`
  - `re`
  - `math`
  - `datetime`
  - `os`

You can install the Python libraries using `pip`:

`pip install numpy matplotlib pandas`
# Visual Results
In this section, I have included images that demonstrates the final product of my code. The green figures were made manually which takes about 2-5 minutes to generate. The figure with a distinct color was made using my code automatically which takes less than 3 seconds to generate.

## Green Color Figure: 
    Corresponds to a manually drawn parcel (2-5 minutes to generate)

## Distinct Color Figure: 
    Corresponds to the parcel generated by my code (less than 3 seconds to generate)

![SC1](https://github.com/user-attachments/assets/3c0bfd72-6d8e-4b8f-aa28-35aa16d593ab)
![SC2](https://github.com/user-attachments/assets/c3d17219-bb17-47b0-8314-94063c762bc5)
![SC3](https://github.com/user-attachments/assets/e114a0a3-1333-4b14-ac6f-54377539b38f)
![CS4](https://github.com/user-attachments/assets/2dfefa65-8b77-49a4-9336-92de5c340583)

## Project Usage

This project is used by the following company:

- **Office of the County Engineer**
    - GIS Team

## Contributor

This project was developed and maintained by:

- **Alexis Gracia**  
  [GitHub Profile](https://github.com/AlexisGGracia)  
  [Email](mailto:agg3455@my.utexas.edu)
