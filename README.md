# 🚀 Interplanetary Mission Design to NEO Asteroid 22099

## 📝 Overview
This project was developed for the **Space Mission Analysis Laboratory** course at Politecnico di Milano. It details the preliminary design and trajectory optimization of an interplanetary mission to the Near-Earth Object (NEO) Asteroid 22099 (2000 EX106). 

The mission design is divided into three main scenarios based on the **Patched-Conics method**, assuming impulsive maneuvers:

* **Scenario 1 - Earth Parking Orbit:** Design and optimization of geocentric orbital transfers from a highly eccentric initial orbit to a designated parking orbit. Different strategies (standard sequence, bi-elliptic, and bitangent transfers) were compared to minimize the required $\Delta V$ and Time of Flight (TOF).
* **Scenario 2 - Interplanetary Transfer:** Optimization of the heliocentric transfer trajectory from Earth to the asteroid. A grid search combined with a Genetic Algorithm (GA) and `fmincon` was implemented to find the optimal transfer orbit parameters, balancing fuel consumption ($\Delta V$) and mission duration.
* **Scenario 3 - Escape and Capture Phases:** Design of the hyperbolic escape trajectory from Earth's Sphere of Influence (SOI) and the arrival hyperbolic trajectory for capture around the asteroid, ensuring velocity matching with the heliocentric transfer orbit.

## 📄 Read the Report
The full theoretical background, mathematical models, maneuver analysis, and final results are detailed in the project report (Note: The report is in Italian):
👉 **[Read the Full PDF Report Here](Report.pdf)**

## 🛠️ Built With
* **MATLAB** (For numerical optimization, orbital mechanics calculations, and 3D trajectory plotting)
* **LaTeX** (For documentation)

## 🚀 Key Results
The globally optimized mission profile requires a total $\Delta V$ of **8.71 km/s** and a total flight time of approximately **233 days**, proving the feasibility of the designed transfer architecture.
