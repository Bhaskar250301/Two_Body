# Newtonian Inspiral Engine

**A numerical simulation of binary black hole dynamics, demonstrating the transition from stable Keplerian orbits to radiative inspiral.**

![Inspiral Demo](results/binary_inspiral_inset.png)

## Project Overview
This repository contains a suite of Python simulations designed to model the orbital dynamics of binary black hole systems. While full General Relativity requires complex tensor calculus, the "chirp" dynamics of a merger can be effectively modeled by augmenting standard Newtonian mechanics with phenomenological dissipative terms.

This project was built to demonstrate:
1.  **Numerical Stability:** Implementation of a 4th-order Runge-Kutta (RK4) integrator.
2.  **Conservation Laws:** Validation of the integrator using stable orbits.
3.  **Radiative Physics:** Modeling energy loss due to gravitational wave emission using drag forces.

## Repository Structure

The project is divided into a "Control" experiment and the main "Inspiral" simulation.

### 1. The Control: `stable_orbit.py`
* **Physics:** Pure Newtonian Gravity ($F = Gm_1m_2/r^2$).
* **Purpose:** Validates the integrator by demonstrating **Energy Conservation**.
* **Result:** The Hamiltonian (Total Energy) remains constant to within numerical precision ($10^{-13}$) over long timescales, proving the stability of the RK4 implementation.
* **Output:** `results/stable_orbit.png`

### 2. The Engine: `inspiral_engine.py`
* **Physics:** Newtonian Gravity + Dissipative Radiation Reaction ($F_{drag} \propto v/r^3$).
* **Purpose:** Simulates the energy loss characteristic of a gravitational wave merger.
* **Key Feature:** The simulation handles high-eccentricity initial conditions, revealing a "staircase" energy loss pattern. This mimics the burst-like emission of gravitational waves during periastron passage (closest approach).
* **Output:** `results/binary_inspiral_inset.png`

## Visualizations

### Radiative Inspiral & Energy Loss
The plot below shows the orbital decay (left) and the total mechanical energy over time (right). Note the **non-linear energy drop**: as the black holes spiral closer, the orbital period decreases and the dissipative force increases, leading to a runaway "chirp" at the merger.

![Inspiral Analysis](results/binary_inspiral_inset.png)

### Integrator Validation (Stable Orbit)
Before modeling dissipation, the engine was tested on a conservative system. The flat energy curve confirms that the orbital precession is physical, not a numerical artifact.

![Stable Orbit](results/stable_orbit.png)

## Technical Implementation
* **Language:** Python 3 (`NumPy`, `Matplotlib`)
* **Integrator:** Runge-Kutta 4 (RK4) with fixed time-step.
* **Safety Mechanisms:** Includes numerical clamping to prevent singularity divergence at $r \to 0$ (The "Merger Plunge").

## Usage
To run the simulations and generate the plots:

```bash
# Run the Control
python stable_orbit.py

# Run the Inspiral Engine
python inspiral_engine.py
