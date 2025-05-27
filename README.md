# LLP-Solver - Dual-Method
This C++ project is a Linear Programming Problem (LPP) solver that applies the Dual Simplex Method to find the optimal solution to minimization problems. It is tailored for problems where the number of constraints equals the number of variables, making it ideal for dual formulations.

## Features :
Accepts user input for:

Number of variables and constraints

Coefficients of constraints (LHS and RHS)

Objective function (minimization)

Implements the Dual Simplex Method

Handles:

Degeneracy (multiple equal ratios)

Multiple optimal solutions

Unbounded solutions

Displays:

Each tableau before and after pivoting

Ratios and pivot elements

Final values of basic variables

Final minimized objective value (Z)

## What’s inside :
Matrix manipulation (transposition, reconstruction)

Pivot operations and ratio testing

Basic variable identification

Detection of degeneracy and multiple optimal solutions

## Constraints :
Only supports problems where the number of variables equals the number of constraints (i.e., dual form applicable).

Designed for educational and demonstration purposes in operations research and optimization.

## How to Use :
Compile the program using a C++ compiler (e.g., g++).

Run the executable and follow prompts to enter:

Coefficients for constraints

Right-hand side values (RHS)

Objective function coefficients

View step-by-step simplex tableau transitions and the final result.
