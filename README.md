# Rigorous Integration
The repository contains code for experiments regarding rigorous integrating os systems with conservations laws using CAPD library [http://capd.ii.uj.edu.pl/].

## Introduction
Rigorous numerical methods are used when we have to deal with uncertainty of coefficients, parameters and initial values. Problems like that are common we use systems of differential equations to model physical phenomena. However, we can overcome this difficulty by enclosing the initial value in a box and performing all computations on intervals. As a result, we do not get a single value but rather the entire set of values containing the exact solution.
On the other hand, there is also a large group of systems with specific properties, such as Hamiltonian systems, which preserve energy (as long as they are autonomous). Moreover, their flow is a symplectic transformation. To the best of our knowledge, no published rigorous algorithm makes use of this additional information.

## Goal
We could adopt existing algorithms by incorporating constraints, such as a condition on symplecticity of the flow in these systems or a condition on energy preservation. This is to achieve more precise results in the case of rigorous methods.

## Experiment
Code was created for one specific system - Henon-Heiles system and we compared the original results with results using additional properties.

## Code
HenonHeiles.cpp contains the implementation of the algorithm and conducted experiments. File HigherOrderEnclosure.cpp contains the implementation of higher order enclosure. Files Var.h and Node.h are are required for the algorithm to work.
