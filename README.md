# Modeling Euler-Bernoulli Beams

## Table of Contents
1. [Overview](#overview)
2. [Introduction](#introduction)
3. [Usage](#usage)
4. [Dependencies](#dependencies)
5. [Authors](#authors)
6. [License](#license)

# Overview

This repository accompanies a paper created for MAD4401 - Numerical Analysis titled: Modeling Euler-Bernoulli Beams. Due to the University of Florida's Student Honor Code, I am unable to share more than the introduction of the original paper publicly. If you are a recruiter or potential employer, please contact me and I can provide said paper.


# Introduction

The Euler-Bernoulli beam equation is used in structural mechanics to model the displacement of 
an elastic beam with a varied load on top of it. It is used in various engineering disciplines, where 
modeling how beams bend under several types of loads is important to safely design buildings, bridges, 
and other types of structures and projects. 

This paper looks at a specific case of the beam equation, being the Clamped-Free Beam, where 
one end of the beam is fixed in place, while the other is allowed to move freely. To model this system, we 
use the Euler–Bernoulli beam equation in its fourth-order differential form. Since analytical solutions are 
limited to idealized scenarios, we apply finite difference methods to discretize and solve the equation 
numerically. This approach allows us to approximate the beam’s displacement under a constant 
distributed load and investigate how solution accuracy and computational behavior change with various 
levels of subdivision. 

The goal of this paper is to evaluate the accuracy, convergence, and numerical stability of the 
finite difference method when applied to both the clamped-free beam as well as the clamped-clamped 
beam. By comparing the numerical results to the known analytical solution and analyzing the error and 
condition number of the system, we gain insight into the effectiveness and limitations of the numerical 
approach.


# Usage

Run the program using `python main.py`
The computational steps within main.py were not designed inherently to be run all at once.
For ease of use, comment out the steps not being used at the bottom of main.py and then run.

# Dependencies

- [matplotlib](https://matplotlib.org/)


# License

This repository is shared publically to showcase the abilities of the author and add to their portfolio. Any usage of this repository for purposes of academic dishonesty is strictly prohibited. The codebase is partially obfuscated in an effort to prevent this from occurring. By continuing to use this repository, the user assumes all liability for any recourse, positive or negative, brought upon them from it.