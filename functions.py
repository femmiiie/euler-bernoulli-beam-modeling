from math import pi, sin

"""
All displacement functions take the same set of inputs for consistency,
even if the functions themselves don't make use of all of the inputs
"""

def unweighted(x, g, w, d, p, L):
  return -480 * w * d * g


def sinusoidal(x, g, w, d, p, L):
  f = unweighted(x, g, w, d, p, L)
  s = -(p * g * sin(pi * x / L))
  return f + s


def diver(x, g, w, d, p, L):
  f = unweighted(x, g, w, d, p, L)
  wght = 0
  if (x >= 1.8 and x <= 2):
    wght = -g * 70 / 0.2

  return f + wght




def y_correct(x, f, w, d, L, E, I):
  return (f / (24*E*I)) * x**2 * (x**2 - (4*L*x) + (6 * (L**2)))