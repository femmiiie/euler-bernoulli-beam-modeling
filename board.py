from scipy.sparse import lil_matrix
from scipy.sparse.linalg import onenormest, inv, spsolve
import numpy as np

class Board:
  """ Representation of the board for use in main

  Attributes:
    disp_func: weight function at a position on the board x
    correct: analytical solution function at a position on the board x
    g: gravitational constant
    w: width
    d: thickness
    L: length
    E: Young Modulus
    I: area moment of inertia
    inorm: lambda for the inf-norm of a vector
    n: number of grid steps for the board
  """
  def __init__(self, f, y, g, w, d, p, L, E, n=None):
    self.disp_func = f
    self.correct = y
    self.g = g
    self.w = w
    self.d = d
    self.L = L
    self.p = p
    self.E = E
    self.I = w * (d**3) / 12

    self.inorm = lambda i : max(abs(min(i)), max(i))

    if n is not None:
      self.set_n(n)


  def set_n(self, n):
    """Defines the number of grid steps to use in calculations"""

    self.n = n
    self.h = self.L / n
    self.b = (self.h**4 / (self.E * self.I)) * np.array([self.f(x) for x in self.get_x_vals()])
    self.create_structure_matrix(n)


  def get_x_vals(self):
    """Returns the set of n equally spaced x values along the board"""

    if self.n is None:
      return None
    
    return [self.h*i for i in range(1, self.n + 1)]


  def get_correct(self):
    """Returns the analytical solution of the board without weight"""

    return [self.correct(
      self.h*i, 
      self.f(self.h*i), 
      self.w, 
      self.d, 
      self.L, 
      self.E, 
      self.I) 
    for i in range(1, self.n + 1)]


  def get_fe(self):
    """Gets the forward error of the solution compared to correct"""

    x_c = self.solve()
    correct = self.get_correct()
    return [x_c[i] - correct[i] for i in range(self.n)]


  def solve(self, n=None):
    """
    
    """

    if n is not None:
      self.set_n(n)

    return spsolve(self.A.tocsr(), self.b)
    
  
  def f(self, x):
    """Wrapper function for disp_func"""

    return self.disp_func(
      x, 
      self.g, 
      self.w, 
      self.d,
      self.p, 
      self.L
    )


  def pow_2_chart(self, k):
    """Iterates through power of 2 grid steps and collects metrics on them
    
    Args:
      k: Amount of iterations to perform
    """

    n = []
    fe = []
    cond = []
    x = []

    for i in range(1, k + 1):      
      x_c = self.solve(10 * 2**i)

      n.append(10 * 2**i)
      x.append(x_c[-1])
      fe.append(self.inorm(self.get_fe()))

      #there is no known way of quickly finding the condition number of a large sparse matrix in python to my knowledge
      #if runtime is a concern, comment out this line and instead just append 1 to cond
      cond.append(onenormest(inv(self.A.tocsc()).dot(self.A.tocsc())))

    print("--------------------------------------------------")
    print("|   n   |     x_n      |      fe      |  cond    |")
    print("--------------------------------------------------")
    for i in range(len(n)):
      print(f"| {n[i]:5d} | {x[i]:.5e} | {fe[i]:1e} | {cond[i]:.2e} |")
    print("--------------------------------------------------")


  def create_structure_matrix(self, n):
    """Creates an nxn structure matrix defined in the paper"""

    A = lil_matrix((n, n))  
    
    #first row boundary case
    A[0, 0] = 16
    A[0, 1] = -9
    A[0, 2] = 8/3
    A[0, 3] = -1/4  
  
    for i in range(1, n-2): #generalized case
      A[i, i-2] = 1
      A[i, i-1] = -4
      A[i, i]   = 6
      A[i, i+1] = -4
      A[i, i+2] = 1 
  
    A[1, -1] = 0 #remove wraparound in second row from loop 
 
    #last 2 rows boundary case
    A[-2, -4:] = [16/17, -60/17, 72/17, -28/17]
    A[-1, -4:] = [-12/17, 96/17, -156/17, 72/17]

    self.A = A.tocsr() #store as CSR for space efficiency