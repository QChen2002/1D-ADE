import numpy as np
import matplotlib.pyplot as plt 

# Define the elements
class Elements:
    def __init__(self, dim, N1, N2, N3):
        if dim == 1:
            # The label of node 1
            self.N1 = N1
            # The label of node 2
            self.N2 = N2
        elif dim == 2:
            # The label of node 1
            self.N1 = N1
            # The label of node 2
            self.N2 = N2
            # The label of node 2
            self.N3 = N3

# Define the problem.
class Problem:
    def __init__(self, Pe, L, bc1, bc2, base_function):
        # Set the cell Peclet number
        self.Pe = Pe
        # Set the computational domain
        self.L = L
        # Compute h by h = Pe/20.
        self.h = Pe/20
        # The Dirichlet boundary condition at x = 0
        self.bc1 = bc1
        # The Dirichlet boundary condition at x = 1
        self.bc2 = bc2
        # Set the order of base function
        self.base_function = base_function
        # Initialize the global stiffness matrix. 
        # The number of nodes = (L/h)*order+1.
        self.N = int(self.base_function*self.L/self.h+1)
        # The number of elements = L/h
        self.N_ele = int(self.L/self.h)
        self.K =  np.zeros((self.N, self.N))
        # Initialize the force matrix.
        self.f = np.zeros((self.N,1))
        # Set up the elements
        self.get_elements()
        # Set up the grids
        self.grids = np.linspace(0,L,self.N)
        
    
    # Set up the elements
    def get_elements(self):
        self.elements = []
        for i in range(self.N_ele):
            # Create a new element
            self.elements.append(Elements(self.base_function, self.base_function*i+1,self.base_function*i+2,self.base_function*i+3))

    # Get the global stiffness matrix
    def get_global_K(self):
        for element in self.elements:
            # Compute the component of element stiffness matrix
            if self.base_function == 1:
                # Use linear base function
                Ke_11 = -1/2 + 0.05/self.h
                Ke_12 = 1/2 - 0.05/self.h
                Ke_21 = -1/2 - 0.05/self.h
                Ke_22 = 1/2 + 0.05/self.h
            elif self.base_function == 2:
                # Use quadratic base function
                Ke_11 = -1/2 + 0.05*7/(3*self.h)
                Ke_12 = -2/3 - 0.05*8/(3*self.h)
                Ke_13 = 1/6 + 0.05/(3*self.h)
                Ke_21 = 2/3 - 0.05*8/(3*self.h)
                Ke_22 = 0.05*16/(3*self.h)
                Ke_23 = -2/3 - 0.05*8/(3*self.h)
                Ke_31 = -1/6 + 0.05/(3*self.h)
                Ke_32 = 2/3 - 0.05*8/(3*self.h) 
                Ke_33 = 1/2 + 0.05*7/(3*self.h)


            # Put the element stiffness matrix into global stiffness matrix
            if self.base_function == 1:
                self.K[element.N1-1][element.N1-1] += Ke_11
                self.K[element.N1-1][element.N2-1] += Ke_12
                self.K[element.N2-1][element.N1-1] += Ke_21
                self.K[element.N2-1][element.N2-1] += Ke_22
            
            elif self.base_function == 2:
                self.K[element.N1-1][element.N1-1] += Ke_11
                self.K[element.N1-1][element.N2-1] += Ke_12
                self.K[element.N1-1][element.N3-1] += Ke_13
                self.K[element.N2-1][element.N1-1] += Ke_21
                self.K[element.N2-1][element.N2-1] += Ke_22
                self.K[element.N2-1][element.N3-1] += Ke_23
                self.K[element.N2-1][element.N1-1] += Ke_31
                self.K[element.N2-1][element.N2-1] += Ke_32
                self.K[element.N2-1][element.N3-1] += Ke_33

    # Get the global force matrix
    def get_global_f(self):
        self.f[0][0] = -self.bc1
        self.f[-1][0] = self.bc2

    # Reduce the global stiffness matrix and force matrix
    def reduce_K_f(self):
        for i in range(self.N):
            self.f[i][0] += (-self.K[i][0]*self.bc1-self.K[i][-1]*self.bc2)
        self.K = np.delete(self.K,i,axis=0)
        self.K = np.delete(self.K,0,axis=0)
        self.K = np.delete(self.K,i,axis=1)
        self.K = np.delete(self.K,0,axis=1)
        self.f = np.delete(self.f,i,axis=0)
        self.f = np.delete(self.f,0,axis=0)
         
    
    # Solve the equation.
    def solver(self):
        # Assemble the global stiffness matrix.
        self.get_global_K()
        # Assemble the global force matrix
        self.get_global_f()
        # Reduce the global stiffness matrix and force matrix.
        self.reduce_K_f()
        # Solve the nodal equation.
        self.x = np.linalg.inv(self.K) @ self.f
        self.x = np.insert(self.x,int(1/self.h)-1,self.bc2,axis=0)
        self.x = np.insert(self.x,0,self.bc1,axis=0)

# Test 4 different Pe
problem1 = Problem(1,1,0,10,1)
problem2 = Problem(1,1,0,10,2)
problem1.solver()
problem2.solver()
# Compare with the analytical solution
grids_analytical = np.linspace(0,problem1.L,41)
sol = []
for i in grids_analytical:
    sol.append(10/(np.exp(20)-1)*np.exp(20*i)-10/(np.exp(20)-10))
# Plot
plt.plot(problem1.grids,problem1.x,label='Numerical Solution for Pe = 1.25, linear basis')
plt.plot(problem2.grids,problem2.x,label='Numerical Solution for Pe = 1.25, quadratic basis')
plt.plot(grids_analytical,sol,label='Analytical Solution')
plt.legend()
plt.show()