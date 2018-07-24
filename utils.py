import numpy as np
import pprint
from matplotlib import pyplot as plt

#define constant
beta = 9.274009994 * 5.03445 * 1e-2 #wavenumber.T^-1 SI

class construct_wigner(object):
    def __init__(self, S):
        state_num = int(2*S+1)
        S_plus = np.zeros((state_num, state_num), dtype = float)
        S_minus = np.zeros((state_num, state_num), dtype = float)
        S_z = np.zeros((state_num, state_num), dtype = complex)
        
        for i in range(state_num-1):
            m = S-i-1
            S_plus[i][i+1] = np.sqrt(S*(S+1) - m*(m+1))
        
        for i in range(state_num-1):
            m = S-i
            S_minus[i+1][i] = np.sqrt(S*(S+1) -m*(m-1))
    
        for i in range(state_num):
            m = S-i
            S_z.real[i][i] = m
        
        S_x = np.zeros((state_num, state_num), dtype = complex)
        S_y = np.zeros((state_num, state_num), dtype = complex)
        
        
        S_x.real = 0.5 *(S_plus + S_minus) #Sx = [(S+) + (S-)]/2
        S_y.imag = (-0.5) *(S_plus - S_minus) #Sy = -[(S+) + (S-)]/2  j
        
        self._S_plus = S_plus
        self._S_minus = S_minus
        self._S_x = S_x
        self._S_y = S_y
        self._S_z = S_z
        
    @property
    def S_x(self):
        return self._S_x
    
    @property
    def S_y(self):
        return self._S_y
    
    @property
    def S_z(self):
        return self._S_z

def construct_D(D, E_D):
    E = E_D * D
    D_xx = -1/3 * D + E
    D_yy = -1/3 * D - E
    D_zz =  2/3 * D
    return D_xx, D_yy, D_zz #unit in wavenum

class electronham:
    def __init__(self, D, E_D, B, spin, g):
        D_xx, D_yy, D_zz = construct_D(D, E_D)
        D_tensor = np.array([[D_xx, 0,     0],
                             [0,    D_yy,  0],
                             [0,    0,  D_zz]])
        g_tensor = np.array([[g[0], 0,     0],
                             [0,    g[1],  0],
                             [0,    0,  g[2]]])
        
        S_x, S_y, S_z = spin.S_x, spin.S_y, spin.S_z
        S_operator = np.array([S_x, S_y, S_z])
        
        spinham = three_dim_multiply(D_tensor, operator1 = S_operator, operator2 = S_operator) + \
                    beta * three_dim_multiply(g_tensor, operator1 = S_operator, vector = B)
        
        # ~ spinham = D_xx * np.dot(S_x, S_x) + D_yy * np.dot(S_y, S_y) + D_zz *np.dot(S_z, S_z) \
                # ~ + beta *(g[0] * B[0] * S_x + g[1] * B[1] * S_y + g[2] * B[2] * S_z)
        self.spinham = spinham
        self.S_x = S_x
        self.S_y = S_y
        self.S_z = S_z

def compute_expect(operator, wav):
    return np.diagonal(np.dot(wav.T, np.dot(operator, wav))).real

def three_dim_multiply(tensor, **ops):
    total_value = np.zeros_like(ops['operator1'][0])
    if 'operator2' in ops:
        for i in range(3):
            for j in range(3):
                total_value += tensor[i,j] * np.dot(ops['operator1'][i], ops['operator2'][j])
        return total_value
    elif 'vector' in ops:
        for i in range(3):
            for j in range(3):
                total_value += tensor[i,j] * ops['operator1'][i] * ops['vector'][j]
        return total_value
    else:
        raise RuntimeError
        
def eulerAnglesToRotationMatrix(theta) :
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
    R = np.dot(R_z, np.dot( R_y, R_x ))
    # three angles are in the list of theta and the rotation operations are z-y-z
    return R



if __name__ == "__main__":
    exit()
    # ~ Bz_acc = np.linspace(0.0, 20.0, 200)
    # ~ eigen_w_acc = []
    # ~ spindef = construct_wigner(2.5)
    # ~ for Bz in Bz_acc:
        # ~ test_ham = electronham(20, 0, [0, 0, Bz], spindef, [2,2,2])
        # ~ eigen_w, eigen_v = np.linalg.eig(test_ham.spinham)
        # ~ eigen_w_acc.append(eigen_w)
        
    # ~ energy = np.array(eigen_w_acc) 
    # ~ plt.plot(Bz_acc, energy.real)
    # ~ plt.show()

    

    
    
