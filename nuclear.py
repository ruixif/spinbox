import numpy as np
import pprint
from matplotlib import pyplot as plt
from utils import construct_wigner, three_dim_multiply
from scipy import signal
from math import cos as Cos
from math import sin as Sin



PI = 3.1415926
SQR2 = np.sqrt(2)
mms_to_Joule = 14.41e3 * 1.6e-19 / 2.998e11  #convert mm/s to joule 
gbetaFE_g = 9.16 / (14.41e1 * 1.6 / 2.998) # unit: mm/s .T^-1, value: 0.1191
gbetaFE_e = -0.0676 #unit: mm/s .T^-1
complex_identity = np.zeros((4, 4), dtype = complex)
complex_identity.real = np.eye(4) 
CG = {(0,1):np.sqrt(3/6), 
      (1,1):-np.sqrt(2/6), 
      (2,1):np.sqrt(1/6),
      (3,1):0.0,
      (0,2):0.0,
      (1,2):np.sqrt(1/6),
      (2,2):-np.sqrt(2/6),
      (3,2):np.sqrt(3/6),
        }


class Nuclearham:
    def __init__(self, IS, dEq, eta, B, theta):
        excited = construct_wigner(1.5)
        ground = construct_wigner(0.5)
        
        #unit mm/s
        # ~ self._excited_ham = IS + dEq * (0.5 * np.dot(excited.S_z, excited.S_z) - 5/8*complex_identity + \
                       # ~ eta / 6 * (np.dot(excited.S_x, excited.S_x) - np.dot(excited.S_y, excited.S_y))) + \
                       # ~ gbetaFE_e * (B[0] * excited.S_x + B[1] * excited.S_y + B[2] * excited.S_z)
       
        V_xx, V_yy, V_zz = construct_V(dEq, eta)
        
        V_tensor = np.array([[V_xx, 0,     0],
                             [0,    V_yy,  0],
                             [0,    0,  V_zz]])
        
        g_tensor = np.eye(3) # only to match the function arguments
                             
        # ~ self._excited_ham = IS + V_xx * np.dot(excited.S_x, excited.S_x) + V_yy * np.dot(excited.S_y, excited.S_y) +\
                            # ~ V_zz * np.dot(excited.S_z, excited.S_z) + \
                            # ~ gbetaFE_e * (B[0] * excited.S_x + B[1] * excited.S_y + B[2] * excited.S_z)
                            
        excited_I_ops = np.array([excited.S_x, excited.S_y, excited.S_z])
        ground_I_ops =  np.array([ground.S_x,  ground.S_y,  ground.S_z])
        
        self._excited_ham = IS + three_dim_multiply(V_tensor, operator1 = excited_I_ops, operator2 = excited_I_ops) + \
                    gbetaFE_e * three_dim_multiply(g_tensor, operator1 = excited_I_ops, vector = B)
        self._ground_ham = + gbetaFE_g * three_dim_multiply(g_tensor, operator1 = ground_I_ops, vector = B)
        self._theta = theta
    
    @property
    def excited_ham(self):
        return self._excited_ham
    
    @property
    def ground_ham(self):
        return self._ground_ham
    
    def compute_state(self):
        excited_w, excited_v = np.linalg.eig(self._excited_ham)
        ground_w, ground_v = np.linalg.eig(self._ground_ham)
        
        
        E_prob = []
        total_prob = 0.0
        
        for g_w, g_v in zip(ground_w, ground_v.T):
            for e_w, e_v in zip(excited_w, excited_v.T):
                each_prob = fermi_golden(e_v, g_v, self._theta)
                E_prob.append(((e_w-g_w).real, each_prob))
                total_prob += each_prob
                
        self.E_prob = E_prob
                

def fermi_golden(e_vec, g_vec_temp, theta):
    trans_p = 0
    g_vec = [0, g_vec_temp[0], g_vec_temp[1], 0]
    
    for m1 in range(4):
        for m2 in range(4):
            for m3 in range(4):
                for m4 in range(4):
                    if (m3 >=1) and (m3<=2) and (m4>=1) and (m4<=2):
                        temp_v =  e_vec[m1]* e_vec[m2] * g_vec[m3] * g_vec[m4]* \
                                CG[(m1,m3)] * CG[(m2,m4)] * angle_function(m = m1-m3, m_prime = m2-m4, theta = theta) * \
                                Sin(theta)
                                #the last Sin(theta) is for integration purpose on a sphere.
                        trans_p += abs(temp_v.real)
                                
    return abs(trans_p.real)

def angle_function(m, m_prime, theta):
    if abs(m_prime) == 2 or abs(m_prime) == 3:
        return 0
    
    if abs(m) == 2 or abs(m) == 3:
        return 0
    
    theta_matrix = np.matrix([[0.5 * (1 + (Cos(theta))**2),  0.25 * SQR2 * (1 + Sin(2*theta)), 0.5 * ((Sin(theta))**2)],
                              [0.25 * SQR2 * (1 + Sin(2*theta)), ((Sin(theta))**2), 0.25 * SQR2 * (1 + Sin(2*theta))],
                              [0.5 * ((Sin(theta))**2), 0.25 * SQR2 * (1 + Sin(2*theta)), 0.5 * (1 + (Cos(theta))**2)]]                             
                             )
    index_x = -m + 1
    index_y = -m_prime + 1
    return theta_matrix[index_x, index_y]
    
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def construct_V(dEq, eta):
    V_zz = 1/3 * dEq / np.sqrt(1+1/3* (eta**2))
    V_xx = 0.5 * (-1+eta) * V_zz
    V_yy = 0.5 * (-1-eta) * V_zz
    
    return V_xx, V_yy, V_zz
    

if __name__ == "__main__":
    
    linewidth = 0.05
    N_pts = 1000
    N_powder = 30
    v_max = 20
    total_spectra = np.zeros((N_pts,))
    velocity_range = np.linspace(-v_max, v_max, N_pts)
    angle_space = np.linspace(0.001, PI+0.001, N_powder)
    
    #angle_space = [PI/3]
    for theta in angle_space:
        nuclearham = Nuclearham(IS = 0.0, dEq = 1, eta = 0, B = [0,0,0], theta=theta)
        nuclearham.compute_state()
        
        #pp = pprint.PrettyPrinter()
        #pp.pprint(nuclearham.E_prob)

        for each_E, each_prob in nuclearham.E_prob:
            total_spectra += each_prob * gaussian(velocity_range, each_E, linewidth) / len(angle_space)
        
        
    plt.xlabel('velocity: mm/s')
    plt.ylabel('Absorption Percentage')
    plt.plot(velocity_range, total_spectra)
        
    
    plt.show()
    
    
    
    
    
        






