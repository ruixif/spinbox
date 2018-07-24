import numpy as np
import pprint
from matplotlib import pyplot as plt
from utils import electronham, construct_wigner, compute_expect
from operator import mul

if __name__ == "__main__":
    spindef = construct_wigner(2.5)
    test_ham = electronham(20, 0.2, [0, 0, 0.1], spindef, [2,2,2])
    eigen_w, eigen_v = np.linalg.eig(test_ham.spinham)

    
    S_x_expect, S_y_expect, S_z_expect = compute_expect(test_ham.S_x, eigen_v), \
                                         compute_expect(test_ham.S_y, eigen_v), \
                                         compute_expect(test_ham.S_z, eigen_v)
    
    A = [-22, -22, -22]
    
    B_int = []
    
    for i in range(len(S_x_expect)):
        B_int.append([S_x_expect[i] * A[0], S_y_expect[i] * A[1], S_z_expect[i] * A[2]])
    
    print(B_int)

