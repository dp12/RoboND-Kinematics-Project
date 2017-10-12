from sympy import *

q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta_i
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
#
#   
# Create Modified DH parameters
s = {alpha0:     0, a0:         0, d1:  0.75,
     alpha1: -pi/2, a1:      0.35, d2:     0, q2: q2-pi/2,
     alpha2:     0, a2:      1.25, d3:     0,
     alpha3: -pi/2, a3:    -0.054, d4:   1.5,
     alpha4:  pi/2, a4:         0, d5:     0,
     alpha5: -pi/2, a5:         0, d6:     0,}

def TF_Matrix(alpha, a, d, q):
    TF = Matrix([[ cos(q),             -sin(q),            0,            a],
                 [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                 [ sin(q)*sin(alpha), cos(q)*sin(alpha), cos(alpha), cos(alpha)*d],
                 [ 0,                   0,                   0,           1]])
    return TF

T3_4 = TF_Matrix(alpha3, a3, d4, q4).subs(s)
T4_5 = TF_Matrix(alpha4, a4, d5, q5).subs(s)
T5_6 = TF_Matrix(alpha5, a5, d6, q6).subs(s)

#R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
R3_6 = T3_4[0:3, 0:3] * T4_5[0:3, 0:3] * T5_6[0:3, 0:3]
print(R3_6)
#R3_6 = R0_3.T * ROT_EE
