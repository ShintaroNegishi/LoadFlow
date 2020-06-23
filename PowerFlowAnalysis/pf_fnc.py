import numpy as np
import pandas as pd
import cmath

def common_parameter():
    # いずれファイル読み込みからの一般化をする！
    r = np.inf * np.ones((4,4)) # 抵抗
    x = np.zeros((4,4)) # インダクタンス
    b = np.zeros((4,4)) # サセプタンス
    bc = np.zeros(4) # 容量サセプタンス
    P = np.zeros(4)
    Q = np.zeros(4)
    
    n = 4 # ノード数
    
    r[0,1] = 0.01
    x[0,1] = 0.5
    b[0,1] = 0.4
    
    r[1,0] = 0.01
    x[1,0] = 0.5
    b[1,0] = 0.4
    
    r[1,2] = 0.005
    x[1,2] = 0.25
    b[1,2] = 0.2
    
    r[2,1] = 0.005
    x[2,1] = 0.25
    b[2,1] = 0.2
    
    r[2,3] = 0.01
    x[2,3] = 0.5
    b[2,3] = 0.4
    
    r[3,2] = 0.01
    x[3,2] = 0.5
    b[3,2] = 0.4
    
    bc[0] = 0
    bc[1] = 0.1
    bc[2] = 0.1
    bc[3] = 0
    
    # ノードに流入する電力
    P[1] = -0.6
    Q[1] = -0.3
    P[2] = -0.6
    Q[2] = -0.3
    P[3] = 0.6
    
    return r, x, b, bc, P, Q, n

def node_calc(V,theta,r,x,b,n):
    I_dash = np.zeros((n,n), dtype=np.complex) # ノードiからノードjに向かい流出する電流
    Power = np.zeros((n,n), dtype=np.complex) # ノードiからノードjに向かい流出する電力潮流
    V_dot = np.zeros(n, dtype=np.complex)
    
    #V_dot = V * np.exp(1.0j * theta) # ノード電圧（複素表示）
    for ii in range(n):
        V_dot[ii] = V[ii] * cmath.rect(1.0, theta[ii])
    
    # ブランチ潮流の計算
    for ii in range(n):
        for jj in range(n):
            I_dash[ii,jj] = -1.0j * b[ii,jj] / 2.0 * V_dot[ii] + (V_dot[ii]-V_dot[jj]) / (r[ii,jj] + 1.0j * x[ii,jj])
            
    for ii in range(n):
        for jj in range(n):
            Power[ii,jj] = V_dot[ii] * I_dash[ii,jj].conjugate()
            
    return I_dash, Power
    
    
    