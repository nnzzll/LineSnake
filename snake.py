import numpy as np
from scipy.ndimage import gaussian_filter


def calculate_Eint1(contour):
    y = contour[:, 1]
    # d为轮廓线上的控制点间的平均距离
    dist = np.abs(np.diff(y)).tolist()
    dist.append(abs(y[-1]-y[0]))
    d = np.mean(dist)

    # 计算内部能量Eint的一阶项Eint1
    Eint1 = np.zeros(len(y), dtype=np.float32)
    for i in range(len(y)):
        dist = []
        for j in range(-1, 2, 1):
            if(i == len(y)-1):
                dist.append(np.power(d-abs(y[i]+j-y[0]), 2))
            else:
                dist.append(np.power(d-abs(y[i]+j-y[i+1]), 2))
        max_dist = max(dist)
        if(i == len(y)-1):
            Eint1[i] = np.power(d-abs(y[i]-y[0]), 2)/max_dist
        else:
            Eint1[i] = np.power(d-abs(y[i]-y[i+1]), 2)/max_dist
    return Eint1


def calculate_Eint2(contour):
    y = contour[:, 1]

    # 计算内部能量Eint的二阶项Eint2
    Eint2 = np.zeros(len(y), dtype=np.float32)
    for i in range(len(y)):
        if i == 0:
            _max = 0
            for j in range(-1, 2, 1):
                temp = np.power(y[-1]-2*(y[i]+j)+y[i+1], 2)
                _max = temp if temp > _max else _max
            Eint2[i] = np.power(y[-1]-2*y[i]+y[i+1], 2)/_max
        elif i == len(y)-1:
            _max = 0
            for j in range(-1, 2, 1):
                temp = np.power(y[i-1]-2*(y[i]+j)+y[0], 2)
                _max = temp if temp > _max else _max
            Eint2[i] = np.power(y[i-1]-2*y[i]+y[0], 2)/_max
        else:
            _max = 0
            for j in range(-1, 2, 1):
                temp = np.power(y[i-1]-2*(y[i]+j)+y[i+1], 2)
                _max = temp if temp > _max else _max
            Eint2[i] = np.power(y[i-1]-2*y[i]+y[i+1], 2)/_max
    return Eint2


def calculate_Eimg(I, contour):
    x = contour[:, 0]
    y = contour[:, 1]
    length = len(contour)
    # 计算边缘灰阶梯度G,由于为极坐标系,灰阶梯度只关注y方向
    G = np.zeros(length,dtype=int)
    for i in range(length):
        G[i] = np.int32(I[y[i], x[i]]-I[y[i]+2, x[i]])

    # 计算背景亮度,用高斯模糊代替
    Ibg = gaussian_filter(I, 1)

    # 计算图像的边缘对比度特征量C
    C = np.zeros(length)
    for i in range(length):
        C[i] = G[i]/(I[y[i], x[i]]+1e-8)
    # 计算图像力Eimage
    _minC = C.min()
    _maxC = C.max()
    Eimg = (_minC-C)/(_maxC-_minC)
    return Eimg


def calculate_Econ(contour):
    y = contour[:, 1]
    length = len(y)
    # 计算各控制点到重心距离的均值D
    D = np.mean(np.abs(y-y.mean()))

    # 计算约束力
    Econ = np.zeros(length)
    for i in range(length):
        _max = 0
        for j in range(-1, 2, 1):
            temp = np.power(D-abs(y[i]+j-y.mean()), 2)
            _max = temp if temp > _max else _max
        Econ[i] = - np.power(D-abs(y[i]-y.mean()), 2)/_max
    return Econ


def calculate_dist(contour):
    y = contour[:, 1]
    dist = np.zeros(len(y))
    for i in range(len(y)):
        if i == len(y)-1:
            dist[i] = (y[0]-y[i])**2
        else:
            dist[i] = (y[i+1]-y[i])**2
    return dist


def calculate_E(I, contour,alpha,beta,gamma,sigma):
    Eint1 = calculate_Eint1(contour)
    Eint2 = calculate_Eint2(contour)
    Eimg = calculate_Eimg(I, contour)
    Econ = calculate_Econ(contour)
    E_dist = calculate_dist(contour)
    return np.mean(alpha*Eint1+beta*Eint2+0.3*E_dist+gamma*Eimg+sigma*Econ)


def iteratation(I, contour,alpha,beta,gamma,sigma):
    length = len(contour)
    new_contour = contour.copy()
    init_E = calculate_E(I,contour,alpha,beta,gamma,sigma)
    # print("Init E:", init_E)
    E_old = init_E
    E_new = E_old+10
    it = 0
    while(abs(E_new-E_old)>0.01):
        it +=1
        E_old = E_new
        for i in range(length):
            minE = E_old if (i==0) else minE
            minEj = 0
            x = new_contour[:, 0]
            y = new_contour[:, 1]
            for j in range(-5, 6, 1):
                y_new = new_contour[:,1].copy()
                y_new[i] += j
                j_contour = np.array([x, y_new], dtype=int).T
                E_new = calculate_E(I, j_contour,alpha,beta,gamma,sigma)
                if E_new < minE:
                    minE = E_new
                    minEj = j
            y[i] += minEj
            y[i] = 23 if y[i]<23 else y[i]
            y[i] = 124 if y[i]>124 else y[i]
            new_contour = np.array([x, y], dtype=int).T
            # print("It {}\tNo.{}\tminE:{}\tminEj:{}".format(it,i+1, minE, minEj))
        E_new = minE
        # print("E_old:",E_old,"E_new:",E_new)
    return new_contour