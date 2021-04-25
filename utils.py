import cv2
import numpy as np
import matplotlib.pyplot as plt

def polarToCart(r:list,theta:list,center:tuple):
    x = [center[0]+r[i]*np.cos(theta[i]) for i in range(len(r))]
    y = [center[1]-r[i]*np.sin(theta[i]) for i in range(len(r))]
    return x,y

def img2polar(image: np.ndarray, size: tuple, center: tuple, radius: int) -> np.ndarray:
    '''直角坐标系下的图像转换到极坐标系下,方向为顺时针方向,x轴正方向为0°
    '''
    flags = cv2.INTER_CUBIC + cv2.WARP_FILL_OUTLIERS + cv2.WARP_POLAR_LINEAR
    linear_polar_image = cv2.warpPolar(image, size, center, radius, flags)
    dst = np.rot90(linear_polar_image, -1)
    return dst

def plot_result(polar_img:np.ndarray,x,y,x_ori,y_ori):
    plt.imshow(polar_img,cmap='gray')
    plt.plot(x,y,'o',c='cyan')
    plt.plot(x_ori,y_ori,'o')
    plt.show()