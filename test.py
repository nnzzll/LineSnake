import cv2
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, median_filter
from utils import *

# Load C++ dll
lib = ctypes.cdll.LoadLibrary('Lib/snake.dll')


def LineSnake(polar_img: np.ndarray, gaussian_img: np.ndarray, y: np.ndarray, alpha: float = 0, beta: float = 0.7, gamma: float = 1.8, sigma: float = 0.05):
    rows, cols = polar_img.shape
    polar_img_ptr = ctypes.cast(
        polar_img.ctypes.data, ctypes.POINTER(ctypes.c_double))
    gauss_ptr = ctypes.cast(gaussian_img.ctypes.data,
                            ctypes.POINTER(ctypes.c_double))
    y_ptr = ctypes.cast(y.ctypes.data, ctypes.POINTER(ctypes.c_double))
    lib.active_contour(polar_img_ptr, gauss_ptr, y_ptr, ctypes.c_int(rows), ctypes.c_int(
        cols), ctypes.c_double(alpha), ctypes.c_double(beta), ctypes.c_double(gamma), ctypes.c_double(sigma))
    return y.copy()

def main():
    image = cv2.imread(
        'C:\\Users\\Administrator\\Desktop\\ivus\\image\\1.png', 0).astype(float)
    polar_img = img2polar(image, (128, 360), (255, 255), 255)
    I = median_filter(polar_img, 5)
    I = median_filter(I, 3)
    gauss = gaussian_filter(I, 1)

    x = np.array([10*i for i in range(36)])
    y = 35*np.ones([36], float)
    y_new = LineSnake(I,gauss,y.copy(),0,1,1.8,0.05)
    plot_result(I,x,y_new,x,y)

    r_ori = y*2
    r_new = y_new*2
    theta = x*np.pi/180
    cart_x,cart_y = polarToCart(r_ori,theta,(255,255))
    cart_x_new,cart_y_new = polarToCart(r_new,theta,(255,255))
    plot_result(image,cart_x_new,cart_y_new,cart_x,cart_y)

if __name__ == '__main__':
    main()
