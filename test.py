import time
import ctypes
import pydicom
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, median_filter
from utils import *

# Load C++ dll
lib = ctypes.cdll.LoadLibrary('Lib/snake.dll')


def LineSnake(polar_img: np.ndarray, y: np.ndarray, alpha: float = 0, beta: float = 0.7, gamma: float = 1.8, sigma: float = 0.05, dist: float = 0.3):
    rows, cols = polar_img.shape
    polar_img_ptr = ctypes.cast(
        polar_img.ctypes.data, ctypes.POINTER(ctypes.c_double))
    y_ptr = ctypes.cast(y.ctypes.data, ctypes.POINTER(ctypes.c_double))
    output_y = np.zeros(len(y)).astype(float)
    output_y_ptr = ctypes.cast(
        output_y.ctypes.data, ctypes.POINTER(ctypes.c_double))
    lib.active_contour(polar_img_ptr, y_ptr, output_y_ptr,ctypes.c_int(len(y)), ctypes.c_int(rows), ctypes.c_int(
        cols), ctypes.c_double(alpha), ctypes.c_double(beta), ctypes.c_double(gamma), ctypes.c_double(sigma))
    return output_y.copy()


def main():
    # dicom = pydicom.read_file(
    #     "D:\\Dataset\\IVUS\\IVUSData\\single\\dicom\\IM000001"
    # )
    dicom = pydicom.read_file(
        "D:\\Dataset\\IVUS+DSA\\example\\data-1\\IVUS\\20201124\\041629\\Run1\\PDNQVOLJ")
    data = dicom.pixel_array

    # fig = plt.figure()
    begin = time.time()

    for i in range( len(data)):
        image = data[i].astype(float)
        polar_img = img2polar(image, (128, 360), (255, 255), 255)
        polar_img[:21,:]=0
        I = median_filter(polar_img, 5)
        for _ in range(2):
            I = median_filter(I, 3)
        gauss = gaussian_filter(I, 1)
        y = 35*np.ones([12], float)
        x = np.array([30*i for i in range(12)])

        y_new_lumen = LineSnake(I, y.copy(), 0, 0.2, 1.8, 0, 0.4)
        fig, ax = plt.subplots(1, 2)
        ax[1].imshow(I, 'gray')
        ax[1].plot(x, y_new_lumen, 'o', c='r', label='lumen')
        ax[1].plot(x, y, 'o')

        r_ori = y*2
        r_new_lumen = y_new_lumen*2
        theta = x*np.pi/180
        cart_x, cart_y = polarToCart(r_ori, theta, (255, 255))
        cart_x_new, cart_y_new = polarToCart(r_new_lumen, theta, (255, 255))

        ax[0].imshow(image, 'gray')
        ax[0].plot(cart_x_new, cart_y_new, '.', c='r', label='lumen')
        ax[0].plot(cart_x, cart_y, '.')
        plt.show()

        # plt.imshow(I,'gray')
        # plt.plot(x,y_new_lumen,'.',c='cyan')
        # plt.axis('off')
        # fig.set_size_inches(I.shape[1]/100,I.shape[0]/100)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        # plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
        #                     hspace=0, wspace=0)
        # plt.margins(0, 0)
        # plt.savefig('result/Polar/Frame_{}.png'.format(i+1))
        # plt.clf()
        # y = y_new_lumen

        if(i+1) % 100 == 0:
            print("已处理{}张图片,用时{:.3f}s.".format(i+1, time.time()-begin))


if __name__ == '__main__':
    main()
