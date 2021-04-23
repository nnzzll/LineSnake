#pragma once
template<typename T = double>
class MyImage {
public:
	T* data;
	long long allocatedSize;
	int cols, rows;

	//����������
	MyImage();
	MyImage(int rows, int cols);
	MyImage(const MyImage& I);
	~MyImage();

	//��ȡ����ֵ�ĺ���
	T getPixel(int i, int j);
};

template<typename T>
MyImage<T>::MyImage() {
	data = nullptr;
	allocatedSize = 0;
	cols = 0;
	rows = 0;
}

template<typename T>
MyImage<T>::MyImage(int height, int width) {
	rows = height;
	cols = width;
	allocatedSize = rows * cols;
	data = nullptr;
}

template<typename T>
MyImage<T>::MyImage(const MyImage& I) {
	cols = I.cols;
	rows = I.rows;
	allocatedSize = I.allocatedSize;
	for (int i = 0; i < allocatedSize; i++)
		data[i] = I.data[i];
}

template<typename T>
MyImage<T>::~MyImage()
{
	if (data)
		delete[] data;
}

template<typename T>
T MyImage<T>::getPixel(int i, int j)
{
	T pixel = data[i * cols + j];
	return pixel;
}