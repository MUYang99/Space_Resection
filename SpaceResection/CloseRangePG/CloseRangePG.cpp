// CloseRangePG.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include<stdio.h>
#include"string.h"
#include<stdlib.h>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include <sstream> 
#include <string>
#include<math.h>
#include <opencv2\opencv.hpp>
#include <opencv2\imgproc\types_c.h>

using namespace std;
#define pixsize 0.0051966

//控制点像方坐标和物方坐标
struct PointData {
	double x;
	double y;
	double X;
	double Y;
	double Z;
};

//矩阵转置函数
void MatrixTranposition(double *MatrixOrigin, double *MatrixFinal, int m, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			MatrixFinal[i*m + j] = MatrixOrigin[j*n + i];
		}
	}
}

//矩阵相乘
void MatrixMultiplication(double *MatrixOrign, double *MatrixOrignCopy, double *MatrixPlus, int m, int n, int s) {
	int i, j, k;
	for (i = 0; i < m; i++)
		for (j = 0; j < s; j++)
		{
			MatrixPlus[i*s + j] = 0.0;
			for (k = 0; k < n; k++)
				MatrixPlus[i*s + j] += MatrixOrign[i*n + k] * MatrixOrignCopy[j + k * s];
		}
}

//矩阵求逆函数
void MatrixInversion(double *MatrixInverse, int n) {
	int i, j, k;
	for (k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			if (i != k)
				MatrixInverse[i*n + k] = -MatrixInverse[i*n + k] / MatrixInverse[k*n + k];
		}
		MatrixInverse[k*n + k] = 1 / MatrixInverse[k*n + k];
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
						MatrixInverse[i*n + j] += MatrixInverse[k*n + j] * MatrixInverse[i*n + k];
				}
			}
		}
		for (j = 0; j < n; j++)
		{
			if (j != k)
				MatrixInverse[k*n + j] *= MatrixInverse[k*n + k];
		}
	}
}


int main() {
	/******************************************设定各种初值********************************************/
	//未知数初值设定（11个）
	double x0 = 0.0;
	double y0 = 0.0;
	double f = 36.0;//内方位元素初值设定

	double Xs = 4210, Ys = -160, Zs = -1120;
	double phi = 0.0;
	double w = 0.0;
	double k = 0.0;//外方位元素初值设定

	double k1 = 0.0;
	double k2 = 0.0;//畸变差初值的设定

	//矩阵大小设定
	double *A = new double[2 * 59 * 11];//未知数改正数向量的系数矩阵

	double *L = new double[2 * 59];	//L矩阵为像点像平面坐标的观测值与其像点坐标计算值之差

	double *AT = new double[11 * 2 * 59];
	double *ATL = new double[11 * 1];
	double *ATA = new double[11 * 11];
	double *AX = new double[2 * 59];
	double *v0 = new double[2 * 59];
	double *v1= new double[2 * 59];// 过渡矩阵

	//改正数初值设定
	double R[11] = { 0 };//未知数改正数
	
	double a1 = 0.0, a2 = 0.0, a3 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;//旋转矩阵参数改正数
	double dx = 0.0;
	double dy = 0.0;//畸变差改正数

	//近似值初值设定
	double *xo = new double[59];
	double *yo = new double[59];//像点坐标近似值的定义

	//其他参数初值设定
	double *X_ = new double[59];
	double *Y_ = new double[59];
	double *Z_ = new double[59];//存放共线条件方程中的分子和分母式

	double vv[1];
	vv[0] = 0.0;//单位权中误差的分子

	double r2 = 0.0;//像径的平方
	int count = 1;//迭代次数

	/*************************读取IMG_2324控制点像方坐标文件和控制点物方坐标文件****************************/
	ifstream ifsImageSpace("IMG_2423控制点像方坐标refine.txt");
	if (!ifsImageSpace) {
		printf("Open Image Space Data Error!");
		system("pause");
	}
	vector<double>ImageSpaceData;
	double i;
	while (ifsImageSpace >> i) {
		ImageSpaceData.push_back(i);
	}
	ifsImageSpace.close();//读取IMG_2423控制点像方坐标文件

	ifstream ifsObjectSpace("控制点坐标.txt");
	if (!ifsObjectSpace) {
		printf("Open Object Space Data Error!");
		system("pause");
	}
	vector<double>ObjectSpaceData;
	double o;
	while (ifsObjectSpace >> o) {
		ObjectSpaceData.push_back(o);
	}
	ifsObjectSpace.close();//读取控制点物方坐标文件
	int ObjectNumber = ObjectSpaceData[0];

	/*********************************************坐标转换**************************************************/
	int ergodic = 0;
	PointData points[59];//共有59个控制点
	for (int i = 1; i < ObjectNumber * 4-1; i+=4) {
		for (int j = 0; j < 59 * 3; j+=3) {
			if (ObjectSpaceData[i] == ImageSpaceData[j]) {
				points[ergodic].x = ImageSpaceData[j + 1] * pixsize;
				points[ergodic].y = ImageSpaceData[j + 2] * pixsize;
				points[ergodic].X = ObjectSpaceData[i + 2];
				points[ergodic].Y = ObjectSpaceData[i + 3];
				points[ergodic].Z = -ObjectSpaceData[i + 1];
				ergodic++;
			}
		}
	}

	/*********************************************迭代计算*************************************************/
	do {
		a1 = cos(phi)*cos(k) - sin(phi)*sin(w)*sin(k);
		a2 = -cos(phi)*sin(k) - sin(phi)*sin(w)*cos(k);
		a3 = -sin(phi)*cos(w);
		b1 = cos(w)*sin(k);
		b2 = cos(w)*cos(k);
		b3 = -sin(w);
		c1 = sin(phi)*cos(k) + cos(phi)*sin(w)*sin(k);
		c2 = -sin(phi)*sin(k) + cos(phi)*sin(w)*cos(k);
		c3 = cos(phi)*cos(w);//旋转矩阵R

		for (int i = 0; i < ergodic; i++) {
			X_[i] = a1 * (points[i].X - Xs) + b1 * (points[i].Y - Ys) + c1 * (points[i].Z - Zs);
			Y_[i] = a2 * (points[i].X - Xs) + b2 * (points[i].Y - Ys) + c2 * (points[i].Z - Zs);
			Z_[i] = a3 * (points[i].X - Xs) + b3 * (points[i].Y - Ys) + c3 * (points[i].Z - Zs);
		}//对每个控制点列共线条件方程，存放分子式和分母式

		for (int i = 0; i < ergodic; i++) {
			A[11 * (i * 2 + 0) + 0] = (a1 * f + a3 * (points[i].x - x0)) / Z_[i];//a11i
			A[11 * (i * 2 + 0) + 1] = (b1 * f + b3 * (points[i].x - x0)) / Z_[i];//a12i
			A[11 * (i * 2 + 0) + 2] = (c1 * f + c3 * (points[i].x - x0)) / Z_[i];//a13i
			A[11 * (i * 2 + 0) + 3] = (points[i].y - y0)*sin(w) - ((points[i].x - x0)*((points[i].x - x0)*cos(k) - (points[i].y - y0)*sin(k)) / f + f * cos(k))*cos(w);//a14i
			A[11 * (i * 2 + 0) + 4] = -f * sin(k) - (points[i].x - x0)*((points[i].x - x0)*sin(k) + (points[i].y - y0)*cos(k)) / f;//a15i
			A[11 * (i * 2 + 0) + 5] = (points[i].y - y0);//a16i
			A[11 * (i * 2 + 0) + 6] = (points[i].x - x0) / f;//a17i
			A[11 * (i * 2 + 0) + 7] = 1.0;//a18i
			A[11 * (i * 2 + 0) + 8] = 0.0;//a19i

			r2 = (points[i].x - x0)*(points[i].x - x0) + (points[i].y - y0)*(points[i].y - y0);//向径的平方
			A[11 * (i * 2 + 0) + 9] = (points[i].x - x0)*(r2);
			A[11 * (i * 2 + 0) + 10] = (points[i].x - x0)*(r2*r2);

			A[11 * (i * 2 + 1) + 0] = (a2 * f + a3 * (points[i].y - y0)) / Z_[i];//a21i
			A[11 * (i * 2 + 1) + 1] = (b2 * f + b3 * (points[i].y - y0)) / Z_[i];//a22i
			A[11 * (i * 2 + 1) + 2] = (c2 * f + c3 * (points[i].y - y0)) / Z_[i];//a23i
			A[11 * (i * 2 + 1) + 3] = -(points[i].x - x0)*sin(w) - ((points[i].y - y0)*((points[i].x - x0)*cos(k) - (points[i].y - y0)*sin(k)) / f - f * sin(k))*cos(w);//a24i
			A[11 * (i * 2 + 1) + 4] = -f * cos(k) - (points[i].y - y0)*((points[i].x - x0)*sin(k) + (points[i].y - y0)*cos(k)) / f;//a25i
			A[11 * (i * 2 + 1) + 5] = -(points[i].x - x0);//a26i
			A[11 * (i * 2 + 1) + 6] = (points[i].y - y0) / f;//a27i
			A[11 * (i * 2 + 1) + 7] = 0.0;//a28i
			A[11 * (i * 2 + 1) + 8] = 1.0;//a29i

			A[11 * (i * 2 + 1) + 9] = (points[i].y - y0)*(r2);
			A[11 * (i * 2 + 1) + 10] = (points[i].y - y0)*(r2*r2);

			dx = (points[i].x - x0)*(k1*r2 + k2 * r2*r2);
			dy = (points[i].y - y0)*(k1*r2 + k2 * r2*r2);//畸变差模型

			xo[i] = -f * X_[i] / Z_[i] + x0 - dx;
			yo[i] = -f * Y_[i] / Z_[i] + y0 - dy;//像点坐标近似值

			L[2 * i + 0] = points[i].x - xo[i];
			L[2 * i + 1] = points[i].y - yo[i];//L矩阵的计算
		}
		
		//计算(AT*A)逆*(AT*L)求得物方空间坐标R[3]
		MatrixTranposition(A, AT, 2*59, 11);
		MatrixMultiplication(AT, A, ATA, 11, 2 * 59, 11);
		MatrixInversion(ATA, 11);
		MatrixMultiplication(AT, L, ATL, 11, 2 * 59, 1);
		MatrixMultiplication(ATA, ATL, R, 11, 11, 1);
		MatrixMultiplication(A, R, AX, 2 * 59, 11, 1);
		Xs += R[0], Ys += R[1], Zs += R[2];
		phi += R[3], w += R[4], k += R[5];
		f += R[6], x0 += R[7], y0 += R[8];
		k1 += R[9], k2 += R[10];//改正数
		count++;
	} while ((fabs(R[3]) > 0.000001) || (fabs(R[4]) > 0.000001) || (fabs(R[5]) > 0.000001));

	/**************************************迭代完成后各种中误差计算***********************************************/
	MatrixMultiplication(A, R, AX, 2 * 59, 11, 1);
	for (int i = 0; i < 2 * 59; i++)
	{
		v0[i] = AX[i] - L[i];
		v1[i] = v0[i] / pixsize;
	}//像点坐标残差

	MatrixMultiplication(v0, v0, vv, 1, 2 * 59, 1);
	double sigma = sqrt(vv[0] / (2 * 59 - 11));//单位权中误差
	
	double SigmaOfUnknown[11];
	for (int i = 0; i < 11; i++)
	{
		SigmaOfUnknown[i] = sigma * sqrt(ATA[i * 11 + i]);
	}//十一个未知数的中误差

	/***********************************************输出结果****************************************************/
	ofstream outfile("IMG_2423后交结果.txt", ios::out);
	if (!outfile)
	{
		cerr << "open error" << endl;
		exit(1);
	}
	outfile << "IMG_2423空间后方交会结果：" << endl << endl;
	outfile << "Xs（毫米）  Ys（毫米）  Zs（毫米）  phi（弧度）  w（弧度）  k（弧度）" << endl << Xs << "  " << Ys << "  " << Zs << "  " << phi << "  " << w << "  " << k << endl << endl;
	cout << "IMG_2423空间后方交会结果：" << endl << endl;
	cout << "Xs（毫米）  Ys（毫米）  Zs（毫米）  phi（弧度）  w（弧度）  k（弧度）" << endl << Xs << "  " << Ys << "  " << Zs << "  " << phi << "  " << w << "  " << k << endl << endl;

	outfile << "内方位元素：f,x0,y0 (毫米)" << endl << f << "  " << x0 << "  " << y0 << endl << endl;
	outfile << "k1,k2:   " << k1 << "   " << k2 << endl << endl;
	outfile << "迭代次数:  " << count << endl << endl;
	outfile << "单位权中误差（毫米） " << sigma << endl;
	double sigmanew=sigma/pixsize;
	outfile << "单位权中误差（像素） " << sigmanew << endl << endl;
	cout << "内方位元素：f,x0,y0 (毫米)" << endl << f << "  " << x0 << "  " << y0 << endl << endl;
	cout << "k1,k2:   " << k1 << "   " << k2 << endl << endl;
	cout << "迭代次数:  " << count << endl << endl;
	cout << "单位权中误差（毫米） " << sigma << endl;
	cout << "单位权中误差（像素） " << sigmanew << endl << endl;

	outfile << "各未知数中误差:  " << endl;
	outfile << "mXs =  " << SigmaOfUnknown[0] << endl;
	outfile << "mYs =  " << SigmaOfUnknown[1] << endl;
	outfile << "mZs =  " << SigmaOfUnknown[2] << endl;
	outfile << "mphi =    " << SigmaOfUnknown[3] << endl;
	outfile << "mw =  " << SigmaOfUnknown[4] << endl;
	outfile << "mk =  " << SigmaOfUnknown[5] << endl;
	outfile << "mf =  " << SigmaOfUnknown[6] << endl;
	outfile << "mx0 =  " << SigmaOfUnknown[7] << endl;
	outfile << "my0 =  " << SigmaOfUnknown[8] << endl;
	outfile << "mk1 = " << SigmaOfUnknown[9] << endl;
	outfile << "mk2 = " << SigmaOfUnknown[10] << endl << endl;
	cout << "各未知数中误差:  " << endl;
	cout << "mXs =  " << SigmaOfUnknown[0] << endl;
	cout << "mYs =  " << SigmaOfUnknown[1] << endl;
	cout << "mZs =  " << SigmaOfUnknown[2] << endl;
	cout << "mphi =    " << SigmaOfUnknown[3] << endl;
	cout << "mw =  " << SigmaOfUnknown[4] << endl;
	cout << "mk =  " << SigmaOfUnknown[5] << endl;
	cout << "mf =  " << SigmaOfUnknown[6] << endl;
	cout << "mx0 =  " << SigmaOfUnknown[7] << endl;
	cout << "my0 =  " << SigmaOfUnknown[8] << endl;
	cout << "mk1 = " << SigmaOfUnknown[9] << endl;
	cout << "mk2 = " << SigmaOfUnknown[10] << endl << endl;

	outfile << "像点坐标观测值残差" << endl;
	outfile << "点号    x残差(毫米)     y残差(毫米)     x残差(像素)     y残差(像素)" << endl;
	cout << "像点坐标观测值残差" << endl;
	cout << "点号    x残差(毫米)     y残差(毫米)     x残差(像素)     y残差(像素)" << endl;
	for (int i = 0; i < 59; i++)
	{
		outfile << ImageSpaceData[i*3] << '\t' << v0[2 * i + 0] << '\t' << v0[2 * i + 1] << '\t' << v1[2 * i + 0] << '\t' << v1[2 * i + 1] << endl;
		cout << ImageSpaceData[i * 3] << '\t' << v0[2 * i + 0] << '\t' << v0[2 * i + 1] << '\t' << v1[2 * i + 0] << '\t' << v1[2 * i + 1] << endl;
	}
	outfile.close();
	return 0;
}
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
