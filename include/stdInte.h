#ifndef INTE_H
#define INTE_H

#define MAXWINDOWSIZE 16

//#define WGT  1.5              /* Gauss型窓関数のボケ半径の係数 */
#define WGSS 1.25             /* Gauss型窓関数を２つの矩形で表すための幅の比 small */
#define WGSL 2.5              /* Gauss型窓関数を２つの矩形で表すための幅の比 large */
#define COEFGS 0.041667       /* Gauss型窓関数を２つの矩形で表すための係数 (2D) */

#define WNNDEsHoGD  0.7          /* NNDEGDを四角で測っているための補正 */
#define WNNDEGD  1.2          /* NNDEGDを四角で測っているための補正 */
#define WGS      1.5          /* Gauss型窓関数を1つの矩形で表すための幅の比 */

#define MU  1.0              /* 緩和係数 of Newton */
#define DIRTHRESH 30.0       /* エッジがどうか判定するしきい値 */
#define ROWINTE (ROW + 2 * MAXWINDOWSIZE + 1)
#define COLINTE (COL + 2 * MAXWINDOWSIZE + 1)
#define DNNL    {0, 1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512};
#define NDNNL   19

#define NI             8       /* For inverse matrix */

/*----------------------------------------------------------------------------*/
void copyNormalGpt(double inGpt[3][3], double outGpt[3][3]);
void multplyMV(double inMat[NI][NI + 1], double v[NI]);
void solveLEq(double inMat[NI][NI + 1]);

void calInte(double g_can[ROW][COL], int g_ang[ROW][COL], int inteAng[ROWINTE][COLINTE][9],
		double inteCanDir[ROWINTE][COLINTE][9], double inteDx2Dir[ROWINTE][COLINTE][9], double inteDy2Dir[ROWINTE][COLINTE][9]);
void calInte64(double g_can[ROW][COL], char sHOG[ROW - 4][COL - 4], int inteAng[ROWINTE][COLINTE][64],
		double inteCanDir[ROWINTE][COLINTE][64], double inteDx2Dir[ROWINTE][COLINTE][64], double inteDy2Dir[ROWINTE][COLINTE][64]);

double winpatInte(int g_ang1[ROW][COL], int inteAng[ROWINTE][COLINTE][9]);
double sHoGpatInte(char sHoG1[ROW - 4][COL - 4], int inteAng[ROWINTE][COLINTE][64]);

void gptcorInte(int g_ang1[ROW][COL], double g_can1[ROW][COL],
		           int g_ang2[ROW][COL], double g_can2[ROW][COL],
							 double gwt[ROW][COL], double inteCanDir[ROWINTE][COLINTE][9],
		           double inteDx2Dir[ROWINTE][COLINTE][9], double inteDy2Dir[ROWINTE][COLINTE][9], double dnn, double gpt[3][3]);

void gptcorsHoGInte(char sHoG1[ROW - 4][COL - 4], double g_can1[ROW][COL],
		                char sHoG2[ROW - 4][COL - 4], double g_can2[ROW][COL], double gwt[ROW][COL], double inteCanDir[ROWINTE][COLINTE][64],
		                double inteDx2Dir[ROWINTE][COLINTE][64], double inteDy2Dir[ROWINTE][COLINTE][64], double dnn, double gpt[3][3]);
#endif
