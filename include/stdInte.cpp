#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "parameter.h"
#include "stdInte.h"
#include "utility.h"
#include "stdGpt.h"
//-----------------------------------------------------------------------------
void copyNormalGpt(double inGpt[3][3], double outGpt[3][3]) {
	int i, j;
	double nf = 1.0 / inGpt[2][2];
	for(i = 0 ; i < 3 ; ++i) {
		for(j = 0 ; j < 3 ; ++j) {
			outGpt[i][j] = nf * inGpt[i][j];
		}
	}
}

void multplyMV(double inMat[NI][NI + 1], double v[NI]) {
	int i, j;
	double sum;
	for (i = 0 ; i < NI ; ++i) {
		sum = 0.0;
		for (j = 0 ; j < NI ; ++j) {
			sum += inMat[i][j] * inMat[j][NI];
		}
		v[i] = sum;
	}
}

void solveLEq(double inMat[NI][NI + 1]) {
	int i, j, j2, maxI;
	double tmp, pivVal;

	//printNxN1(inMat);
	for (j = 0 ; j < NI ; ++j) {
		/* Search pivot */
		maxI = j;
		for (i = j + 1 ; i < NI ; ++i) {
			if (fabs(inMat[i][j]) > fabs(inMat[maxI][j])) maxI = i;
		}
		pivVal = inMat[maxI][j];
		if (maxI != j) {
			/* Swap j th row and maxI th row */
			for (j2 = 0 ; j2 < NI + 1 ; ++j2) {
				tmp = inMat[j][j2];  inMat[j][j2]  = inMat[maxI][j2] / pivVal;  inMat[maxI][j2]  = tmp;
			}
		} else {
			for (j2 = 0 ; j2 < NI + 1; ++j2) inMat[j][j2] /= pivVal;
		}
		for (i = 0 ; i < NI ; ++i) {
			if (i == j) continue;
			pivVal = inMat[i][j];
			for (j2 = 0 ; j2 < NI + 1 ; ++j2) inMat[i][j2]  -= pivVal * inMat[j][j2];
		}
		//printNxN(inMat);
	}
}

//-----------------------------------------------------------------------------

void calInte64(double g_can[ROW][COL], char sHOG[ROW - 4][COL - 4], int inteAng[ROWINTE][COLINTE][64],
		double inteCanDir[ROWINTE][COLINTE][64], double inteDx2Dir[ROWINTE][COLINTE][64], double inteDy2Dir[ROWINTE][COLINTE][64]) {
	int x, y, d, xInte, yInte, angInte[64];
	int maxWinP = MAXWINDOWSIZE + 1;
	double canDirInte[64], dx2DirInte[64], dy2DirInte[64];

	/* Set init data */
	for (y = 0 ; y < ROWINTE; ++y) {
		for (x = 0 ; x < COLINTE ; ++x) {
			for (d = 0 ; d < 64 ; ++d) {
				inteAng[y][x][d] = 0;
				inteCanDir[y][x][d] = 0.0;
				inteDx2Dir[y][x][d] = 0.0;
				inteDy2Dir[y][x][d] = 0.0;
			}
		}
	}

	for (y = 2 ; y < ROW - 2 ; ++y) {
		for (x = 2 ; x < COL - 2 ; ++x) {
      if (sHOG[y - 2][x - 2] == -1)
        continue;
      else
        d = sHoG2Idx(sHOG[y - 2][x - 2]);
			xInte = x + maxWinP;
			yInte = y + maxWinP;
			inteAng[yInte][xInte][d] = 1;
			inteCanDir[yInte][xInte][d] = g_can[y][x];
			inteDx2Dir[yInte][xInte][d] = g_can[y][x] * (x - CX);
			inteDy2Dir[yInte][xInte][d] = g_can[y][x] * (y - CY);
		}
	}

	/* Calculate Integral for x direction */
	for (yInte = maxWinP ; yInte < ROW + maxWinP ; ++yInte) {
		for (d = 0 ; d < 64 ; ++d) {
			angInte[d] = 0;
			canDirInte[d] = 0.0;
			dx2DirInte[d] = 0.0;
			dy2DirInte[d] = 0.0;
		}
		for (xInte = maxWinP ; xInte < COLINTE ; ++xInte) {
			for (d = 0 ; d < 64 ; ++d) {
				angInte[d]    = inteAng[yInte][xInte][d]    = inteAng[yInte][xInte][d]    + angInte[d];
				canDirInte[d] = inteCanDir[yInte][xInte][d] = inteCanDir[yInte][xInte][d] + canDirInte[d];
				dx2DirInte[d] = inteDx2Dir[yInte][xInte][d] = inteDx2Dir[yInte][xInte][d] + dx2DirInte[d];
				dy2DirInte[d] = inteDy2Dir[yInte][xInte][d] = inteDy2Dir[yInte][xInte][d] + dy2DirInte[d];
			}
		}
	}

	/* Calculate Integral for y direction */
	for (xInte = maxWinP ; xInte < COLINTE ; ++xInte) {  /* 後側余裕も積分した後で0ではないので，COLINTE */
		for (d = 0 ; d < 64 ; ++d) {
			angInte[d] = 0;
			canDirInte[d] = 0.0;
			dx2DirInte[d] = 0.0;
			dy2DirInte[d] = 0.0;
		}
		for (yInte = maxWinP ; yInte < ROWINTE ; ++yInte) {
			for (d = 0 ; d < 64 ; ++d) {
				angInte[d]    = inteAng[yInte][xInte][d]    = inteAng[yInte][xInte][d]    + angInte[d];
				canDirInte[d] = inteCanDir[yInte][xInte][d] = inteCanDir[yInte][xInte][d] + canDirInte[d];
				dx2DirInte[d] = inteDx2Dir[yInte][xInte][d] = inteDx2Dir[yInte][xInte][d] + dx2DirInte[d];
				dy2DirInte[d] = inteDy2Dir[yInte][xInte][d] = inteDy2Dir[yInte][xInte][d] + dy2DirInte[d];
			}
		}
	}
}

void calInte(double g_can[ROW][COL], int g_ang[ROW][COL], int inteAng[ROWINTE][COLINTE][9],
		double inteCanDir[ROWINTE][COLINTE][9], double inteDx2Dir[ROWINTE][COLINTE][9], double inteDy2Dir[ROWINTE][COLINTE][9]) {
	int x, y, d, xInte, yInte, angInte[9];
	int maxWinP = MAXWINDOWSIZE + 1;
	double canDirInte[9], dx2DirInte[9], dy2DirInte[9];

	/* Set init data */
	for (y = 0 ; y < ROWINTE; ++y) {
		for (x = 0 ; x < COLINTE ; ++x) {
			for (d = 0 ; d < 9 ; ++d) {
				inteAng[y][x][d] = 0;
				inteCanDir[y][x][d] = 0.0;
				inteDx2Dir[y][x][d] = 0.0;
				inteDy2Dir[y][x][d] = 0.0;
			}
		}
	}
	for (y = 0 ; y < ROW ; ++y) {
		for (x = 0 ; x < COL ; ++x) {
			d = g_ang[y][x];
			xInte = x + maxWinP;
			yInte = y + maxWinP;
			inteAng[yInte][xInte][d] = 1;
			inteCanDir[yInte][xInte][d] = g_can[y][x];
			inteDx2Dir[yInte][xInte][d] = g_can[y][x] * (x - CX);
			inteDy2Dir[yInte][xInte][d] = g_can[y][x] * (y - CY);
		}
	}

	/* Calculate Integral for x direction */
	for (yInte = maxWinP ; yInte < ROW + maxWinP ; ++yInte) {
		for (d = 0 ; d < 9 ; ++d) {
			angInte[d] = 0;
			canDirInte[d] = 0.0;
			dx2DirInte[d] = 0.0;
			dy2DirInte[d] = 0.0;
		}
		for (xInte = maxWinP ; xInte < COLINTE ; ++xInte) {
			for (d = 0 ; d < 9 ; ++d) {
				angInte[d]    = inteAng[yInte][xInte][d]    = inteAng[yInte][xInte][d]    + angInte[d];
				canDirInte[d] = inteCanDir[yInte][xInte][d] = inteCanDir[yInte][xInte][d] + canDirInte[d];
				dx2DirInte[d] = inteDx2Dir[yInte][xInte][d] = inteDx2Dir[yInte][xInte][d] + dx2DirInte[d];
				dy2DirInte[d] = inteDy2Dir[yInte][xInte][d] = inteDy2Dir[yInte][xInte][d] + dy2DirInte[d];
			}
		}
	}

	/* Calculate Integral for y direction */
	for (xInte = maxWinP ; xInte < COLINTE ; ++xInte) {  /* 後側余裕も積分した後で0ではないので，COLINTE */
		for (d = 0 ; d < 9 ; ++d) {
			angInte[d] = 0;
			canDirInte[d] = 0.0;
			dx2DirInte[d] = 0.0;
			dy2DirInte[d] = 0.0;
		}
		for (yInte = maxWinP ; yInte < ROWINTE ; ++yInte) {
			for (d = 0 ; d < 9 ; ++d) {
				angInte[d]    = inteAng[yInte][xInte][d]    = inteAng[yInte][xInte][d]    + angInte[d];
				canDirInte[d] = inteCanDir[yInte][xInte][d] = inteCanDir[yInte][xInte][d] + canDirInte[d];
				dx2DirInte[d] = inteDx2Dir[yInte][xInte][d] = inteDx2Dir[yInte][xInte][d] + dx2DirInte[d];
				dy2DirInte[d] = inteDy2Dir[yInte][xInte][d] = inteDy2Dir[yInte][xInte][d] + dy2DirInte[d];
			}
		}
	}
}

double winpatInte(int g_ang1[ROW][COL], int inteAng[ROWINTE][COLINTE][9]) {
	int x1, y1, wN, ang1, dnn = 0, count = 0;
	int maxWinP = MAXWINDOWSIZE + 1;
	int pPos, mPos, sectInte, nDnnL = NDNNL;
	double ddnn;
	int dnnL[] = DNNL;

	for (wN = 0 ; wN < NDNNL ; ++wN) {
		if (dnnL[wN] >= MAXWINDOWSIZE) {
			nDnnL = wN + 1;
			break;
		}
	}

	for (y1 = MARGINE ; y1 < ROW - MARGINE ; ++y1) {
		for (x1 = MARGINE ; x1 < COL - MARGINE ; ++x1) {
			if ((ang1 = g_ang1[y1][x1]) != 8) {
				for (wN = 0 ; wN < nDnnL ; ++wN) {
					pPos = maxWinP + dnnL[wN];
					mPos = MAXWINDOWSIZE - dnnL[wN];
					sectInte = inteAng[y1 + pPos][x1 + pPos][ang1] - inteAng[y1 + pPos][x1 + mPos][ang1]
																																												 -inteAng[y1 + mPos][x1 + pPos][ang1] + inteAng[y1 + mPos][x1 + mPos][ang1];
					if (sectInte > 0) {
//						printf("(%d, %d) sectInte = %d dnn = %d \n", x1, y1, sectInte, dnnL[wN]);
						++count;
						dnn += dnnL[wN];
						break;
					}
				}
			}
		}
	}
	//printf("count = %d dnn = %d \n", count, dnn);
	if (count == 0) ddnn = MAXWINDOWSIZE;
	else ddnn = (double) dnn / count;
	return ddnn;
}

double sHoGpatInte(char sHoG1[ROW - 4][COL - 4], int inteAng[ROWINTE][COLINTE][64]) {
	int x1, y1, wN, ang1, dnn = 0, count = 0;
	int maxWinP = MAXWINDOWSIZE + 1;
	int pPos, mPos, sectInte, nDnnL = NDNNL;
	double ddnn;
	int dnnL[] = DNNL;

	for (wN = 0 ; wN < NDNNL ; ++wN) {
		if (dnnL[wN] >= MAXWINDOWSIZE) {
			nDnnL = wN + 1;
			break;
		}
	}

	for (y1 = MARGINE + 2 ; y1 < ROW - MARGINE - 2; ++y1) {
		for (x1 = MARGINE + 2 ; x1 < COL - MARGINE - 2 ; ++x1) {
			if ((ang1 = sHoG2Idx(sHoG1[y1][x1])) != -1) {
				for (wN = 0 ; wN < nDnnL ; ++wN) {
					pPos = maxWinP + dnnL[wN];
					mPos = MAXWINDOWSIZE - dnnL[wN];
					sectInte = inteAng[y1 + pPos][x1 + pPos][ang1] - inteAng[y1 + pPos][x1 + mPos][ang1]
																																												 -inteAng[y1 + mPos][x1 + pPos][ang1] + inteAng[y1 + mPos][x1 + mPos][ang1];
					if (sectInte > 0) {
//						printf("(%d, %d) sectInte = %d dnn = %d \n", x1, y1, sectInte, dnnL[wN]);
						++count;
						dnn += dnnL[wN];
						break;
					}
				}
			}
		}
	}
	//printf("count = %d dnn = %d \n", count, dnn);
	if (count == 0) ddnn = MAXWINDOWSIZE;
	else ddnn = (double) dnn / count;
	return ddnn;
}

void gptcorInte(int g_ang1[ROW][COL], double g_can1[ROW][COL],
		int g_ang2[ROW][COL], double g_can2[ROW][COL], double gwt[ROW][COL], double inteCanDir[ROWINTE][COLINTE][9],
		double inteDx2Dir[ROWINTE][COLINTE][9], double inteDy2Dir[ROWINTE][COLINTE][9], double dnn, double gpt[3][3]) {
	/* determination of optimal GAT components */
	/* that yield the maximal correlation value */
	int x1, y1;
	int i, j, loop;
	double g0, gx1, gy1, gx2, gy2;
	double gx1x1, gx1y1, gy1y1, gx1x2, gx1y2, gy1x2, gy1y2;
	double t0, tx2, ty2;
	double gx1x1x1, gx1x1y1, gx1y1y1, gy1y1y1;
	double gx1x1x2, gx1x1y2, gx1y1x2, gx1y1y2, gy1y1x2, gy1y1y2;
	double gx1x1x1x1, gx1x1x1y1, gx1x1y1y1, gx1y1y1y1, gy1y1y1y1;
	double tx1, ty1, tx1x1, tx1y1, ty1y1, tx1x1x1, tx1x1y1, tx1y1y1, ty1y1y1;
	double U[NI][NI + 1], U0[NI][NI+1];
	double V[NI][NI + 1], v[NI];
	double dx1, dy1;
	double tGpt1[3][3], tGpt2[3][3];
	double detA, r, rg0;
	double Ainv11, Ainv12, Ainv21, Ainv22, grAinv11, grAinv12, grAinv21, grAinv22;
	double var = WGT * WGT * dnn * dnn;

	int ang1;
#ifdef TWOWINDOW
	int windowSSS = (int) (WGSS * WGT * dnn + 0.9999);
	int windowSSL = (int) (WGSL * WGT * dnn + 0.9999);
	int pPosS, mPosS, pPosL, mPosL;
	if (windowSSL > MAXWINDOWSIZE) {
		windowSSS = MAXWINDOWSIZE / 2;
		windowSSL = MAXWINDOWSIZE;
	}
	pPosS = MAXWINDOWSIZE + 1 + windowSSS;
	mPosS = MAXWINDOWSIZE - windowSSS;
	pPosL = MAXWINDOWSIZE + 1 + windowSSL;
	mPosL = MAXWINDOWSIZE - windowSSL;
#else
	int windowS = (int) (WGS * WGT * dnn + 0.9999);
	int pPos, mPos;
	if (windowS > MAXWINDOWSIZE) windowS = MAXWINDOWSIZE;
	pPos = MAXWINDOWSIZE + 1 + windowS;
	mPos = MAXWINDOWSIZE - windowS;
#endif

	/* Gaussian weigthed mean values */
	g0 = gx1 = gy1 = gx2 = gy2 = 0.0;
	gx1x1 = gx1y1 = gy1y1 = 0.0;
	gx1x2 = gx1y2 = gy1x2 = gy1y2 = 0.0;
	gx1x1x1 = gx1x1y1 = gx1y1y1 = gy1y1y1 = 0.0;
	gx1x1x2 = gx1x1y2 = gx1y1x2 = gx1y1y2 = gy1y1x2 = gy1y1y2 = 0.0;
	gx1x1x1x1 = gx1x1x1y1 = gx1x1y1y1 = gx1y1y1y1 = gy1y1y1y1 = 0.0;

	for (y1 = MARGINE ; y1 < ROW - MARGINE ; y1++) {
		dy1 = y1 - CY;
		for (x1 = MARGINE ; x1 < COL - MARGINE ; x1++) {
			dx1 = x1 - CX;

			ang1 = g_ang1[y1][x1];
#ifdef TWOWINDOW
			t0 = tx2 = ty2 = COEFGS;
			t0  *= inteCanDir[y1 + pPosL][x1 + pPosL][ang1] - inteCanDir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteCanDir[y1 + mPosL][x1 + pPosL][ang1] + inteCanDir[y1 + mPosL][x1 + mPosL][ang1];
			tx2 *= inteDx2Dir[y1 + pPosL][x1 + pPosL][ang1] - inteDx2Dir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteDx2Dir[y1 + mPosL][x1 + pPosL][ang1] + inteDx2Dir[y1 + mPosL][x1 + mPosL][ang1];
			ty2 *= inteDy2Dir[y1 + pPosL][x1 + pPosL][ang1] - inteDy2Dir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteDy2Dir[y1 + mPosL][x1 + pPosL][ang1] + inteDy2Dir[y1 + mPosL][x1 + mPosL][ang1];
			t0  += inteCanDir[y1 + pPosS][x1 + pPosS][ang1] - inteCanDir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteCanDir[y1 + mPosS][x1 + pPosS][ang1] + inteCanDir[y1 + mPosS][x1 + mPosS][ang1];
			tx2 += inteDx2Dir[y1 + pPosS][x1 + pPosS][ang1] - inteDx2Dir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteDx2Dir[y1 + mPosS][x1 + pPosS][ang1] + inteDx2Dir[y1 + mPosS][x1 + mPosS][ang1];
			ty2 += inteDy2Dir[y1 + pPosS][x1 + pPosS][ang1] - inteDy2Dir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteDy2Dir[y1 + mPosS][x1 + pPosS][ang1] + inteDy2Dir[y1 + mPosS][x1 + mPosS][ang1];
#else
			t0  = inteCanDir[y1 + pPos][x1 + pPos][ang1] - inteCanDir[y1 + pPos][x1 + mPos][ang1]
					 - inteCanDir[y1 + mPos][x1 + pPos][ang1] + inteCanDir[y1 + mPos][x1 + mPos][ang1];
			tx2 = inteDx2Dir[y1 + pPos][x1 + pPos][ang1] - inteDx2Dir[y1 + pPos][x1 + mPos][ang1]
					 - inteDx2Dir[y1 + mPos][x1 + pPos][ang1] + inteDx2Dir[y1 + mPos][x1 + mPos][ang1];
			ty2 = inteDy2Dir[y1 + pPos][x1 + pPos][ang1] - inteDy2Dir[y1 + pPos][x1 + mPos][ang1]
					 - inteDy2Dir[y1 + mPos][x1 + pPos][ang1] + inteDy2Dir[y1 + mPos][x1 + mPos][ang1];
#endif
			t0 *= g_can1[y1][x1]; tx2 *= g_can1[y1][x1]; ty2 *= g_can1[y1][x1];


			g0    += t0;
			gx1   += tx1 = t0  * dx1;
			gy1   += ty1 = t0  * dy1;
			gx1x1 += tx1x1 = tx1 * dx1;
			gx1y1 += tx1y1 = tx1 * dy1;
			gy1y1 += ty1y1 = ty1 * dy1;
			gx1x1x1 += tx1x1x1 = tx1x1 * dx1;
			gx1x1y1 += tx1x1y1 = tx1x1 * dy1;
			gx1y1y1 += tx1y1y1 = tx1y1 * dy1;
			gy1y1y1 += ty1y1y1 = ty1y1 * dy1;
			gx1x1x1x1 += tx1x1x1 * dx1;
			gx1x1x1y1 += tx1x1x1 * dy1;
			gx1x1y1y1 += tx1x1y1 * dy1;
			gx1y1y1y1 += tx1y1y1 * dy1;
			gy1y1y1y1 += ty1y1y1 * dy1;

			gx2   += tx2;
			gy2   += ty2;
			gx1x2 += tx2 * dx1;
			gx1y2 += ty2 * dx1;
			gy1x2 += tx2 * dy1;
			gy1y2 += ty2 * dy1;
			gx1x1x2 += tx2 * dx1 * dx1;
			gx1x1y2 += ty2 * dx1 * dx1;
			gx1y1x2 += tx2 * dx1 * dy1;
			gx1y1y2 += ty2 * dx1 * dy1;
			gy1y1x2 += tx2 * dy1 * dy1;
			gy1y1y2 += ty2 * dy1 * dy1;
		}
	}
	 // printf("dnn = %f g0 = %f gx1 = %f  gy1 = %f gx2 = %f  gy2 = %f\n", dnn, g0, gx1, gy1, gx2, gy2);
	 // printf("dnn = %f normal gx1x1x1x1 = %f gx1x1x1y1 = %f  gx1x1y1y1 = %f gx1y1y1y1 = %f  gy1y1y1y1y1 = %f\n", dnn, gx1x1x1x1/g0, gx1x1x1y1/g0, gx1x1y1y1/g0, gx1y1y1y1/g0, gy1y1y1y1/g0);

	if (g0 == 0.0) {
#ifdef PRINT
		printf("GAT calculation failure by zero sum!!!\n");
#endif
	}
	/* Newton method */
	r = 0.5 * var;
	rg0 = r * g0;

	V[0][0] = /* D11 */  gx1x1;
	V[1][0] = /* D21 */  gx1y1;
	V[2][0] = /* 0   */  0.0;
	V[3][0] = /* 0   */  0.0;
	V[4][0] = /* m1  */  gx1;
	V[5][0] = /* 0   */  0.0;
	V[6][0] = /* E11 */  gx1x1x1;
	V[7][0] = /* E31 */  gx1x1y1;
	V[0][1] = /* D12 */  gx1y1;
	V[1][1] = /* D22 */  gy1y1;
	V[2][1] = /* 0   */  0.0;
	V[3][1] = /* 0   */  0.0;
	V[4][1] = /* m2  */  gy1;
	V[5][1] = /* 0   */  0.0;
	V[6][1] = /* E12 */  gx1x1y1;
	V[7][1] = /* E32 */  gx1y1y1;
	V[0][2] = /* 0   */  0.0;
	V[1][2] = /* 0   */  0.0;
	V[2][2] = /* D11 */  gx1x1;
	V[3][2] = /* D21 */  gx1y1;
	V[4][2] = /* 0   */  0.0;
	V[5][2] = /* m1  */  gx1;
	V[6][2] = /* E21 */  gx1x1y1;
	V[7][2] = /* E41 */  gx1y1y1;
	V[0][3] = /* 0   */  0.0;
	V[1][3] = /* 0   */  0.0;
	V[2][3] = /* D12 */  gx1y1;
	V[3][3] = /* D22 */  gy1y1;
	V[4][3] = /* 0   */  0.0;
	V[5][3] = /* m2  */  gy1;
	V[6][3] = /* E22 */  gx1y1y1;
	V[7][3] = /* E42 */  gy1y1y1;
	V[0][4] = /* m1  */  gx1;
	V[1][4] = /* m2  */  gy1;
	V[2][4] = /* 0   */  0.0;
	V[3][4] = /* 0   */  0.0;
	V[4][4] = /* 1   */  g0;
	V[5][4] = /* 0   */  0.0;
	V[6][4] = /* D11 */  gx1x1;
	V[7][4] = /* D21 */  gx1y1;
	V[0][5] = /* 0   */  0.0;
	V[1][5] = /* 0   */  0.0;
	V[2][5] = /* m1  */  gx1;
	V[3][5] = /* m2  */  gy1;
	V[4][5] = /* 0   */  0.0;
	V[5][5] = /* 1   */  g0;
	V[6][5] = /* D12 */  gx1y1;
	V[7][5] = /* D22 */  gy1y1;

	V[0][6] = /*-E11 */ -gx1x1x1;
	V[1][6] = /*-E21 */ -gx1x1y1;
	V[2][6] = /*-E31 */ -gx1x1y1;
	V[3][6] = /*-E41 */ -gx1y1y1;
	V[4][6] = /*-D11 */ -gx1x1;
	V[5][6] = /*-D21 */ -gx1y1;
	V[6][6] = /*-F11 - r * D11 */  -gx1x1x1x1 - gx1x1y1y1 - r * gx1x1;
	V[7][6] = /*-F21 - r * D21 */  -gx1x1x1y1 - gx1y1y1y1 - r * gx1y1;

	V[0][7] = /*-E12 */ -gx1x1y1;
	V[1][7] = /*-E22 */ -gx1y1y1;
	V[2][7] = /*-E32 */ -gx1y1y1;
	V[3][7] = /*-E42 */ -gy1y1y1;
	V[4][7] = /*-D12 */ -gx1y1;
	V[5][7] = /*-D22 */ -gy1y1;
	V[6][7] = /*-F12 - r * D12 */ -gx1x1x1y1 - gx1y1y1y1 - r * gx1y1;
	V[7][7] = /*-F22 - r * D22 */ -gx1x1y1y1 - gy1y1y1y1 - r * gy1y1;

	for (i = 0 ; i < NI ; ++i) {
		for (j = 0 ; j < NI ; ++j) {
			U0[i][j] = V[i][j];
		}
	}

	initGpt(tGpt1);
	//print3x3(tGpt1);
	for (loop = 0 ; loop < MAXNR ; ++loop) {

		V[0][8] = /* a11 */  tGpt1[0][0];
		V[1][8] = /* a12 */  tGpt1[0][1];
		V[2][8] = /* a21 */  tGpt1[1][0];
		V[3][8] = /* a22 */  tGpt1[1][1];
		V[4][8] = /* b1  */  tGpt1[0][2];
		V[5][8] = /* b2  */  tGpt1[1][2];
		V[6][8] = /* c1  */  tGpt1[2][0];
		V[7][8] = /* c2  */  tGpt1[2][1];
		multplyMV(V, v);

		detA = tGpt1[0][0] * tGpt1[1][1] - tGpt1[1][0] * tGpt1[0][1];
		Ainv11 =  tGpt1[1][1] / detA;
		Ainv21 = -tGpt1[1][0] / detA;
		Ainv12 = -tGpt1[0][1] / detA;
		Ainv22 =  tGpt1[0][0] / detA;
		grAinv11 = rg0 * Ainv11;
		grAinv21 = rg0 * Ainv21;
		grAinv12 = rg0 * Ainv12;
		grAinv22 = rg0 * Ainv22;

		for (i = 0 ; i < NI ; ++i) {
			for (j = 0 ; j < NI ; ++j) {
				U[i][j] = U0[i][j];
			}
		}

		U[0][0] += /* D11 + rT11 */     grAinv11 * Ainv11;
		U[1][0] += /* D21 + rT21 */     grAinv11 * Ainv12;
		U[2][0] += /* rT31 */           grAinv21 * Ainv11;
		U[3][0] += /* rT41 */           grAinv21 * Ainv12;

		U[0][1] += /* D12 + rT12 */     grAinv11 * Ainv21;
		U[1][1] += /* D22 + rT22 */     grAinv11 * Ainv22;
		U[2][1] += /* rT32 */           grAinv21 * Ainv21;
		U[3][1] += /* rT42 */           grAinv21 * Ainv22;

		U[0][2] += /* rT13 */           grAinv21 * Ainv11;
		U[1][2] += /* rT23 */           grAinv21 * Ainv12;
		U[2][2] += /* D11 + rT33 */     grAinv22 * Ainv11;
		U[3][2] += /* D21 + rT43 */     grAinv22 * Ainv12;

		U[0][3] += /* rT14 */           grAinv21 * Ainv21;
		U[1][3] += /* rT24 */           grAinv21 * Ainv22;
		U[2][3] += /* D12 + rT34 */     grAinv22 * Ainv21;
		U[3][3] += /* D22 + rT44 */     grAinv22 * Ainv22;

		U[0][8] = /* G'11 - v[0] */    gx1x2 + grAinv11 - v[0];
		U[1][8] = /* G'12 - v[1] */    gy1x2 + grAinv12 - v[1];
		U[2][8] = /* G'21 - v[2] */    gx1y2 + grAinv21 - v[2];
		U[3][8] = /* G'22 - v[3] */    gy1y2 + grAinv22 - v[3];
		U[4][8] = /* n1   - v[4] */    gx2 - v[4];
		U[5][8] = /* n2   - v[5] */    gy2 - v[5];
		U[6][8] = /* h1   - v[6] */    gx1x1x2 + gx1y1y2 - r * gx1 - v[6];
		U[7][8] = /* h2   - v[7] */    gx1y1x2 + gy1y1y2 - r * gy1 - v[7];

		//printf("gx2 = %lf gy2 = %lf  v[4] = %lf  v[5] = %lf \n", gx2, gy2, v[4], v[5]);
		solveLEq(U);

		tGpt1[0][0] += MU * U[0][8];
		tGpt1[0][1] += MU * U[1][8];
		tGpt1[1][0] += MU * U[2][8];
		tGpt1[1][1] += MU * U[3][8];
		tGpt1[0][2] += MU * U[4][8];
		tGpt1[1][2] += MU * U[5][8];
		tGpt1[2][0] += MU * U[6][8];
		tGpt1[2][1] += MU * U[7][8];

		// tGpt1[2][0] = 	tGpt1[2][1] = 0.0; /* Let c1 = c2 = 0.0 */
	}

	// print3x3(tGpt1);

	/* update of GAT components */
	multiply3x3(tGpt1, gpt, tGpt2);
	//print3x3(tGpt2);
	copyNormalGpt(tGpt2, gpt);
}

void gptcorsHoGInte(char sHoG1[ROW - 4][COL - 4], double g_can1[ROW][COL],
		                char sHoG2[ROW - 4][COL - 4], double g_can2[ROW][COL], double gwt[ROW][COL], double inteCanDir[ROWINTE][COLINTE][64],
		                double inteDx2Dir[ROWINTE][COLINTE][64], double inteDy2Dir[ROWINTE][COLINTE][64], double dnn, double gpt[3][3]) {
	/* determination of optimal GAT components */
	/* that yield the maximal correlation value */
	int x1, y1;
	int i, j, loop;
	double g0, gx1, gy1, gx2, gy2;
	double gx1x1, gx1y1, gy1y1, gx1x2, gx1y2, gy1x2, gy1y2;
	double t0, tx2, ty2;
	double gx1x1x1, gx1x1y1, gx1y1y1, gy1y1y1;
	double gx1x1x2, gx1x1y2, gx1y1x2, gx1y1y2, gy1y1x2, gy1y1y2;
	double gx1x1x1x1, gx1x1x1y1, gx1x1y1y1, gx1y1y1y1, gy1y1y1y1;
	double tx1, ty1, tx1x1, tx1y1, ty1y1, tx1x1x1, tx1x1y1, tx1y1y1, ty1y1y1;
	double U[NI][NI + 1], U0[NI][NI+1];
	double V[NI][NI + 1], v[NI];
	double dx1, dy1;
	double tGpt1[3][3], tGpt2[3][3];
	double detA, r, rg0;
	double Ainv11, Ainv12, Ainv21, Ainv22, grAinv11, grAinv12, grAinv21, grAinv22;
	double var = WGT * WGT * dnn * dnn;

	int ang1;
#ifdef TWOWINDOW
	int windowSSS = (int) (WGSS * WGT * dnn + 0.9999);
	int windowSSL = (int) (WGSL * WGT * dnn + 0.9999);
	int pPosS, mPosS, pPosL, mPosL;
	if (windowSSL > MAXWINDOWSIZE) {
		windowSSS = MAXWINDOWSIZE / 2;
		windowSSL = MAXWINDOWSIZE;
	}
	pPosS = MAXWINDOWSIZE + 1 + windowSSS;
	mPosS = MAXWINDOWSIZE - windowSSS;
	pPosL = MAXWINDOWSIZE + 1 + windowSSL;
	mPosL = MAXWINDOWSIZE - windowSSL;
#else
	int windowS = (int) (WGS * WGT * dnn + 0.9999);
	int pPos, mPos;
	if (windowS > MAXWINDOWSIZE) windowS = MAXWINDOWSIZE;
	pPos = MAXWINDOWSIZE + 1 + windowS;
	mPos = MAXWINDOWSIZE - windowS;
#endif

	/* Gaussian weigthed mean values */
	g0 = gx1 = gy1 = gx2 = gy2 = 0.0;
	gx1x1 = gx1y1 = gy1y1 = 0.0;
	gx1x2 = gx1y2 = gy1x2 = gy1y2 = 0.0;
	gx1x1x1 = gx1x1y1 = gx1y1y1 = gy1y1y1 = 0.0;
	gx1x1x2 = gx1x1y2 = gx1y1x2 = gx1y1y2 = gy1y1x2 = gy1y1y2 = 0.0;
	gx1x1x1x1 = gx1x1x1y1 = gx1x1y1y1 = gx1y1y1y1 = gy1y1y1y1 = 0.0;

	for (y1 = MARGINE + 2 ; y1 < ROW - MARGINE - 2 ; y1++) {
		dy1 = y1 - CY;
		for (x1 = MARGINE + 2 ; x1 < COL - MARGINE - 2 ; x1++) {
			dx1 = x1 - CX;

			ang1 = sHoG2Idx(sHoG1[y1 - 2][x1 - 2]);
      if (ang1 == -1) continue;
      // printf("ang1 = %d\n", ang1);
#ifdef TWOWINDOW
			t0 = tx2 = ty2 = COEFGS;
			t0  *= inteCanDir[y1 + pPosL][x1 + pPosL][ang1] - inteCanDir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteCanDir[y1 + mPosL][x1 + pPosL][ang1] + inteCanDir[y1 + mPosL][x1 + mPosL][ang1];
			tx2 *= inteDx2Dir[y1 + pPosL][x1 + pPosL][ang1] - inteDx2Dir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteDx2Dir[y1 + mPosL][x1 + pPosL][ang1] + inteDx2Dir[y1 + mPosL][x1 + mPosL][ang1];
			ty2 *= inteDy2Dir[y1 + pPosL][x1 + pPosL][ang1] - inteDy2Dir[y1 + pPosL][x1 + mPosL][ang1]
					 - inteDy2Dir[y1 + mPosL][x1 + pPosL][ang1] + inteDy2Dir[y1 + mPosL][x1 + mPosL][ang1];
			t0  += inteCanDir[y1 + pPosS][x1 + pPosS][ang1] - inteCanDir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteCanDir[y1 + mPosS][x1 + pPosS][ang1] + inteCanDir[y1 + mPosS][x1 + mPosS][ang1];
			tx2 += inteDx2Dir[y1 + pPosS][x1 + pPosS][ang1] - inteDx2Dir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteDx2Dir[y1 + mPosS][x1 + pPosS][ang1] + inteDx2Dir[y1 + mPosS][x1 + mPosS][ang1];
			ty2 += inteDy2Dir[y1 + pPosS][x1 + pPosS][ang1] - inteDy2Dir[y1 + pPosS][x1 + mPosS][ang1]
					 - inteDy2Dir[y1 + mPosS][x1 + pPosS][ang1] + inteDy2Dir[y1 + mPosS][x1 + mPosS][ang1];
#else
			t0  = inteCanDir[y1 + pPos][x1 + pPos][ang1] - inteCanDir[y1 + pPos][x1 + mPos][ang1]
					 - inteCanDir[y1 + mPos][x1 + pPos][ang1] + inteCanDir[y1 + mPos][x1 + mPos][ang1];
			tx2 = inteDx2Dir[y1 + pPos][x1 + pPos][ang1] - inteDx2Dir[y1 + pPos][x1 + mPos][ang1]
					 - inteDx2Dir[y1 + mPos][x1 + pPos][ang1] + inteDx2Dir[y1 + mPos][x1 + mPos][ang1];
			ty2 = inteDy2Dir[y1 + pPos][x1 + pPos][ang1] - inteDy2Dir[y1 + pPos][x1 + mPos][ang1]
					 - inteDy2Dir[y1 + mPos][x1 + pPos][ang1] + inteDy2Dir[y1 + mPos][x1 + mPos][ang1];
#endif
			t0 *= g_can1[y1][x1]; tx2 *= g_can1[y1][x1]; ty2 *= g_can1[y1][x1];


			g0    += t0;
			gx1   += tx1 = t0  * dx1;
			gy1   += ty1 = t0  * dy1;
			gx1x1 += tx1x1 = tx1 * dx1;
			gx1y1 += tx1y1 = tx1 * dy1;
			gy1y1 += ty1y1 = ty1 * dy1;
			gx1x1x1 += tx1x1x1 = tx1x1 * dx1;
			gx1x1y1 += tx1x1y1 = tx1x1 * dy1;
			gx1y1y1 += tx1y1y1 = tx1y1 * dy1;
			gy1y1y1 += ty1y1y1 = ty1y1 * dy1;
			gx1x1x1x1 += tx1x1x1 * dx1;
			gx1x1x1y1 += tx1x1x1 * dy1;
			gx1x1y1y1 += tx1x1y1 * dy1;
			gx1y1y1y1 += tx1y1y1 * dy1;
			gy1y1y1y1 += ty1y1y1 * dy1;

			gx2   += tx2;
			gy2   += ty2;
			gx1x2 += tx2 * dx1;
			gx1y2 += ty2 * dx1;
			gy1x2 += tx2 * dy1;
			gy1y2 += ty2 * dy1;
			gx1x1x2 += tx2 * dx1 * dx1;
			gx1x1y2 += ty2 * dx1 * dx1;
			gx1y1x2 += tx2 * dx1 * dy1;
			gx1y1y2 += ty2 * dx1 * dy1;
			gy1y1x2 += tx2 * dy1 * dy1;
			gy1y1y2 += ty2 * dy1 * dy1;
		}
	}
	 // printf("dnn = %f g0 = %f gx1 = %f  gy1 = %f gx2 = %f  gy2 = %f\n", dnn, g0, gx1, gy1, gx2, gy2);
	 // printf("dnn = %f normal gx1x1x1x1 = %f gx1x1x1y1 = %f  gx1x1y1y1 = %f gx1y1y1y1 = %f  gy1y1y1y1y1 = %f\n", dnn, gx1x1x1x1/g0, gx1x1x1y1/g0, gx1x1y1y1/g0, gx1y1y1y1/g0, gy1y1y1y1/g0);

	if (g0 == 0.0) {
#ifdef PRINT
		printf("GAT calculation failure by zero sum!!!\n");
#endif
	}
	/* Newton method */
	r = 0.5 * var;
	rg0 = r * g0;

	V[0][0] = /* D11 */  gx1x1;
	V[1][0] = /* D21 */  gx1y1;
	V[2][0] = /* 0   */  0.0;
	V[3][0] = /* 0   */  0.0;
	V[4][0] = /* m1  */  gx1;
	V[5][0] = /* 0   */  0.0;
	V[6][0] = /* E11 */  gx1x1x1;
	V[7][0] = /* E31 */  gx1x1y1;
	V[0][1] = /* D12 */  gx1y1;
	V[1][1] = /* D22 */  gy1y1;
	V[2][1] = /* 0   */  0.0;
	V[3][1] = /* 0   */  0.0;
	V[4][1] = /* m2  */  gy1;
	V[5][1] = /* 0   */  0.0;
	V[6][1] = /* E12 */  gx1x1y1;
	V[7][1] = /* E32 */  gx1y1y1;
	V[0][2] = /* 0   */  0.0;
	V[1][2] = /* 0   */  0.0;
	V[2][2] = /* D11 */  gx1x1;
	V[3][2] = /* D21 */  gx1y1;
	V[4][2] = /* 0   */  0.0;
	V[5][2] = /* m1  */  gx1;
	V[6][2] = /* E21 */  gx1x1y1;
	V[7][2] = /* E41 */  gx1y1y1;
	V[0][3] = /* 0   */  0.0;
	V[1][3] = /* 0   */  0.0;
	V[2][3] = /* D12 */  gx1y1;
	V[3][3] = /* D22 */  gy1y1;
	V[4][3] = /* 0   */  0.0;
	V[5][3] = /* m2  */  gy1;
	V[6][3] = /* E22 */  gx1y1y1;
	V[7][3] = /* E42 */  gy1y1y1;
	V[0][4] = /* m1  */  gx1;
	V[1][4] = /* m2  */  gy1;
	V[2][4] = /* 0   */  0.0;
	V[3][4] = /* 0   */  0.0;
	V[4][4] = /* 1   */  g0;
	V[5][4] = /* 0   */  0.0;
	V[6][4] = /* D11 */  gx1x1;
	V[7][4] = /* D21 */  gx1y1;
	V[0][5] = /* 0   */  0.0;
	V[1][5] = /* 0   */  0.0;
	V[2][5] = /* m1  */  gx1;
	V[3][5] = /* m2  */  gy1;
	V[4][5] = /* 0   */  0.0;
	V[5][5] = /* 1   */  g0;
	V[6][5] = /* D12 */  gx1y1;
	V[7][5] = /* D22 */  gy1y1;

	V[0][6] = /*-E11 */ -gx1x1x1;
	V[1][6] = /*-E21 */ -gx1x1y1;
	V[2][6] = /*-E31 */ -gx1x1y1;
	V[3][6] = /*-E41 */ -gx1y1y1;
	V[4][6] = /*-D11 */ -gx1x1;
	V[5][6] = /*-D21 */ -gx1y1;
	V[6][6] = /*-F11 - r * D11 */  -gx1x1x1x1 - gx1x1y1y1 - r * gx1x1;
	V[7][6] = /*-F21 - r * D21 */  -gx1x1x1y1 - gx1y1y1y1 - r * gx1y1;

	V[0][7] = /*-E12 */ -gx1x1y1;
	V[1][7] = /*-E22 */ -gx1y1y1;
	V[2][7] = /*-E32 */ -gx1y1y1;
	V[3][7] = /*-E42 */ -gy1y1y1;
	V[4][7] = /*-D12 */ -gx1y1;
	V[5][7] = /*-D22 */ -gy1y1;
	V[6][7] = /*-F12 - r * D12 */ -gx1x1x1y1 - gx1y1y1y1 - r * gx1y1;
	V[7][7] = /*-F22 - r * D22 */ -gx1x1y1y1 - gy1y1y1y1 - r * gy1y1;

	for (i = 0 ; i < NI ; ++i) {
		for (j = 0 ; j < NI ; ++j) {
			U0[i][j] = V[i][j];
		}
	}

	initGpt(tGpt1);
	//print3x3(tGpt1);
	for (loop = 0 ; loop < MAXNR ; ++loop) {

		V[0][8] = /* a11 */  tGpt1[0][0];
		V[1][8] = /* a12 */  tGpt1[0][1];
		V[2][8] = /* a21 */  tGpt1[1][0];
		V[3][8] = /* a22 */  tGpt1[1][1];
		V[4][8] = /* b1  */  tGpt1[0][2];
		V[5][8] = /* b2  */  tGpt1[1][2];
		V[6][8] = /* c1  */  tGpt1[2][0];
		V[7][8] = /* c2  */  tGpt1[2][1];
		multplyMV(V, v);

		detA = tGpt1[0][0] * tGpt1[1][1] - tGpt1[1][0] * tGpt1[0][1];
		Ainv11 =  tGpt1[1][1] / detA;
		Ainv21 = -tGpt1[1][0] / detA;
		Ainv12 = -tGpt1[0][1] / detA;
		Ainv22 =  tGpt1[0][0] / detA;
		grAinv11 = rg0 * Ainv11;
		grAinv21 = rg0 * Ainv21;
		grAinv12 = rg0 * Ainv12;
		grAinv22 = rg0 * Ainv22;

		for (i = 0 ; i < NI ; ++i) {
			for (j = 0 ; j < NI ; ++j) {
				U[i][j] = U0[i][j];
			}
		}

		U[0][0] += /* D11 + rT11 */     grAinv11 * Ainv11;
		U[1][0] += /* D21 + rT21 */     grAinv11 * Ainv12;
		U[2][0] += /* rT31 */           grAinv21 * Ainv11;
		U[3][0] += /* rT41 */           grAinv21 * Ainv12;

		U[0][1] += /* D12 + rT12 */     grAinv11 * Ainv21;
		U[1][1] += /* D22 + rT22 */     grAinv11 * Ainv22;
		U[2][1] += /* rT32 */           grAinv21 * Ainv21;
		U[3][1] += /* rT42 */           grAinv21 * Ainv22;

		U[0][2] += /* rT13 */           grAinv21 * Ainv11;
		U[1][2] += /* rT23 */           grAinv21 * Ainv12;
		U[2][2] += /* D11 + rT33 */     grAinv22 * Ainv11;
		U[3][2] += /* D21 + rT43 */     grAinv22 * Ainv12;

		U[0][3] += /* rT14 */           grAinv21 * Ainv21;
		U[1][3] += /* rT24 */           grAinv21 * Ainv22;
		U[2][3] += /* D12 + rT34 */     grAinv22 * Ainv21;
		U[3][3] += /* D22 + rT44 */     grAinv22 * Ainv22;

		U[0][8] = /* G'11 - v[0] */    gx1x2 + grAinv11 - v[0];
		U[1][8] = /* G'12 - v[1] */    gy1x2 + grAinv12 - v[1];
		U[2][8] = /* G'21 - v[2] */    gx1y2 + grAinv21 - v[2];
		U[3][8] = /* G'22 - v[3] */    gy1y2 + grAinv22 - v[3];
		U[4][8] = /* n1   - v[4] */    gx2 - v[4];
		U[5][8] = /* n2   - v[5] */    gy2 - v[5];
		U[6][8] = /* h1   - v[6] */    gx1x1x2 + gx1y1y2 - r * gx1 - v[6];
		U[7][8] = /* h2   - v[7] */    gx1y1x2 + gy1y1y2 - r * gy1 - v[7];

		//printf("gx2 = %lf gy2 = %lf  v[4] = %lf  v[5] = %lf \n", gx2, gy2, v[4], v[5]);
		solveLEq(U);

		tGpt1[0][0] += MU * U[0][8];
		tGpt1[0][1] += MU * U[1][8];
		tGpt1[1][0] += MU * U[2][8];
		tGpt1[1][1] += MU * U[3][8];
		tGpt1[0][2] += MU * U[4][8];
		tGpt1[1][2] += MU * U[5][8];
		tGpt1[2][0] += MU * U[6][8];
		tGpt1[2][1] += MU * U[7][8];

		// tGpt1[2][0] = 	tGpt1[2][1] = 0.0; /* Let c1 = c2 = 0.0 */
	}

	// print3x3(tGpt1);

	/* update of GAT components */
	multiply3x3(tGpt1, gpt, tGpt2);
	//print3x3(tGpt2);
	copyNormalGpt(tGpt2, gpt);
}
