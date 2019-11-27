/* make */
/* ./execGpt 9 11 */
/* ./execGpt 7 3 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "include/parameter.h"
#include "include/utility.h"
#include "include/stdGpt.h"
#include "include/acclGpt.h"
#include "include/stdInte.h"

/* Image storage arrays */
unsigned char image1[MAX_IMAGESIZE][MAX_IMAGESIZE];
unsigned char image2[MAX_IMAGESIZE][MAX_IMAGESIZE];
int x_size1 = COL, y_size1 = ROW; /* width & height of image1*/
int x_size2, y_size2; /* width & height of image2 */

int    inteAng[ROWINTE][COLINTE][9];
double inteCanDir[ROWINTE][COLINTE][9];
double inteDx2Dir[ROWINTE][COLINTE][9];
double inteDy2Dir[ROWINTE][COLINTE][9];

int    inteAng64[ROWINTE][COLINTE][64];
double inteCanDir64[ROWINTE][COLINTE][64];
double inteDx2Dir64[ROWINTE][COLINTE][64];
double inteDy2Dir64[ROWINTE][COLINTE][64];

double H1[ROW][COL * 162], Ht1[ROW][COL * 27];
double H2[ROW - 4][(COL - 4) * 6 * 64 * 3], Ht2[ROW - 4][(COL - 4) * 64 * 3];
double H3[ROW - 4][(COL - 4) * 6 * 64 * 3], Ht3[ROW - 4][(COL - 4) * 64 * 3];
double D1[ROW][COL * 8];
double D2[ROW - 4][(COL - 4) * 64];

double ndis[(2 * ROW - 1) * (2 * COL - 1)];
int coor[(2 * ROW - 1) * (2 * COL - 1)][2];

int main(int argc, char* argv[]) {
	if (argc != 1) {
#undef	MATCHMETHOD
#undef  DISTANCETYPE
#define MATCHMETHOD atoi(argv[1])
#define DISTANCETYPE atoi(argv[2])
	}
  int image3[ROW2][COL2], image4[ROW][COL];					// image3: test image	image4: training image
	int margine = CANMARGIN / 2;
	int x, y, iter;
	char csvname[MAX_FILENAME], foldername[MAX_FILENAME];	// GAT, NGAT, GPT, NGPT, name of .csv file, foldername
	double gk[ROW][COL], gwt[ROW][COL], dnn, temp_dnn, var;			// Gaussian window initial, Gaussian window, window size, variance
	int g_ang1[ROW][COL], g_ang2[ROW][COL];					// direction of gradients
	char g_HoG1[ROW][COL][8], g_HoG2[ROW][COL][8];			// HoG feature of the images
	char sHoG1[ROW - 4][COL - 4], sHoG2[ROW - 4][COL - 4];
	double g_nor1[ROW][COL], g_nor2[ROW][COL];				// norm of gradients
	double g_can1[ROW][COL], g_can2[ROW][COL];				// canonicalized images
	double g_can11[ROW - CANMARGIN][COL - CANMARGIN], g_can22[ROW - CANMARGIN][COL - CANMARGIN];
															// canonicalized images center
	double old_cor0, old_cor1, new_cor1;					//
	double org_cor, gat_corf, gat_corb;
	double gpt0[3][3], gpt1[3][3], gptInv[3][3];

	clock_t start, end, start0, end0, end1;
	double elapse, elapse1 = 0.0, elapse2 = 0.0;

	char fileName[128];
#ifdef SAVEDATAFORMATLAB
    // sprintf(fileName, "%s/%s_%d_%d.csv", IMGDIR, TsIMAGE, MATCHMETHOD, DATATYPE);
    sprintf(fileName, "%s/%s_%d_%d_%d.csv", IMGDIR, TsIMAGE, MATCHMETHOD, DATATYPE, DISTANCETYPE);
    FILE *fp;
    if((fp = fopen(fileName, "w")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
#endif

/* Initialize the GPT matrix */
initGpt(gpt0);
initGpt2(gpt1, ZOOM, ZOOM * BETA, B1, B2, ROT);

/* initialize Gauss window function */
for (y = 0; y < ROW; y++)
	for (x = 0; x < COL; x++)
		gk[y][x] = exp(-(x*x+y*y)/2.0);

/* Load template image and save it to image4, the local memory */
sprintf(fileName, "%s/%s.pgm", IMGDIR, RgIMAGE);
load_image_file(fileName, image1, COL, ROW);
for (y = 0; y < ROW; y++)
	for (x = 0; x < COL; x++)
		image4[y][x] = image1[y][x];
procImg(g_can2, g_ang2, g_nor2, g_HoG2, sHoG2, image1);

start = clock();
calInte(g_can2, g_ang2, inteAng, inteCanDir, inteDx2Dir, inteDy2Dir);
calInte64(g_can2, sHoG2, inteAng64, inteCanDir64, inteDx2Dir64, inteDy2Dir64);
end   = clock();
printf("time of making integral images = %f \n", (double)(end - start) / CLOCKS_PER_SEC);

/* Make template tables if required */

#if MAKETEMP == 1
	start = clock();
	sprintf(fileName, "%s/%s", IMGDIR, RgIMAGE);
	makeTemp(g_ang2, g_can2, gk, H1, fileName);
	makeTemp64(sHoG2, g_can2, gk, H2, fileName);
	makeTemp64_far(sHoG2, g_can2, gk, H3, fileName);
	winTbl(g_ang2, D1, fileName);
	winTbl64(sHoG2, D2, fileName);
	searchTbl(ROW, COL, fileName);
	end   = clock();
	printf("time of making templates = %f \n", (double)(end - start) / CLOCKS_PER_SEC);
	return 0;
#elif MAKETEMP == 0
	loadTbls(D1, D2, ndis, coor);
	loadTemp(H1);
	loadTemp64(H2);
	loadTemp64_far(H3);
#endif

/* Load test image and save it to image3, the local memory */
sprintf(fileName, "%s/%s.pgm", IMGDIR, TsIMAGE);
load_image_file(fileName, image1, COL2, ROW2);
for (y = 0; y < ROW2; y++)
	for (x = 0; x < COL2; x++)
		image3[y][x] = image1[y][x];

/* save the initial image */
for (y = 0; y < ROW2; y++)
	for (x = 0; x < COL2; x++)
		image2[y][x] = image1[y][x];
bilinear_normal_projection(gpt1, COL, ROW, COL2, ROW2, image1, image2);
sprintf(fileName, "%s/%s_init.pgm", IMGDIR, RgIMAGE);
save_image_file(fileName, image2, COL, ROW);
procImg(g_can1, g_ang1, g_nor1, g_HoG1, sHoG1, image2);

/***************Pre-setting finish***************/

/* calculate the initial correlation */
old_cor1 = 0.0;
for (y = margine ; y < ROW - margine ; y++) 
	for (x = margine ; x < COL - margine ; x++) 
		old_cor1 += g_can1[y][x] * g_can2[y][x];
org_cor = old_cor1;
printf("Original cor. = %f\n", org_cor);
old_cor0 = old_cor1;

/* calculate the initial dnn */

double d2 = 0.0;

switch (DISTANCETYPE) {
case 0:
	dnn = winpat(g_ang1, g_ang2);
	if (dnn > DNNSWITCHTHRE)
		dnn = sHoGpat(sHoG1, sHoG2);
	break;
case 1:
	dnn = winpat(g_ang1, g_ang2);
	break;
case 2:
	dnn = fwinpat(g_ang1, g_ang2, D1, ndis, coor);
	break;
case 3:
	dnn = sHoGpat(sHoG1, sHoG2);
	// d2  =
	break;
case 4:
	dnn = fsHoGpat(sHoG1, sHoG2, D2, ndis, coor);
	break;
case 5:
	dnn = WNNDEGD * winpatInte(g_ang1, inteAng);
	break;
case 6:
	dnn = WNNDEsHoGD * sHoGpatInte(sHoG1, inteAng64);
	break;
case 10:
	dnn = fsHoGpat(sHoG1, sHoG2, D2, ndis, coor);
	if (dnn <= DNNSWITCHTHRE)
		dnn = fwinpat(g_ang1, g_ang2, D1, ndis, coor);
	break;
case 11:
	dnn = WNNDEGD * winpatInte(g_ang1, inteAng);
	if (dnn > DNNSWITCHTHRE)
		dnn = WNNDEsHoGD * sHoGpatInte(sHoG1, inteAng64);
	break;
case 12:
	// dnn = sHoGpatInte(sHoG1, inteAng64) * 1.2;
	break;
}

/***************Main iteration loop*************/
/* lap the start time */
start = clock();
for (iter = 0 ; iter < MAXITER ; iter++) {
	if (iter >= 100)
	{
		fprintf(fp, "%d,%f,%f,%f,%f\n", iter, new_cor1, dnn, 1 / var, (double)(end1 - start) / CLOCKS_PER_SEC);
		continue;
	}
	/* update gauss window function */
	var = pow(WGT * dnn, 2);
	for (y = 0; y < ROW; y++)
		for (x = 0; x < COL; x++)
			gwt[y][x] = pow(gk[y][x], 1.0 / var);

	/* select matching method */
	start0 = clock();
	// switch (MATCHMETHOD) {
	switch (MATCHMETHOD) {
	case 1:
	case 2:
		if (iter <= 5 || iter % 2 == 0)
			gatcor(g_ang1, g_can1, g_ang2, g_can2, gwt, gpt1);
		else
			pptcor(g_ang1, g_can1, g_ang2, g_can2, gwt, gpt1);
		break;
	case 5:
		sgptcor(g_ang1, g_can1, g_ang2, g_can2, gwt, gpt1);
		break;
	case 6:
		nsgptcor(g_ang1, g_can1, g_ang2, g_can2, gwt, gpt1, dnn);
		break;
	case 7:
		nsgptcorSpHOG5x5(g_ang1, sHoG1, g_can1,	g_ang2, sHoG2, g_can2, gwt, gpt1, dnn);
		break;
	case 8:
		gptcorInte(g_ang1, g_can1, g_ang2, g_can2, gwt, inteCanDir, inteDx2Dir, inteDy2Dir, dnn, gpt1);
		break;
	case 9:
		gptcorsHoGInte(sHoG1, g_can1, sHoG2, g_can2, gwt, inteCanDir64, inteDx2Dir64, inteDy2Dir64, dnn, gpt1);
		break;
	case 16:
		fnsgptcor(g_ang1, g_can1, gpt1, dnn, H1, Ht1);
		break;
	case 17:
		if (dnn > 8.0)
			fnsgptcorSpHOG5x5_far(g_ang1, sHoG1, g_can1, gpt1, dnn, H3, Ht3);
		else
			fnsgptcorSpHOG5x5(g_ang1, sHoG1, g_can1, gpt1, dnn, H2, Ht2);
	}
	end0 = clock();
	elapse1 += (double)(end0 - start0) / CLOCKS_PER_SEC;

	/* transform the test image and update g_can1, g_ang1, g_nor1, g_HoG1, sHoG1 */
	for (y = 0; y < ROW2; y++)
		for (x = 0; x < COL2; x++)
			image1[y][x] = (unsigned char)image3[y][x];
	bilinear_normal_projection(gpt1, COL, ROW, COL2, ROW2, image1, image2);
	procImg(g_can1, g_ang1, g_nor1, g_HoG1, sHoG1, image2);

	/* update correlation */
	new_cor1 = 0.0;
	for (y = margine ; y < ROW - margine ; y++) 
		for (x = margine ; x < COL - margine ; x++) 
			new_cor1 += g_can1[y][x] * g_can2[y][x];

	/* Calculation distance */
	start0 = clock();
	switch (DISTANCETYPE) {
	case 0:
		dnn = winpat(g_ang1, g_ang2);
		if (dnn > DNNSWITCHTHRE)
				dnn = sHoGpat(sHoG1, sHoG2);
		break;
	case 1:
		dnn = winpat(g_ang1, g_ang2);
		break;
	case 2:
		dnn = fwinpat(g_ang1, g_ang2, D1, ndis, coor);
		break;
	case 3:
		dnn = sHoGpat(sHoG1, sHoG2);
		break;
	case 4:
		dnn = fsHoGpat(sHoG1, sHoG2, D2, ndis, coor);
		break;
	case 5:
		dnn = WNNDEGD * winpatInte(g_ang1, inteAng);
		break;
	case 6:
		dnn = WNNDEsHoGD * sHoGpatInte(sHoG1, inteAng64);
		break;
	case 10:
		dnn = fwinpat(g_ang1, g_ang2, D1, ndis, coor);
		if (dnn > DNNSWITCHTHRE)
			dnn = fsHoGpat(sHoG1, sHoG2, D2, ndis, coor);
		break;
	case 11:
		dnn = WNNDEGD * winpatInte(g_ang1, inteAng);
		if (dnn > DNNSWITCHTHRE)
			dnn = WNNDEsHoGD * sHoGpatInte(sHoG1, inteAng64);
		break;
	case 12:
		// dnn = sHoGpatInte(sHoG1, inteAng64) * 1.2;
		break;
	}
	end0 = clock();
	elapse2 = (double)(end0 - start) / CLOCKS_PER_SEC;
	end1 = clock();
	/* display message */
	printf("iter = %d, new col. = %f dnn = %f  var = %f (d2 = %f) (time = %f)\n", iter, new_cor1, dnn, 1 / var, d2, elapse2);
#ifdef SAVEDATAFORMATLAB
			fprintf(fp, "%d,%f,%f,%f,%f\n", iter, new_cor1, dnn, 1 / var, (double)(end1 - start) / CLOCKS_PER_SEC);
#endif

}
/* display the calculation time */
end = clock();
elapse = (double)(end - start) / CLOCKS_PER_SEC;
printf("\nelapsed time = %.3f sec\nelapse1 = %.3f, elapse2 = %.3f \n", elapse, elapse1, elapse2);
#ifdef SAVEDATAFORMATLAB
	fprintf(fp, "%f,%f,%f\n", gpt1[0][0], gpt1[0][1], gpt1[0][2]);
	fprintf(fp, "%f,%f,%f\n", gpt1[1][0], gpt1[1][1], gpt1[1][2]);
	fprintf(fp, "%f,%f,%f\n", gpt1[2][0], gpt1[2][1], gpt1[2][2]);
	fprintf(fp, "%f\n", elapse);
	fclose(fp);
#endif
sprintf(fileName, "%s/%s_out.pgm", IMGDIR, RgIMAGE);
printf("%s\n", fileName);
save_image_file(fileName, image2, COL, ROW);

	printf("OK\n");
}
