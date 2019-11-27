#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "parameter.h"
#include "utility.h"
#include"accGpt.h"

void winTbl(int g_ang2[ROW][COL], double D[ROW][COL * 8], char *);
// void winTbl64(char sHoG2[ROW - 4][COL - 4], char *);
void searchTbl(int row, int col, char *);
void loadTbls(double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]);
// void loadTbls64(double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]);
// void loadTbls64_far(double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]);
void makeTemp(int g_ang2[ROW][COL], double g_can2[ROW][COL], double gwt[ROW][COL], double H[ROW][COL * 162], char *);
// void makeTemp64(char sHoG2[ROW - 4][COL - 4], double g_can2[ROW][COL], double gwt[ROW][COL], double H[ROW - 4][(COL - 4) * 6 * 64 * 3], char *);
// void makeTemp64_far(char sHoG2[ROW - 4][COL - 4], double g_can2[ROW][COL], double gwt[ROW][COL], double H[ROW - 4][(COL - 4) * 6 * 64 * 3], char *);

double fwinpat(int g_ang1[ROW][COL], int g_ang2[ROW][COL], double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]);
// double fsHoGpat(char sHoG1[ROW - 4][COL - 4], char sHoG2[ROW - 4][COL - 4], double D[ROW][COL * 64], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]);

void fgatcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fngatcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fpptcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fnpptcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fsgptcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fnsgptcor(int g_ang1[ROW][COL], double g_can1[ROW][COL], double gpt[3][3], double dnn, Ht[ROW][COL * 27], char *);
void fnsgptcor5x5Spl(int g_ang1[ROW][COL], char sHoG1[ROW - 4][COL - 4], double g_can1[ROW][COL],
                     double gpt[3][3], double dnn, Ht[ROW - 4][(COL - 4) * 64 * 3] char *);

/****************************************************************/
void winTbl(int g_ang2[ROW][COL], double D[ROW][COL * 8], char *) {
    char fnt[128];
    int x2, y2, tx2, ty2, s;
    double D[ROW][COL * 8];
    double minInit, min, delta;
    
    sprintf(fnt, "%sDTbl", fn);
    
    minInit = (ROW - 2 * MARGINE) * (ROW - 2 * MARGINE) + (COL - 2 * MARGINE) * (COL - 2 * MARGINE);
    for (s = 0 ; s < 8 ; s++) {
        for (y2 = MARGINE ; y2 < ROW - MARGINE; y2++) {
            for (x2 = MARGINE ; x2 < COL - MARGINE ; x2++) {
                min = minInit;
                for (ty2 = MARGINE ; ty2 < ROW - MARGINE ; ty2++) {
                    for (tx2 = MARGINE ; tx2 < COL - MARGINE ; tx2++) {
                        if (g_ang2[ty2][tx2] == s) {
                            delta = (y2 - ty2) * (y2 - ty2) + (x2 - tx2) * (x2 - tx2);
                            if (delta < min) min = delta;
                        }
                    }
                }
                D[y2][x2 + s * COL] = sqrt(min);
            }
        }
    }
    
    /* Make file pointer */
    FILE *fp;
    if((fp = fopen(fnt, "wb")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    fwrite(D,  sizeof(double), COL * ROW * 8, fp);
    fclose(fp);
    
    sprintf(fnt, "%s.csv", fnt);
    if((fp = fopen(fnt, "w")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    for (y2 = 0 ; y2 < ROW ; y2++) {
        for (x2 = 0 ; x2 < COL ; x2++) {
            fprintf(fp, "%d,", g_ang2[y2][x2]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (y2 = 0 ; y2 < ROW ; y2++) {
        for (x2 = 0 ; x2 < 8 * COL ; x2++) {
            fprintf(fp, "%f,", D[y2][x2]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("Finish \n");
}

void searchTbl(int row, int col, char *fn) {
    int coor[(2 * row - 1) * (2 * col - 1)][2], ix, iy, count;
    double ndis[(2 * row - 1) * (2 * col - 1)];
    
    
    count = 0;
    /* submit */
    for (iy = 0 ; iy < 2 * row - 1 ; iy++) {
        for (ix = 0 ; ix < 2 * col - 1 ; ix++) {
            // yCoor[iy][ix] = iy - row + 1;
            // xCoor[iy][ix] = ix - col + 1;
            coor[count][0] = iy - row + 1;
            coor[count][1] = ix - col + 1;
            ndis[count]  = sqrt(coor[count][0] * coor[count][0] + coor[count][1] * coor[count][1]);
            count++;
        }
    }
    
    int indx[(2 * row - 1) * (2 * col - 1)], tmp[(2 * row - 1) * (2 * col - 1)][2];
    for (ix = 0 ; ix < (2 * row - 1) * (2 * col - 1) ; ix++) indx[ix] = ix;
    quickSortAsc(ndis, indx, 0, (2 * row - 1) * (2 * col - 1) - 1);
    
    for (ix = 0 ; ix < (2 * row - 1) * (2 * col - 1) ; ix++) {
        tmp[ix][0] = coor[indx[ix]][0];
        tmp[ix][1] = coor[indx[ix]][1];
    }
    
    for (ix = 0 ; ix < (2 * row - 1) * (2 * col - 1) ; ix++) {
        coor[ix][0] = tmp[ix][0];
        coor[ix][1] = tmp[ix][1];
    }
    
    /* Make file pointer */
    FILE *fp;
    if((fp = fopen(fn, "wb")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    fwrite(coor,  sizeof(int), (2 * row - 1) * (2 * col - 1) * 2, fp);
    fwrite(ndis,  sizeof(double), (2 * row - 1) * (2 * col - 1), fp);
    
    
    fclose(fp);
    
    sprintf(fn, "%s.csv", fn);
    if((fp = fopen(fn, "w")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    
    for (ix = 0 ; ix < 2 ; ix++) {
        for (count = 0 ; count < (2 * row - 1) * (2 * col - 1) ; count++) {
            fprintf(fp, "%d,", coor[count][ix]);
        }
        fprintf(fp, "\n");
    }
    
    for (count = 0 ; count < (2 * row - 1) * (2 * col - 1) ; count++) {
        fprintf(fp, "%f,", ndis[count]);
    }
    
    fclose(fp);
}

void loadTbls(double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]) {
    /* Load table */
    char fn[128];
    sprintf(fn, "%s/%s_tempDTbl", IMGDIR, RgIMAGE);
    FILE *fp;
    if((fp = fopen(fn, "rb")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    fread(D, sizeof(double), ROW * COL * 8, fp);
    fclose(fp);
    
    sprintf(fn, "%s/%s_SearchTbl", IMGDIR, RgIMAGE);
    if((fp = fopen(fn, "rb")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    fread(coor, sizeof(int), (2 * ROW - 1) * (2 * COL - 1) * 2, fp);
    fread(ndis, sizeof(double), (2 * ROW - 1) * (2 * COL - 1), fp);
    fclose(fp);
}

void makeTemp(int g_ang2[ROW][COL], double g_can2[ROW][COL], double gwt[ROW][COL], double H[ROW][COL * 162], char *fn) {
    int x1, y1, x2, y2, s, i, ix, iy, thre1, thre2;
    double tv, t0, tx2, ty2;
    double dx1, dy1, dx2, dy2;
    double gwt[ROW][COL];
    int count = 0;
    char fnt[128];
    sprintf(fnt, "%s", fn);
    
    // double H0[ROW][COL], H1x[ROW][COL], H1y[ROW][COL];
    // double H[ROW][COL * 162];
    
    /* Loop for making template */
    for (i = 0 ; i < 6 ; i++) {
        thre1 = 27 * i * COL;
        /* update gauss */
        for (iy = 0; iy < ROW; iy++) {
            for (ix = 0; ix < COL; ix++) {
                gwt[iy][ix] = pow(gt[iy][ix], VARTABLE[i]);
            }
        }
        
        for (s = -1 ; s < 8 ; s++) {
            thre2 = (s + 1) * 3 * COL;
            
            for (y1 = MARGINE ; y1 < ROW - MARGINE ; y1++) {
                dy1 = y1 - CY;
                for (x1 = MARGINE ; x1 < COL - MARGINE ; x1++) {
                    dx1 = x1 - CX;
                    t0  = 0.0; tx2 = 0.0; ty2 = 0.0;
                    for (y2 = MARGINE ; y2 < ROW - MARGINE ; y2++) {
                        dy2 = y2 - CY;
                        for (x2 = MARGINE ; x2 < COL - MARGINE ; x2++) {
                            dx2 = x2 - CX;
                            
                            if (s == g_ang2[y2][x2]) {
                                tv     = gwt[abs(y2 - y1)][abs(x2 - x1)] * g_can2[y2][x2];
                                t0    += tv;
                                tx2   += tv * dx2;
                                ty2   += tv * dy2;
                            }
                            
                        }
                    }
                    
                    H[y1][thre1 + thre2 + x1]           = t0;
                    H[y1][thre1 + thre2 + COL + x1]     = tx2;
                    H[y1][thre1 + thre2 + 2 * COL + x1] = ty2;
                    
                }
            }
            
            /* write file */
            count++;
            printf("thre1 + thre2 = %d\n", thre1 + thre2);
            /*
             fwrite(H0,  sizeof(double), COL * ROW, fp);
             fwrite(H1x, sizeof(double), COL * ROW, fp);
             fwrite(H1y, sizeof(double), COL * ROW, fp);
             */
            
        }
    }
    printf("All %d times\n", count);
    /* Make file pointer */
    FILE *fp;
    if((fp = fopen(fnt, "wb")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    fwrite(H,  sizeof(double), COL * ROW * 162, fp);
    fclose(fp);
    
    sprintf(fnt, "%s.csv", fnt);
    if((fp = fopen(fnt, "w")) == NULL ) {
        printf("\nCannot open the file! \n");
        exit(EXIT_FAILURE);
    }
    for (y1 = 0 ; y1 < ROW ; y1++) {
        for (x1 = 0 ; x1 < 162 * COL ; x1++) {
            fprintf(fp, "%f,", H[y1][x1]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

double fwinpat(int g_ang1[ROW][COL], int g_ang2[ROW][COL], double D[ROW][COL * 8], double ndis[(2 * ROW - 1) * (2 * COL - 1)], int coor[(2 * ROW - 1) * (2 * COL - 1)][2]) {
    /* calculation of mean of nearest-neighbor interpoint distances */
    /* with the same angle code between two images */
    double min, minInit, delta, dnn1, dnn2;
    int x1, y1, x2, y2;
    int angcode;
    int count1, count2;
    
    minInit = sqrt((ROW - 2 * MARGINE) * (ROW - 2 * MARGINE) + (COL - 2 * MARGINE) * (COL - 2 * MARGINE));
    /* from the 1st image */
    count1 = 0;
    dnn1 = 0.0;
    for (y1 = MARGINE ; y1 < ROW - MARGINE; y1++) {
        for (x1 = MARGINE ; x1 < COL - MARGINE ; x1++) {
            angcode = g_ang1[y1][x1];
            if (angcode == -1) continue;
            count1++;
            min = minInit;
            /*
             for (y2 = MARGINE ; y2 < ROW - MARGINE ; y2++) {
             for (x2 = MARGINE ; x2 < COL - MARGINE ; x2++) {
             if (g_ang2[y2][x2] != angcode) continue;
             delta = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
             if (delta < min) min = delta;
             }
             }
             */
            delta = D[y1][x1 + COL * angcode];
            if (delta < min) min = delta;
            // printf("angCode = %d, (%d, %d) nn1 = %f \n", angcode, x1, y1, min);
            dnn1 += min;
        }
    }
    dnn1 /= (double)count1;
    // printf("  count1  %d  ", count1);
    
    /* from the 2nd image */
    count2 = 0;
    dnn2 = 0.0;
    for (y2 = MARGINE ; y2 < ROW - MARGINE ; y2++) {
        for (x2 = MARGINE ; x2 < COL - MARGINE ; x2++) {
            angcode = g_ang2[y2][x2];
            if (angcode == -1) continue;
            count2++;
            min = minInit;
            /*
             for (y1 = MARGINE ; y1 < ROW - MARGINE ; y1++) {
             for (x1 = MARGINE ; x1 < COL - MARGINE ; x1++) {
             if (g_ang1[y1][x1] != angcode) continue;
             delta = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
             if (delta < min) min = delta;
             }
             }
             */
            for (y1 = 0 ; y1 < (2 * ROW - 1) * (2 * COL - 1) ; y1++) {
                if (y2 + coor[y1][0] < 0 || y2 + coor[y1][0] >= ROW || x2 + coor[y1][1] < 0 || x2 + coor[y1][1] >= COL ) continue;
                if (g_ang1[y2 + coor[y1][0]][x2 + coor[y1][1]] != angcode) continue;
                delta = ndis[y1];
                // printf("y1 = %d nn1 = %f \n", y1, ndis[y1]);
                if (delta < min) min = delta;
                break;
            }
            dnn2 += min;
        }
    }
    dnn2 /= (double)count2;
    // printf("  count2  %d  ", count2);
    
    /* printf("Gauss parameter %f  %f  \n", dnn1, dnn2); */
    return (dnn1 + dnn2)/2.0;
}





