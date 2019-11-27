#ifndef ACCLGPT_H
#define ACCLGPT_H

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

#endif
