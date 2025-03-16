#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>
#include "functions.h"
using namespace std;
#define PI 3.14159265
// function to initialize and output the starting position.
//double** initialization(double** r, int npart, double xlength, double ylength, double dl, double ds, double bl) {
double** initialization(double** r, int npart, double xlength, double ylength, double dl, double ds, double bl, int coin, int num_cent, double dc) {

    int i;
    int j;
    double x, y, w, v;
    double a, b, c;
    double dis1, dis2, dis3, dis4, discc; // those are the 'distance' in 2D plane, should not be used at the force calculation, but ok to use at initialization


    i = 0;
    //srand((unsigned)time(NULL));// seed, if not set, then can fix rand
    int t = std::chrono::steady_clock::now().time_since_epoch().count();

    srand(t); // try more accurate seed in ms
    x = rand() / double(RAND_MAX) * xlength / 10 + xlength / 2;
    y = rand() / double(RAND_MAX) * ylength / 10 + ylength / 2;

    r[i][0] = x;
    r[i][1] = y;
    r[i][2] = x;
    r[i][3] = y;

    for (i = 1; i < num_cent; i++)
    {
    begin_c:
        x = rand() / double(RAND_MAX) * xlength;
        y = rand() / double(RAND_MAX) * ylength;
        r[i][0] = x;
        r[i][1] = y;
        r[i][2] = x;
        r[i][3] = y;
        for (j = 0; j < i; j++)
        {
            discc = pow(pow(r[i][0] - r[j][0] - xlength * int((r[i][0] - r[j][0]) / (xlength / 2)), 2) + pow(r[i][1] - r[j][1] - ylength * int((r[i][1] - r[j][1]) / (ylength / 2)), 2), 0.5);
            if (discc < 6 * dc) {
                goto begin_c;
            }
        }
    }

    for (i = num_cent; i < npart; i++)
    {
    begin:
            x = rand() / double(RAND_MAX) * xlength;
            y = rand() / double(RAND_MAX) * ylength;
            a = rand() / double(RAND_MAX) - 0.5;
            b = rand() / double(RAND_MAX) - 0.5;
            c = pow((bl * bl - pow(dl - ds, 2) / 4) / (a * a + b * b), 0.5);
        w = x + a * c;
        v = y + b * c;
        r[i][0] = x;
        r[i][1] = y;
        r[i][2] = w;
        r[i][3] = v;
          for (j = 0; j < num_cent; j++)
        {
            dis1 = pow(pow(r[i][0] - r[j][0] - xlength * int((r[i][0] - r[j][0]) / (xlength / 2)), 2) + pow(r[i][1] - r[j][1] - ylength * int((r[i][1] - r[j][1]) / (ylength / 2)), 2), 0.5);
            dis4 = pow(pow(r[i][2] - r[j][0] - xlength * int((r[i][2] - r[j][0]) / (xlength / 2)), 2) + pow(r[i][3] - r[j][1] - ylength * int((r[i][3] - r[j][1]) / (ylength / 2)), 2), 0.5);

            if (dis1 < (dl + dc) / 2 || dis4 < (ds + dc) / 2) {
                goto begin;
            }
        }
        for (j = num_cent; j < i; j++)
        {
            dis1 = pow(pow(r[i][0] - r[j][0] - xlength * int((r[i][0] - r[j][0]) / (xlength / 2)), 2) + pow(r[i][1] - r[j][1] - ylength * int((r[i][1] - r[j][1]) / (ylength / 2)), 2), 0.5);
            dis2 = pow(pow(r[i][0] - r[j][2] - xlength * int((r[i][0] - r[j][2]) / (xlength / 2)), 2) + pow(r[i][1] - r[j][3] - ylength * int((r[i][1] - r[j][3]) / (ylength / 2)), 2), 0.5);
            dis3 = pow(pow(r[i][2] - r[j][2] - xlength * int((r[i][2] - r[j][2]) / (xlength / 2)), 2) + pow(r[i][3] - r[j][3] - ylength * int((r[i][3] - r[j][3]) / (ylength / 2)), 2), 0.5);
            dis4 = pow(pow(r[i][2] - r[j][0] - xlength * int((r[i][2] - r[j][0]) / (xlength / 2)), 2) + pow(r[i][3] - r[j][1] - ylength * int((r[i][3] - r[j][1]) / (ylength / 2)), 2), 0.5);

            if (dis1 < dl || dis3 < ds || dis2 < (dl + ds) / 2 || dis4 < (dl + ds) / 2) {
                goto begin;
            }

            if (w < x) {
                r[i][4] = atan((v - y) / (w - x)) + PI;
            }
            else
            {
                r[i][4] = atan((v - y) / (w - x));
            }
            r[i][5] = r[i][4] + PI;
        }
    }
    ////1110beta
    //for (i = 0; i < num_cent; i++) {
    //    r[i][6] = ds / 2  + 1000;
    //    r[i][7] = ds / 2 + ds / 3 * 2.44949;
    //}
    //for (i = num_cent; i < npart; i++) {
    //    r[i][6] = dl / 2;
    //    r[i][7] = ds / 2;
    //}
 
    //function to read unfinished particle position
    //int i;
    //std::ifstream fp;
    //if (coin == 0) {
    //    fp.open("./particle position0_unfinished.txt", std::ios::out); // remember to change filename here if needed
    //}
    //if (coin == 1) {
    //    fp.open("./particle position1_unfinished.txt", std::ios::out); // remember to change filename here if needed
    //}
    //if (coin == 2) {
    //    fp.open("./particle position2_unfinished.txt", std::ios::out); // remember to change filename here if needed
    //}
    //// check whether can open
    //if (!fp.is_open()) {
    //    std::cout << "fail to open the unfinished file\n";
    //}
    //else
    //    fp.seekg(-1, ios_base::end);

    //for (i = 0; i < 3; i++) //read the third last line, to avoid corruption in the last cycle
    //{while (fp.peek()!='\n')
    //    {
    //    fp.seekg(-1, fp.cur);
    //    }
    //fp.seekg(-1, fp.cur);
    //}
    //fp.seekg(3, fp.cur); //from the last byte at the end to the first byte at the beginning, note '\n' is 2 bytes in UFT-8

    //int indext;
    //fp >> indext;
    //indext = indext - 100; //careful, this one is related to the frequency of data storage in main.cpp [if (t % 100 == 0)] 
    //int indext1;
    //finding_last_frame:
    //while (fp.peek() != '\n')
    //{
    //    fp.seekg(-1, fp.cur);
    //}
    //fp.seekg(-1, fp.cur);
    //while (fp.peek() != '\n')
    //{
    //    fp.seekg(-1, fp.cur);
    //}
    //fp.seekg(-1, fp.cur);
    //fp.seekg(3, fp.cur);
    //fp >> indext1;
    //if (indext1 != indext)
    //    goto finding_last_frame;

    //for (i = 0; i < npart; i++) //go to the beginning of the second last frame 
    //{
    //    while (fp.peek() != '\n')
    //    {
    //        fp.seekg(-1, fp.cur);
    //    }
    //    fp.seekg(-1, fp.cur);
    //}
    //fp.seekg(3, fp.cur);

    //for (i = 0; i < npart; i++)
    //    fp >> indext >> r[i][0] >> r[i][1] >> r[i][2] >> r[i][3] >> r[i][4] >> r[i][5] >> r[i][6] >> r[i][7];




    //fp.close();




    return r;
}

