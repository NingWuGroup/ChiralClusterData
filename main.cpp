#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <ctime>
#include <cstdlib>
//#include "global_variables.h"
#include "functions.h"
#include <thread>
#include <vector>
#include <Windows.h>
using namespace std;

void dosomething(int cycle) {
    double t1 = clock();
    double config[29];
    int npart;
    double dl;
    double ds;
    double Tr;
    double lambda;
    double Ex;
    double Ey;
    double Ez;
    double Dt;
    double runtime;
    double bl;
    double Al;
    double As;
    double CMl;
    double CMs;
    double xlength;
    double ylength;
    double k1;
    double k2;
    double CMMl_i;
    double CMMs_i;
    double H;
    double lambdaMT;
    double lambdaMF;
    double omega;
    double CMMl_r;
    double CMMs_r;
    int num_cent;
    double dc;
    std::ifstream fp;
    fp.open("./config.txt", std::ios::in);
    // check whether can open config.txt
    if (!fp.is_open()) {
        std::cout << "fail to open file\n";
    }
    short loop = 0;
    while (!fp.eof())
    {
        fp >> config[loop];
        loop++;
    }
    fp.close();

    npart = config[0];
    dl = config[1];
    ds = config[2];
    Tr = config[3];
    lambda = config[4];
    Ex = config[5];
    Ey = config[6];
    Ez = config[7];
    Dt = config[8];
    runtime = config[9];
    bl = config[10];
    Al = config[11];
    As = config[12];
    CMl = config[13];
    CMs = config[14];
    xlength = config[15];
    ylength = config[16];
    k1 = config[17];
    k2 = config[18];
    CMMl_i = config[19];
    CMMs_i = config[20];
    H = config[21];
    //lambdaMT = config[22];
    lambdaMF = -config[23];
    omega = -config[24];
    CMMl_r = config[25];
    CMMs_r = config[26];
    num_cent = config[27];
    dc = config[28];
    #define PI 3.14159265
    // c++ can only define arrays with constant size, e.g. r[111][222] instead of r[m][n]
    // allocate space for position array, each row is the x and y coordinates and the orientation (angle from x axis, classic definition) for large and small lobes in a dimer (xl, yl, xs, ys, phil, phis)
    double** p = new double* [npart];
    for (int i = 0; i < npart; i++) {
        p[i] = new double[8];
    }
    for (int i = 0; i < npart; i++)
        for (int j = 0; j < 8; j++)
            p[i][j] = 0;
    double** f = new double* [npart];
    for (int i = 0; i < npart; i++) {
        f[i] = new double[8];
    }
    double** d = new double* [npart];
    for (int i = 0; i < npart; i++) {
        d[i] = new double[8];
    }
    double* dx = new double[npart];
    double* dy = new double[npart];
    double* dist = new double[npart];
    double* dx1 = new double[npart * (npart - 1) / 2];
    double* dy1 = new double[npart * (npart - 1) / 2];
    double* dx2 = new double[npart * (npart - 1) / 2];
    double* dy2 = new double[npart * (npart - 1) / 2];
    double* dx3 = new double[npart * (npart - 1) / 2];
    double* dy3 = new double[npart * (npart - 1) / 2];
    double* dx4 = new double[npart * (npart - 1) / 2];
    double* dy4 = new double[npart * (npart - 1) / 2];
    double* dz1 = new double[npart * (npart - 1) / 2];
    double* dz2 = new double[npart * (npart - 1) / 2];
    double* dz3 = new double[npart * (npart - 1) / 2];
    double* dz4 = new double[npart * (npart - 1) / 2];
    double* dist1 = new double[npart * (npart - 1) / 2];
    double* dist2 = new double[npart * (npart - 1) / 2];
    double* dist3 = new double[npart * (npart - 1) / 2];
    double* dist4 = new double[npart * (npart - 1) / 2];
    int* itemp = new int[npart * (npart - 1) / 2];
    int* jtemp = new int[npart * (npart - 1) / 2];
    double* zs = new double[npart];
    double* zl = new double[npart];



    lambdaMF = 0;
    //lambdaMF = 118 - 6 * cycle;
    lambdaMT = lambdaMF * 24 * CMMl_i / (pow(CMMl_i, 2) + pow(CMMl_r, 2));
    int coin = 0;
    //if (cycle == 0)
    //    coin = 0;
    //if (cycle == 1)
    //    coin = 1;
    //if (cycle == 2)
    //    coin = 2;
    p = initialization(p, npart, xlength, ylength, dl, ds, bl, coin,num_cent,dc);
    cout << "initialization complete" << endl;

    ofstream fo;
    ofstream debug;

    fo.open("./particle position" + to_string(cycle) + ".txt", std::ios::out);
    debug.open("./force" + to_string(cycle) + ".txt", std::ios::out);
    double theta = 0;
    for (int t = 0; t < runtime / Dt; t++) {
      //for (int t = 0; t < 1; t++) {
        for (int i = 0; i < npart; i++)
            for (int j = 0; j < 8; j++)
            {
                f[i][j] = 0;
                d[i][j] = 0;
            }
        //if (t * Dt > 0.2) {
        //    lambda = 1500;
        //}
        //if (t * Dt > 0.4) {
        //    lambda = 1400;
        //}
        //if (t * Dt > 0.6) {
        //    lambda = 1300;
        //}
        //if (t * Dt > 0.8) {
        //    lambda = 1200;
        //}
        //if (t * Dt > 1) {
        //    lambda = 1100;
        //}
        //if (t * Dt > 1.2) {
        //    lambda = 1000;
        //}
        //if (t * Dt > 1.4) {
        //    lambda = 900;
        //}
        //if (t * Dt > 1.6) {
        //    lambda = 800;
        //}
        //if (t * Dt > 1.8) {
        //    lambda = 700;
        //}
        //if (t * Dt > 2) {
        //    lambda = 600;
        //}
        //if (t * Dt > 2.2) {
        //    lambda = 500;
        //}
        //if (t * Dt > 2.4) {
        //    lambda = 400;
        //}
        Al = config[11];
        As = config[12];

        // function to fix the orientation
            //for (int i = 0; i < npart; i++) {
            //    if (p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2)) > 0) {
            //        p[i][4] = atan((p[i][3] - p[i][1] - ylength * int((p[i][3] - p[i][1]) / (ylength / 2))) / (p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2))));
            //            p[i][5] = atan((p[i][3] - p[i][1] - ylength * int((p[i][3] - p[i][1]) / (ylength / 2))) / (p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2)))) + PI;
            //    }
            //    else {
            //        p[i][4] =PI+ atan((p[i][3] - p[i][1] - ylength * int((p[i][3] - p[i][1]) / (ylength / 2))) / (p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2))));
            //        p[i][5] = PI+ atan((p[i][3] - p[i][1] - ylength * int((p[i][3] - p[i][1]) / (ylength / 2))) / (p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2)))) + PI;

            //    }
            //}

            //d = debugforce(d, p,  dx,  dy,  dist,  dx1,  dy1,  dx2,  dy2,  dx3,  dy3,  dx4,  dy4,  dist1,  dist2,  dist3,  dist4,  itemp,  jtemp, npart, xlength, ylength, dl, ds, bl, Tr, lambda, CMl, CMs, Al, As, k1, k2, Ex, Ey, Ez, Dt, CMMl_i, CMMs_i, CMMl_r, CMMs_r, H, lambdaMT, lambdaMF, theta);
           
            f = force_displacement(f, p,  dx,  dy,  dist,  dx1,  dy1,  dx2,  dy2,  dx3,  dy3,  dx4,  dy4, dz1,dz2,dz3,dz4,zl,zs,  dist1,  dist2,  dist3,  dist4,  itemp,  jtemp, npart, xlength, ylength, dl, ds, bl, Tr, lambda, CMl, CMs, Al, As, k1, k2, Ex, Ey, Ez, Dt, CMMl_i, CMMs_i, CMMl_r, CMMs_r, H, lambdaMT, lambdaMF, theta, num_cent, dc);
  
        theta = theta + 2 * PI * omega * Dt;
        for (int i = 0; i < npart; i++) {
            for (int j = 0; j < 8; j++)
                p[i][j] = p[i][j] + f[i][j];

            while (p[i][0] < 0)
                p[i][0] = p[i][0] + xlength;
            while (p[i][0] > xlength)
                p[i][0] = p[i][0] - xlength;
            while (p[i][2] < 0)
                p[i][2] = p[i][2] + xlength;
            while (p[i][2] > xlength)
                p[i][2] = p[i][2] - xlength;
            while (p[i][1] < 0)
                p[i][1] = p[i][1] + ylength;
            while (p[i][1] > ylength)
                p[i][1] = p[i][1] - ylength;
            while (p[i][3] < 0)
                p[i][3] = p[i][3] + ylength;
            while (p[i][3] > ylength)
                p[i][3] = p[i][3] - ylength;
        }
        if (t % 100 == 0) {
            for (int i = 0; i < npart; i++)
                fo << t + 1 << " " << p[i][0] << " " << p[i][1] << " " << p[i][2] << " " << p[i][3] << " " << p[i][4] << " " << p[i][5] << " " << p[i][6] << " " << p[i][7] << endl;
            fo.flush();
        }
        
        if (t % 100 == 0) {
            for (int i = 0; i < npart; i++)
                debug << t + 1 << " " << d[i][0] << " " << d[i][1] << " " << d[i][2] << " " << d[i][3] << endl;
            debug.flush();
        }
   }
    fo.close();
    debug.close();
    // release the space allocated for the position array
        for (int i = 0; i < npart; i++) {
            delete[] p[i];
            delete[] f[i];
            delete[] d[i];
        }
    delete[] p;
    delete[] f;
    delete[] d;
    delete[] dx;
    delete[] dy;
    delete[] dist;
    delete[] dx1;
    delete[] dy1;
    delete[] dz1;
    delete[] dx2;
    delete[] dy2;
    delete[] dz2;
    delete[] dx3;
    delete[] dy3;
    delete[] dz3;
    delete[] dx4;
    delete[] dy4;
    delete[] dz4;
    delete[] zl;
    delete[] zs;
    delete[] dist1;
    delete[] dist2;
    delete[] dist3;
    delete[] dist4;
    delete[] itemp;
    delete[] jtemp;
    double t2 = clock();
    double trun = t2 - t1;
    cout << "cycle" << cycle << " takes" << ":" << trun / CLOCKS_PER_SEC << " seconds" << endl;


}

void main()
{

    int Numthread = 1;
    std::vector<std::thread> Threads;
    for (int i = 0; i < Numthread; i++) {
        Threads.push_back(std::thread(dosomething,i));
    }

    for (int i = 0; i < Numthread; i++) {
        Threads[i].join();
    }
    Beep(523, 100);
    Beep(578, 100);
    Beep(659, 100);
    Beep(698, 100);
    Beep(784, 100);
    system("PAUSE");
}

