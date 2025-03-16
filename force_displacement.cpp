#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <random>
#define PI 3.14159265
using namespace std;
double** force_displacement(double** r, double** p, double* dx, double* dy, double* dist, double* dx1, double* dy1, double* dx2, double* dy2, double* dx3, double* dy3, double* dx4, double* dy4, double* dz1, double* dz2, double* dz3, double* dz4, double* zl, double* zs, double* dist1, double* dist2, double* dist3, double* dist4, int* itemp, int* jtemp, int npart, double xlength, double ylength, double dl, double ds, double bl, double Tr, double lambda, double CMl, double CMs, double Al, double As, double k1, double k2, double Ex, double Ey, double Ez, double Dt, double CMMl_i, double CMMs_i, double CMMl_r, double CMMs_r, double H, double lambdaMT, double lambdaMF, double theta,int num_cent, double dc) {
	double Ksil = dl / dl;
	double Ksis = ds / dl;
	double sigma_ls = (dl + ds) / 2;
	//double z_c_s = dc / 2;
	//double z_c_s = ds / 2 + ds / 3 * 2.44949;
	//double z_c_s = ds / 2 + ds / 2 * 1.4;
	//double z_c_l = z_c_s + bl; // try to hide one of the lobe in the center, since this code was originally a standing dimer
	//double z_c_l = 1000;
	//double z_s = ds / 2;
	//double z_l = dl / 2;


	// define the distance variables between two lobes in the same dimer, big first, small second (see main.cpp)
	for (int i = 0; i < npart; i++) {
		dx[i] = p[i][2] - p[i][0] - xlength * int((p[i][2] - p[i][0]) / (xlength / 2));
		dy[i] = p[i][3] - p[i][1] - ylength * int((p[i][3] - p[i][1]) / (ylength / 2));
		dist[i] = pow(pow(dx[i], 2) + pow(dy[i], 2) + (dl - ds) * (dl - ds) / 4, 0.5);
		//if (abs(dist[i] - bl) > 0.1)
		//	cout << "alert" << abs(dist[i] - bl) << '\n';
		//else
		//	cout << abs(dist[i] - bl)<<'\n';
	}
	//for (int i = 0; i < num_cent; i++) {
	//	zs[i] = z_c_s;
	//	zl[i] = z_c_l;
	//}
	//for (int i = num_cent; i < npart; i++) {
	//	zs[i] = z_s;
	//	zl[i] = z_l;
	//}
	for (int i = 0; i < npart; i++) {
		zs[i] = p[i][7];
		zl[i] = p[i][6];
	}


	// define the distance variables between each two lobes in different dimers
	int h = 0; // one counter used later
	for (int i = 1; i < npart; i++) {
		for (int j = 0; j < i; j++) {
			dx1[h] = p[i][0] - p[j][0] - xlength * int((p[i][0] - p[j][0]) / (xlength / 2));
			dx2[h] = p[i][0] - p[j][2] - xlength * int((p[i][0] - p[j][2]) / (xlength / 2));
			dx3[h] = p[i][2] - p[j][2] - xlength * int((p[i][2] - p[j][2]) / (xlength / 2));
			dx4[h] = p[i][2] - p[j][0] - xlength * int((p[i][2] - p[j][0]) / (xlength / 2));
			dy1[h] = p[i][1] - p[j][1] - ylength * int((p[i][1] - p[j][1]) / (ylength / 2));
			dy2[h] = p[i][1] - p[j][3] - ylength * int((p[i][1] - p[j][3]) / (ylength / 2));
			dy3[h] = p[i][3] - p[j][3] - ylength * int((p[i][3] - p[j][3]) / (ylength / 2));
			dy4[h] = p[i][3] - p[j][1] - ylength * int((p[i][3] - p[j][1]) / (ylength / 2));
			dz1[h] = zl[i] - zl[j];
			dz2[h] = zl[i] - zs[j];
			dz3[h] = zs[i] - zs[j];
			dz4[h] = zs[i] - zl[j];
			dist1[h] = pow(pow(dx1[h],2) + pow(dy1[h],2) + pow(dz1[h], 2), 0.5);
			dist2[h] = pow(pow(dx2[h],2) + pow(dy2[h],2) + pow(dz2[h], 2), 0.5);
			dist3[h] = pow(pow(dx3[h],2) + pow(dy3[h],2) + pow(dz3[h], 2), 0.5);
			dist4[h] = pow(pow(dx4[h],2) + pow(dy4[h],2) + pow(dz4[h], 2), 0.5);
			itemp[h] = i;
			jtemp[h] = j;
			h++;
		}
	}
	 //BD random force
	unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> SND(0.0, 1.0);
	for (int i = 0; i < npart; i++) {
		for (int j = 0; j < 2; j++) {
			double a = SND(generator);
			r[i][j] = r[i][j] + 1.0 * a / pow(Dt / Ksil, 0.5);
		}
		for (int j = 2; j < 4; j++) {
			double b = SND(generator);
			r[i][j] = r[i][j] + 1.0 * b / pow(Dt / Ksis, 0.5);
		}
	}

	 //EHD force classic

	for (int i = num_cent; i < npart; i++) {
		r[i][0] = r[i][0] - dl * pow(ds, 3) * dl / 2 * As * (dx[i]) / (pow(sigma_ls, 5) * pow(ds, 2) / 4);
		r[i][1] = r[i][1] - dl * pow(ds, 3) * dl / 2 * As * (dy[i]) / (pow(sigma_ls, 5) * pow(ds, 2) / 4);
		r[i][2] = r[i][2] - ds * pow(dl, 3) * ds / 2 * Al * (-dx[i]) / (pow(sigma_ls, 5) * pow(dl, 2) / 4);
		r[i][3] = r[i][3] - ds * pow(dl, 3) * ds / 2 * Al * (-dy[i]) / (pow(sigma_ls, 5) * pow(dl, 2) / 4);
	}

	for (int i = 0; i < npart * (npart - 1) / 2; i++) {

			r[itemp[i]][2] = r[itemp[i]][2] - ds * pow(ds, 3) * ds / 2 * As * (-dx3[i]) / (pow(dist3[i], 5) * pow(ds, 2) / 4) - ds * pow(dl, 3) * ds / 2 * Al * (-dx4[i]) / (pow(dist4[i], 5) * pow(dl, 2) / 4);
			r[jtemp[i]][2] = r[jtemp[i]][2] - ds * pow(ds, 3) * ds / 2 * As * (dx3[i]) / (pow(dist3[i], 5) * pow(ds, 2) / 4) - ds * pow(dl, 3) * ds / 2 * Al * (dx2[i]) / (pow(dist2[i], 5) * pow(dl, 2) / 4);
			r[itemp[i]][3] = r[itemp[i]][3] - ds * pow(ds, 3) * ds / 2 * As * (-dy3[i]) / (pow(dist3[i], 5) * pow(ds, 2) / 4) - ds * pow(dl, 3) * ds / 2 * Al * (-dy4[i]) / (pow(dist4[i], 5) * pow(dl, 2) / 4);
			r[jtemp[i]][3] = r[jtemp[i]][3] - ds * pow(ds, 3) * ds / 2 * As * (dy3[i]) / (pow(dist3[i], 5) * pow(ds, 2) / 4) - ds * pow(dl, 3) * ds / 2 * Al * (dy2[i]) / (pow(dist2[i], 5) * pow(dl, 2) / 4);
			r[itemp[i]][0] = r[itemp[i]][0] -  dl * pow(ds, 3) * dl / 2 * As * (-dx2[i]) / (pow(dist2[i], 5) * pow(ds, 2) / 4) - dl * pow(dl, 3) * dl / 2 * Al * (-dx1[i]) / (pow(dist1[i], 5) * pow(dl, 2) / 4);
			r[jtemp[i]][0] = r[jtemp[i]][0] -  dl * pow(ds, 3) * dl / 2 * As * (dx4[i]) / (pow(dist4[i], 5) * pow(ds, 2) / 4) - dl * pow(dl, 3) * dl / 2 * Al * (dx1[i]) / (pow(dist1[i], 5) * pow(dl, 2) / 4);
			r[itemp[i]][1] = r[itemp[i]][1] -  dl * pow(ds, 3) * dl / 2 * As * (-dy2[i]) / (pow(dist2[i], 5) * pow(ds, 2) / 4) - dl * pow(dl, 3) * dl / 2 * Al * (-dy1[i]) / (pow(dist1[i], 5) * pow(dl, 2) / 4);
			r[jtemp[i]][1] = r[jtemp[i]][1] -  dl * pow(ds, 3) * dl / 2 * As * (dy4[i]) / (pow(dist4[i], 5) * pow(ds, 2) / 4) - dl * pow(dl, 3) * dl / 2 * Al * (dy1[i]) / (pow(dist1[i], 5) * pow(dl, 2) / 4);

	}




	// DIP force (note that this part only contributes to the force, not the torque)
	double dipole_l_x = Ex * pow(dl, 3);
	double dipole_s_x = Ex * pow(ds, 3) * CMs / CMl;
	double dipole_l_y = Ey * pow(dl, 3);
	double dipole_s_y = Ey * pow(ds, 3) * CMs / CMl;
	double dipole_l_z = Ez * pow(dl, 3);
	double dipole_s_z = Ez * pow(ds, 3) * CMs / CMl;
	double dot_prdct_ir_ls = dipole_l_z * (ds - dl) / 2; // i is large j is small
	double dot_prdct_jr_ls = dipole_s_z * (ds - dl) / 2;
	//double dot_prdct_ir_sl = -dot_prdct_jr_ls; // i is small j is large
	//double dot_prdct_jr_sl = -dot_prdct_ir_ls;
	//double dot_prdct_ir_ll = 0; // both are large
	//double dot_prdct_jr_ll = 0;
	//double dot_prdct_ir_ss = dipole_s_z; // both are small
	//double dot_prdct_jr_ss = dipole_s_z;
	//// above products are only in z-direction
	double dot_prdct_ij_ls = dipole_l_z * dipole_s_z + dipole_l_y * dipole_s_y + dipole_l_x * dipole_s_x;
	double dot_prdct_ij_sl = dot_prdct_ij_ls;
	double dot_prdct_ij_ll = pow(dipole_l_z, 2) + pow(dipole_l_y, 2) + pow(dipole_l_x, 2);
	double dot_prdct_ij_ss = pow(dipole_s_z, 2) + pow(dipole_s_y, 2) + pow(dipole_s_x, 2);


	//for (int i = 0; i < npart * (npart - 1) / 2; i++) {
	//	r[itemp[i]][0] = r[itemp[i]][0] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (-dx1[i])) - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (-dx2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dx2[i]));
	//	r[jtemp[i]][0] = r[jtemp[i]][0] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (dx1[i])) - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (dx4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dx4[i]));
	//	r[itemp[i]][1] = r[itemp[i]][1] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (-dy1[i])) - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (-dy2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dy2[i]));
	//	r[jtemp[i]][1] = r[jtemp[i]][1] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (dy1[i])) - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (dy4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dy4[i]));
	//	r[itemp[i]][2] = r[itemp[i]][2] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dx3[i])) - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (-dx4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dx4[i]));
	//	r[jtemp[i]][2] = r[jtemp[i]][2] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (dx3[i])) - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (dx2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dx2[i]));
	//	r[itemp[i]][3] = r[itemp[i]][3] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dy3[i])) - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (-dy4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dy4[i]));
	//	r[jtemp[i]][3] = r[jtemp[i]][3] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (dy3[i])) - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (dy2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dy2[i]));
	//}
	
	for (int i = 0; i < npart * (npart - 1) / 2; i++) {
		if (jtemp[i] > num_cent-1 && itemp[i] > num_cent - 1) {
			r[itemp[i]][0] = r[itemp[i]][0] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (-dx1[i]));
			r[jtemp[i]][0] = r[jtemp[i]][0] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (dx1[i]));
			r[itemp[i]][1] = r[itemp[i]][1] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (-dy1[i]));
			r[jtemp[i]][1] = r[jtemp[i]][1] - lambda / pow(dist1[i], 5) * (dot_prdct_ij_ll * (dy1[i]));

			r[itemp[i]][0] = r[itemp[i]][0] - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (-dx2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dx2[i]));
			r[itemp[i]][1] = r[itemp[i]][1] - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (-dy2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dy2[i]));
			r[jtemp[i]][2] = r[jtemp[i]][2] - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (dx2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dx2[i]));
			r[jtemp[i]][3] = r[jtemp[i]][3] - lambda / pow(dist2[i], 5) * (dot_prdct_ij_ls * (dy2[i]) - 5 / pow(dist2[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dy2[i]));

			r[itemp[i]][2] = r[itemp[i]][2] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dx3[i]));
			r[jtemp[i]][2] = r[jtemp[i]][2] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (dx3[i]));
			r[itemp[i]][3] = r[itemp[i]][3] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dy3[i]));
			r[jtemp[i]][3] = r[jtemp[i]][3] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (dy3[i]));

			r[jtemp[i]][0] = r[jtemp[i]][0] - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (dx4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dx4[i]));
			r[jtemp[i]][1] = r[jtemp[i]][1] - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (dy4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (dy4[i]));
			r[itemp[i]][2] = r[itemp[i]][2] - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (-dx4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dx4[i]));
			r[itemp[i]][3] = r[itemp[i]][3] - lambda / pow(dist4[i], 5) * (dot_prdct_ij_ls * (-dy4[i]) - 5 / pow(dist4[i], 2) * dot_prdct_ir_ls * dot_prdct_jr_ls * (-dy4[i]));
		}
		if (jtemp[i] < num_cent && itemp[i] > num_cent - 1) {
			double dot_prdct_ir_ss = dipole_s_z * (-dz3[i]); // both are small
			double dot_prdct_jr_ss = dipole_s_z * (-dz3[i]);
			r[itemp[i]][2] = r[itemp[i]][2] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dx3[i]) - 5 / pow(dist3[i], 2) * dot_prdct_ir_ss * dot_prdct_jr_ss * (-dx3[i]));
			r[itemp[i]][3] = r[itemp[i]][3] - lambda / pow(dist3[i], 5) * (dot_prdct_ij_ss * (-dy3[i]) - 5 / pow(dist3[i], 2) * dot_prdct_ir_ss * dot_prdct_jr_ss * (-dy3[i]));
		}
	}
	

	// magnetic DIP force assuming all particles have same dipole directions
	//double Hx = H * cos(theta);
	//double Hy = H * sin(theta);
	double p_l_x = H * pow(dl, 3) * (CMMl_r * cos(theta) - CMMl_i * sin(theta)); 
	double p_s_x = H * pow(ds, 3) * (CMMs_r * cos(theta) - CMMs_i * sin(theta));
	double p_l_y = H * pow(dl, 3) * (CMMl_r * sin(theta) + CMMl_i * cos(theta));
	double p_s_y = H * pow(ds, 3) * (CMMs_r * sin(theta) + CMMs_i * cos(theta));
	//double d_p_ir_ls = 0; // i is large j is small
	//double d_p_jr_ls = 0;
	//double d_p_ir_sl = 0; // i is small j is large
	//double d_p_jr_sl = 0;
	//double d_p_ir_ll = 0; // both are large
	//double d_p_jr_ll = 0;
	//double d_p_ir_ss = 0; // both are small
	//double d_p_jr_ss = 0;
	double d_p_ij_ls = p_l_y * p_s_y + p_l_x * p_s_x;
	double d_p_ij_sl = d_p_ij_ls;
	double d_p_ij_ll = pow(p_l_y, 2) + pow(p_l_x, 2);
	double d_p_ij_ss = pow(p_s_y, 2) + pow(p_s_x, 2);
	for (int i = 0; i < npart * (npart - 1) / 2; i++) {
		r[itemp[i]][0] = r[itemp[i]][0] - lambdaMF / pow(dist1[i], 5) * ((+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * p_l_x + (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * p_l_x + d_p_ij_ll * (-dx1[i]) - 5 * (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) / pow(dist1[i], 2) * (-dx1[i])) - lambdaMF / pow(dist2[i], 5) * ((+p_l_x * (-dx2[i]) + p_l_y * (-dy2[i])) * p_s_x + (+p_s_x * (-dx2[i]) + p_s_y * (-dy2[i])) * p_l_x + d_p_ij_ls * (-dx2[i]) - 5 * (+p_l_x * (-dx2[i]) + p_l_y * (-dy2[i])) * (+p_s_x * (-dx2[i]) + p_s_y * (-dy2[i])) / pow(dist2[i], 2) * (-dx2[i]));
		r[jtemp[i]][0] = r[jtemp[i]][0] - lambdaMF / pow(dist1[i], 5) * ((+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * p_l_x + (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * p_l_x + d_p_ij_ll * (dx1[i]) - 5 * (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) / pow(dist1[i], 2) * (dx1[i])) - lambdaMF / pow(dist4[i], 5) * ((+p_s_x * (dx4[i]) + p_s_y * (dy4[i])) * p_l_x + (+p_l_x * (dx4[i]) + p_l_y * (dy4[i])) * p_s_x + d_p_ij_ls * (dx4[i]) - 5 * (+p_s_x * (dx4[i]) + p_s_y * (dy4[i])) * (+p_l_x * (dx4[i]) + p_l_y * (dy4[i])) / pow(dist4[i], 2) * (dx4[i]));
		r[itemp[i]][1] = r[itemp[i]][1] - lambdaMF / pow(dist1[i], 5) * ((+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * p_l_y + (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * p_l_y + d_p_ij_ll * (-dy1[i]) - 5 * (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) * (+p_l_x * (-dx1[i]) + p_l_y * (-dy1[i])) / pow(dist1[i], 2) * (-dy1[i])) - lambdaMF / pow(dist2[i], 5) * ((+p_l_x * (-dx2[i]) + p_l_y * (-dy2[i])) * p_s_y + (+p_s_x * (-dx2[i]) + p_s_y * (-dy2[i])) * p_l_y + d_p_ij_ls * (-dy2[i]) - 5 * (+p_l_x * (-dx2[i]) + p_l_y * (-dy2[i])) * (+p_s_x * (-dx2[i]) + p_s_y * (-dy2[i])) / pow(dist2[i], 2) * (-dy2[i]));
		r[jtemp[i]][1] = r[jtemp[i]][1] - lambdaMF / pow(dist1[i], 5) * ((+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * p_l_y + (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * p_l_y + d_p_ij_ll * (dy1[i]) - 5 * (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) * (+p_l_x * (dx1[i]) + p_l_y * (dy1[i])) / pow(dist1[i], 2) * (dy1[i])) - lambdaMF / pow(dist4[i], 5) * ((+p_s_x * (dx4[i]) + p_s_y * (dy4[i])) * p_l_y + (+p_l_x * (dx4[i]) + p_l_y * (dy4[i])) * p_s_y + d_p_ij_ls * (dy4[i]) - 5 * (+p_s_x * (dx4[i]) + p_s_y * (dy4[i])) * (+p_l_x * (dx4[i]) + p_l_y * (dy4[i])) / pow(dist4[i], 2) * (dy4[i]));
		r[itemp[i]][2] = r[itemp[i]][2] - lambdaMF / pow(dist3[i], 5) * ((+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * p_s_x + (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * p_s_x + d_p_ij_ss * (-dx3[i]) - 5 * (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) / pow(dist3[i], 2) * (-dx3[i])) - lambdaMF / pow(dist4[i], 5) * ((+p_s_x * (-dx4[i]) + p_s_y * (-dy4[i])) * p_l_x + (+p_l_x * (-dx4[i]) + p_l_y * (-dy4[i])) * p_s_x + d_p_ij_sl * (-dx4[i]) - 5 * (+p_s_x * (-dx4[i]) + p_s_y * (-dy4[i])) * (+p_l_x * (-dx4[i]) + p_l_y * (-dy4[i])) / pow(dist4[i], 2) * (-dx4[i]));
		r[jtemp[i]][2] = r[jtemp[i]][2] - lambdaMF / pow(dist3[i], 5) * ((+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * p_s_x + (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * p_s_x + d_p_ij_ss * (dx3[i]) - 5 * (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) / pow(dist3[i], 2) * (dx3[i])) - lambdaMF / pow(dist2[i], 5) * ((+p_l_x * (dx2[i]) + p_l_y * (dy2[i])) * p_s_x + (+p_s_x * (dx2[i]) + p_s_y * (dy2[i])) * p_l_x + d_p_ij_ls * (dx2[i]) - 5 * (+p_l_x * (dx2[i]) + p_l_y * (dy2[i])) * (+p_s_x * (dx2[i]) + p_s_y * (dy2[i])) / pow(dist2[i], 2) * (dx2[i]));
		r[itemp[i]][3] = r[itemp[i]][3] - lambdaMF / pow(dist3[i], 5) * ((+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * p_s_y + (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * p_s_y + d_p_ij_ss * (-dy3[i]) - 5 * (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) * (+p_s_x * (-dx3[i]) + p_s_y * (-dy3[i])) / pow(dist3[i], 2) * (-dy3[i])) - lambdaMF / pow(dist4[i], 5) * ((+p_s_x * (-dx4[i]) + p_s_y * (-dy4[i])) * p_l_y + (+p_l_x * (-dx4[i]) + p_l_y * (-dy4[i])) * p_s_y + d_p_ij_sl * (-dy4[i]) - 5 * (+p_s_x * (-dx4[i]) + p_s_y * (-dy4[i])) * (+p_l_x * (-dx4[i]) + p_l_y * (-dy4[i])) / pow(dist4[i], 2) * (-dy4[i]));
		r[jtemp[i]][3] = r[jtemp[i]][3] - lambdaMF / pow(dist3[i], 5) * ((+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * p_s_y + (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * p_s_y + d_p_ij_ss * (dy3[i]) - 5 * (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) * (+p_s_x * (dx3[i]) + p_s_y * (dy3[i])) / pow(dist3[i], 2) * (dy3[i])) - lambdaMF / pow(dist2[i], 5) * ((+p_l_x * (dx2[i]) + p_l_y * (dy2[i])) * p_s_y + (+p_s_x * (dx2[i]) + p_s_y * (dy2[i])) * p_l_y + d_p_ij_ls * (dy2[i]) - 5 * (+p_l_x * (dx2[i]) + p_l_y * (dy2[i])) * (+p_s_x * (dx2[i]) + p_s_y * (dy2[i])) / pow(dist2[i], 2) * (dy2[i]));
	}
	for (int i = 0; i < npart; i++) {
		r[i][0] = r[i][0] - lambdaMF / pow(sigma_ls, 5) * ((+p_l_x * (dx[i]) + p_l_y * (dy[i])) * p_s_x + (+p_s_x * (dx[i]) + p_s_y * (dy[i])) * p_l_x + d_p_ij_ls * (dx[i]) - 5 * (+p_l_x * (dx[i]) + p_l_y * (dy[i])) * (+p_s_x * (dx[i]) + p_s_y * (dy[i])) / pow(sigma_ls, 2) * (dx[i]));
		r[i][1] = r[i][1] - lambdaMF / pow(sigma_ls, 5) * ((+p_l_x * (dx[i]) + p_l_y * (dy[i])) * p_s_y + (+p_s_x * (dx[i]) + p_s_y * (dy[i])) * p_l_y + d_p_ij_ls * (dy[i]) - 5 * (+p_l_x * (dx[i]) + p_l_y * (dy[i])) * (+p_s_x * (dx[i]) + p_s_y * (dy[i])) / pow(sigma_ls, 2) * (dy[i]));
		r[i][2] = r[i][2] - lambdaMF / pow(sigma_ls, 5) * ((+p_s_x * (-dx[i]) + p_s_y * (-dy[i])) * p_l_x + (+p_l_x * (-dx[i]) + p_l_y * (-dy[i])) * p_s_x + d_p_ij_ls * (-dx[i]) - 5 * (+p_s_x * (-dx[i]) + p_s_y * (-dy[i])) * (+p_l_x * (-dx[i]) + p_l_y * (-dy[i])) / pow(sigma_ls, 2) * (-dx[i]));
		r[i][3] = r[i][3] - lambdaMF / pow(sigma_ls, 5) * ((+p_s_x * (-dx[i]) + p_s_y * (-dy[i])) * p_l_y + (+p_l_x * (-dx[i]) + p_l_y * (-dy[i])) * p_s_y + d_p_ij_ls * (-dy[i]) - 5 * (+p_s_x * (-dx[i]) + p_s_y * (-dy[i])) * (+p_l_x * (-dx[i]) + p_l_y * (-dy[i])) / pow(sigma_ls, 2) * (-dy[i]));
	}
	
	// LJ force
	for (int i = 0; i < npart * (npart - 1) / 2; i++) {	
		if (dist2[i] < sigma_ls *1.01) {
			r[itemp[i]][0] = r[itemp[i]][0] - 2 * Tr / pow(sigma_ls, 2) * (-dx2[i]) * (12 * pow(1.01 * sigma_ls / dist2[i], 14) - 6 * pow(1.01 * sigma_ls / dist2[i], 8));
			r[itemp[i]][1] = r[itemp[i]][1] - 2 * Tr / pow(sigma_ls, 2) * (-dy2[i]) * (12 * pow(1.01 * sigma_ls / dist2[i], 14) - 6 * pow(1.01 * sigma_ls / dist2[i], 8));
			r[jtemp[i]][2] = r[jtemp[i]][2] - 2 * Tr / pow(sigma_ls, 2) * (dx2[i]) * (12 * pow(1.01 * sigma_ls / dist2[i], 14) - 6 * pow(1.01 * sigma_ls / dist2[i], 8));
			r[jtemp[i]][3] = r[jtemp[i]][3] - 2 * Tr / pow(sigma_ls, 2) * (dy2[i]) * (12 * pow(1.01 * sigma_ls / dist2[i], 14) - 6 * pow(1.01 * sigma_ls / dist2[i], 8));
		}
		if (dist1[i] < dl * 1.01) {
			r[itemp[i]][0] = r[itemp[i]][0] - 2 * Tr / pow(dl, 2) * (-dx1[i]) * (12 * pow(1.01 * dl / dist1[i], 14) - 6 * pow(1.01 * dl / dist1[i], 8));
			r[itemp[i]][1] = r[itemp[i]][1] - 2 * Tr / pow(dl, 2) * (-dy1[i]) * (12 * pow(1.01 * dl / dist1[i], 14) - 6 * pow(1.01 * dl / dist1[i], 8));
			r[jtemp[i]][0] = r[jtemp[i]][0] - 2 * Tr / pow(dl, 2) * (dx1[i]) * (12 * pow(1.01 * dl / dist1[i], 14) - 6 * pow(1.01 * dl / dist1[i], 8));
			r[jtemp[i]][1] = r[jtemp[i]][1] - 2 * Tr / pow(dl, 2) * (dy1[i]) * (12 * pow(1.01 * dl / dist1[i], 14) - 6 * pow(1.01 * dl / dist1[i], 8));
		}
		if (dist4[i] < sigma_ls * 1.01) {
			r[itemp[i]][2] = r[itemp[i]][2] - 2 * Tr / pow(sigma_ls, 2) * (-dx4[i]) * (12 * pow(1.01 * sigma_ls / dist4[i], 14) - 6 * pow(1.01 * sigma_ls / dist4[i], 8));
			r[itemp[i]][3] = r[itemp[i]][3] - 2 * Tr / pow(sigma_ls, 2) * (-dy4[i]) * (12 * pow(1.01 * sigma_ls / dist4[i], 14) - 6 * pow(1.01 * sigma_ls / dist4[i], 8));
			r[jtemp[i]][0] = r[jtemp[i]][0] - 2 * Tr / pow(sigma_ls, 2) * (dx4[i]) * (12 * pow(1.01 * sigma_ls / dist4[i], 14) - 6 * pow(1.01 * sigma_ls / dist4[i], 8));
			r[jtemp[i]][1] = r[jtemp[i]][1] - 2 * Tr / pow(sigma_ls, 2) * (dy4[i]) * (12 * pow(1.01 * sigma_ls / dist4[i], 14) - 6 * pow(1.01 * sigma_ls / dist4[i], 8));
		}
		if (dist3[i] < ds * 1.0) {
			r[itemp[i]][2] = r[itemp[i]][2] - 2 * Tr / pow(ds, 2) * (-dx3[i]) * (12 * pow(1.0 * ds / dist3[i], 14) - 6 * pow(1.0 * ds / dist3[i], 8));
			r[itemp[i]][3] = r[itemp[i]][3] - 2 * Tr / pow(ds, 2) * (-dy3[i]) * (12 * pow(1.0 * ds / dist3[i], 14) - 6 * pow(1.0 * ds / dist3[i], 8));
			r[jtemp[i]][2] = r[jtemp[i]][2] - 2 * Tr / pow(ds, 2) * (dx3[i]) * (12 * pow(1.0 * ds / dist3[i], 14) - 6 * pow(1.0 * ds / dist3[i], 8));
			r[jtemp[i]][3] = r[jtemp[i]][3] - 2 * Tr / pow(ds, 2) * (dy3[i]) * (12 * pow(1.0 * ds / dist3[i], 14) - 6 * pow(1.0 * ds / dist3[i], 8));
		}
	}

	//// Try Yukawa instead of LJ force
	//double epslnY = 10;
	//double kappa = 1 / 0.01;
	//for (int i = 0; i < npart * (npart - 1) / 2; i++) {

	//	r[itemp[i]][0] = r[itemp[i]][0] - epslnY * exp(-kappa * sigma_ls * (dist2[i] / sigma_ls - 1)) / pow(dist2[i] / sigma_ls, 2) / sigma_ls  * (1 + kappa * dist2[i]) * (-dx2[i])/ dist2[i];
	//	r[itemp[i]][1] = r[itemp[i]][1] - epslnY * exp(-kappa * sigma_ls * (dist2[i] / sigma_ls - 1)) / pow(dist2[i] / sigma_ls, 2) / sigma_ls  * (1 + kappa * dist2[i]) * (-dy2[i]) / dist2[i];
	//	r[jtemp[i]][2] = r[jtemp[i]][2] - epslnY * exp(-kappa * sigma_ls * (dist2[i] / sigma_ls - 1)) / pow(dist2[i] / sigma_ls, 2) / sigma_ls  * (1 + kappa * dist2[i]) * (dx2[i]) / dist2[i];
	//	r[jtemp[i]][3] = r[jtemp[i]][3] - epslnY * exp(-kappa * sigma_ls * (dist2[i] / sigma_ls - 1)) / pow(dist2[i] / sigma_ls, 2) / sigma_ls  * (1 + kappa * dist2[i]) * (dy2[i]) / dist2[i];

	//		r[itemp[i]][0] = r[itemp[i]][0] - epslnY * exp(-kappa * dl * (dist1[i] / dl - 1)) / pow(dist1[i] / dl, 2) /dl  * (1 + kappa * dist1[i]) * (-dx1[i]) / dist1[i];
	//		r[itemp[i]][1] = r[itemp[i]][1] - epslnY * exp(-kappa * dl * (dist1[i] / dl - 1)) / pow(dist1[i] / dl, 2) / dl * (1 + kappa * dist1[i]) * (-dy1[i]) / dist1[i];
	//		r[jtemp[i]][0] = r[jtemp[i]][0] - epslnY * exp(-kappa * dl * (dist1[i] / dl - 1)) / pow(dist1[i] / dl, 2) / dl * (1 + kappa * dist1[i]) * (dx1[i]) / dist1[i];
	//		r[jtemp[i]][1] = r[jtemp[i]][1] - epslnY * exp(-kappa * dl * (dist1[i] / dl - 1)) / pow(dist1[i] / dl, 2) / dl * (1 + kappa * dist1[i]) * (dy1[i]) / dist1[i];

	//		r[itemp[i]][2] = r[itemp[i]][2] - epslnY * exp(-kappa * sigma_ls * (dist4[i] / sigma_ls - 1)) / pow(dist4[i] / sigma_ls, 2) / sigma_ls * (1 + kappa * dist4[i]) * (-dx4[i]) / dist4[i];
	//		r[itemp[i]][3] = r[itemp[i]][3] - epslnY * exp(-kappa * sigma_ls * (dist4[i] / sigma_ls - 1)) / pow(dist4[i] / sigma_ls, 2) / sigma_ls * (1 + kappa * dist4[i]) * (-dy4[i]) / dist4[i];
	//		r[jtemp[i]][0] = r[jtemp[i]][0] - epslnY * exp(-kappa * sigma_ls * (dist4[i] / sigma_ls - 1)) / pow(dist4[i] / sigma_ls, 2) / sigma_ls * (1 + kappa * dist4[i]) * (dx4[i]) / dist4[i];
	//		r[jtemp[i]][1] = r[jtemp[i]][1] - epslnY * exp(-kappa * sigma_ls * (dist4[i] / sigma_ls - 1)) / pow(dist4[i] / sigma_ls, 2) / sigma_ls * (1 + kappa * dist4[i]) * (dy4[i]) / dist4[i];

	//		r[itemp[i]][2] = r[itemp[i]][2] - epslnY * exp(-kappa * ds * (dist3[i] / ds - 1)) / pow(dist3[i] / ds, 2) / ds * (1 + kappa * dist3[i]) * (-dx3[i]) / dist3[i];
	//		r[itemp[i]][3] = r[itemp[i]][3] - epslnY * exp(-kappa * ds * (dist3[i] / ds - 1)) / pow(dist3[i] / ds, 2) / ds * (1 + kappa * dist3[i]) * (-dy3[i]) / dist3[i];
	//		r[jtemp[i]][2] = r[jtemp[i]][2] - epslnY * exp(-kappa * ds * (dist3[i] / ds - 1)) / pow(dist3[i] / ds, 2) / ds * (1 + kappa * dist3[i]) * (dx3[i]) / dist3[i];
	//		r[jtemp[i]][3] = r[jtemp[i]][3] - epslnY * exp(-kappa * ds * (dist3[i] / ds - 1)) / pow(dist3[i] / ds, 2) / ds * (1 + kappa * dist3[i]) * (dy3[i]) / dist3[i];

	//}




	// Torque and rotation
	// Magnetic (this is mechanism 2 in the langmuir paper 2021, not 3:lag yet)
	for (int i = 0; i < npart; i++) {
		r[i][4] = r[i][4] + lambdaMT * pow(dl, 3) * CMMl_i/ CMMl_i;
		r[i][5] = r[i][5] + lambdaMT * pow(ds, 3) * CMMs_i/ CMMl_i;
	}
	 
	//BD random torque ,only consider theta, not phi
	for (int i = 0; i < npart; i++) {
			double c = SND(generator);
			r[i][4] = r[i][4] + 1.0 * c / pow(Dt / pow(dl, 3), 0.5);
			r[i][5] = r[i][5] + 1.0 * c / pow(Dt / pow(ds, 3), 0.5);
	}



	double mcx;
	double mcy;
	double mcx_a;
	double mcy_a;
	double armlx; 
	double armly;
	double armsx;
	double armsy;
	double angle_torque;
	for (int i = 0; i < npart; i++){
		if (p[i][2] - p[i][0] > xlength / 2)
			mcx = (pow(dl, 3) * (p[i][0] + xlength) + pow(ds, 3) * p[i][2]) / (pow(dl, 3) + pow(ds, 3));
		if (p[i][0] - p[i][2] > xlength / 2)
			mcx = (pow(dl, 3) * p[i][0] + pow(ds, 3) * (p[i][2] + xlength)) / (pow(dl, 3) + pow(ds, 3));
		if (abs(p[i][0] - p[i][2]) < xlength / 2)
			mcx = (pow(dl, 3) * p[i][0] + pow(ds, 3) * p[i][2]) / (pow(dl, 3) + pow(ds, 3));

		if (p[i][3] - p[i][1] > ylength / 2)
			mcy = (pow(dl, 3) * (p[i][1] + ylength) + pow(ds, 3) * p[i][3]) / (pow(dl, 3) + pow(ds, 3));
		if (p[i][1] - p[i][3] > ylength / 2)
			mcy = (pow(dl, 3) * p[i][1] + pow(ds, 3) * (p[i][3] + ylength)) / (pow(dl, 3) + pow(ds, 3));
		if (abs(p[i][1] - p[i][3]) < ylength / 2)
			mcy = (pow(dl, 3) * p[i][1] + pow(ds, 3) * p[i][3]) / (pow(dl, 3) + pow(ds, 3));

		armlx = mcx - p[i][0] - xlength * int((mcx - p[i][0]) / (xlength / 2));
		armly = mcy - p[i][1] - ylength * int((mcy - p[i][1]) / (ylength / 2));
		armsx = mcx - p[i][2] - xlength * int((mcx - p[i][2]) / (xlength / 2));
		armsy = mcy - p[i][3] - ylength * int((mcy - p[i][3]) / (ylength / 2));

		r[i][4] = r[i][4] - r[i][1] * armlx + r[i][0] * armly;
		r[i][5] = r[i][5] - r[i][3] * armsx + r[i][2] * armsy;


		angle_torque = p[i][4] + (r[i][4] + r[i][5]) * Dt / (pow(Ksil, 3) + pow(Ksis, 3));

		mcx_a = mcx + (r[i][0] + r[i][2]) * Dt / (abs(dy[i]) + sigma_ls);
		mcy_a = mcy + (r[i][1] + r[i][3]) * Dt / (abs(dx[i]) + sigma_ls);
		r[i][0] = mcx_a - bl * pow(ds, 3) / (pow(dl, 3) + pow(ds, 3)) * cos(angle_torque) - p[i][0];
		r[i][1] = mcy_a - bl * pow(ds, 3) / (pow(dl, 3) + pow(ds, 3)) * sin(angle_torque) - p[i][1];
		r[i][2] = mcx_a - bl * pow(ds, 3) / (pow(dl, 3) + pow(ds, 3)) * cos(angle_torque) + bl * cos(angle_torque) - p[i][2];
		r[i][3] = mcy_a - bl * pow(ds, 3) / (pow(dl, 3) + pow(ds, 3)) * sin(angle_torque) + bl * sin(angle_torque) - p[i][3];
		r[i][4] = angle_torque - p[i][4];
		r[i][5] = angle_torque - p[i][5] + PI;
	}
	for (int i = 0; i < num_cent; i++) {
		r[i][0] = 0;
		r[i][1] = 0;
		r[i][2] = 0;
		r[i][3] = 0;
		r[i][4] = 0;
		r[i][5] = 0;
	}


	return r;


}
