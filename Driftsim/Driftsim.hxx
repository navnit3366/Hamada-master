/*
 * ThreeDdiffusion.hxx
 *
 *  Created on: Nov 11, 2015
 *      Author: aske
 */

#ifndef DRIFTSIM_HXX_
#define DRIFTSIM_HXX_

#include </home/aske/BOUT/include/bout/physicsmodel.hxx>
#include <math.h>
#include </home/aske/BOUT/include/derivs.hxx>
#include </home/aske/BOUT/include/invert_laplace.hxx>
#include </home/aske/BOUT/include/field_factory.hxx>


class Driftsim :
public PhysicsModel {

private:
	Laplacian* phiSolver; // Laplacian solver for vort -> phi

	Field3D N, phi, vort, x_tor, y_tor, z_tor;



	BoutReal Dx, Dy, Dz;
	BoutReal Lx, Ly, Lz;
	BoutReal r, q, R, Pi;
    BoutReal Dn, alpha, n0;
    BoutReal Omega0, DOmega;
    Field3D omega;
    Field3D xField, yField, zField;
    Field2D Xi, Psibar;
    Vector3D B, b, u;

protected:
	int init(bool restarting);
	int rhs(BoutReal t);


};





#endif /* DRIFTSIM_HXX_ */
