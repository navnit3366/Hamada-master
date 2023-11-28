/*
 * ThreeDdiffusion.hxx
 *
 *  Created on: Nov 11, 2015
 *      Author: aske
 */

#ifndef THREEDDIFFUSION_HXX_
#define THREEDDIFFUSION_HXX_

#include </home/aske/BOUT/include/bout/physicsmodel.hxx>
#include <math.h>
#include </home/aske/BOUT/include/derivs.hxx>

class ThreeDdiffusion :
public PhysicsModel {

private:
	Field3D N;

	BoutReal Dx, Dy, Dz;
	BoutReal Lx, Ly, Lz;

protected:
	int init(bool restarting);
	int rhs(BoutReal t);

};





#endif /* THREEDDIFFUSION_HXX_ */
