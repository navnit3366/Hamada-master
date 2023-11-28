/*
 * ThreeDdiffusion.cxx
 *
 *  Created on: Oct 22, 2015
 *      Author: aske
 */

#include "Driftsim.hxx"

int Driftsim::init(bool restarting){

    // Get the options
    Options *options = Options::getRoot()->getSection("Driftsim");
    Options  *meshoptions = Options::getRoot()->getSection("mesh");

    phiSolver = Laplacian::create();
    phi = 0.; // Starting phi

      meshoptions->get("Lx",Lx,1.0);
      meshoptions->get("Ly",Ly,1.0);


      /*this assumes equidistant grid*/
      mesh->dx = Lx/(mesh->GlobalNx - 2*mesh->xstart);

      mesh->dy = Ly/(mesh->GlobalNy - 2*mesh->ystart);
      GRID_LOAD4(q,r,R,Pi);


      xField = FieldFactory::get()->create3D("xSec:x_", Options::getRoot(), mesh, CELL_CENTRE, 0);
      xField.setBoundary("xSec");
      xField.applyBoundary();

      yField = FieldFactory::get()->create3D("ySec:y_", Options::getRoot(), mesh, CELL_CENTRE, 0);
      yField.setBoundary("ySec");
      yField.applyBoundary();

      zField = FieldFactory::get()->create3D("zSec:z_", Options::getRoot(), mesh, CELL_CENTRE, 0);
      zField.setBoundary("zSec");
      zField.applyBoundary();

      output.write("SIZES: %d, %d, %e\n", mesh->GlobalNy, (mesh->GlobalNy - 2*mesh->ystart), mesh->dy(0,0));

      SAVE_ONCE2(Lx,Ly);


      Options *cytooptions = Options::getRoot()->getSection("cyto");
      OPTION(cytooptions, Dx, 0.0);
      OPTION(cytooptions, Dy, 0.0);
      OPTION(cytooptions, Dz, 0.0);

      SAVE_ONCE3(Dx, Dy, Dz);


      //set mesh
      output<<mesh->ngz<<std::endl;
      for(int x = 0; x < mesh->ngx-1; ++x){
    	  output<<x<<std::endl;
            for(int y = 0; y < mesh->ngy-1; ++y){
                output << "mesh->Bxy("<<x<<","<<y<<")=   " << mesh->Bxy(x,y) << std::endl;
              }
          }
      // Tell BOUT++ to solve N
      SAVE_REPEAT(phi);
      SOLVE_FOR2(N,vort);
      phi.setBoundary("phi");






    return 0;
  }

int Driftsim::rhs(BoutReal t) {
    mesh->communicate(N,vort); // Communicate guard cells
    //mesh->communicate(x_tor);
  //  mesh->communicate(y_tor);
//mesh->communicate(z_tor);
    GRID_LOAD4(q,r,R,Pi);
    BRACKET_METHOD bm;
    bm = BRACKET_ARAKAWA;


    phi = phiSolver->solve(vort, phi);
    mesh->communicate(phi);
    Dn = 1.2;
    alpha = 0.8;
    n0 = 1.1;
    Omega0 = 1.;
    DOmega = 1.3;
    Xi = mesh->Bxy/(2*Pi*q*R);
	Psibar = -(r*mesh->Bxy)/q;
	//B = (Xi*Psibar/(mesh->J))*Grad(z_tor)^Grad(x_tor);
	b = (mesh->J/mesh->g22)*(Grad(zField)^Grad(xField));
	//u = (b^(Grad(phi)))/mesh->Bxy;
	phi.applyBoundary();
    ddt(N) = Dn*Laplace(N) - N*(2/(mesh->Bxy))*bracket(mesh->Bxy, phi, bm) - bracket(phi,N,bm);
    ddt(vort) = (Omega0/n0)*(alpha*(Laplace_par(N) - Laplace_par(phi)) + DOmega*Laplace(vort)) - bracket(phi,vort,bm); //Div(b*Grad_par(N-phi))

    return 0;
  }

int main(int argc, char **argv) {                   \
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0)				      \
    return 0;                                       \
  else if (init_err > 0) 			      \
    return init_err;				      \
  try {                                             \
    Driftsim *model = new Driftsim();           \
    Solver *solver = Solver::create();              \
    solver->setModel(model);                        \
    solver->addMonitor(bout_monitor, Solver::BACK); \
    solver->outputVars(dump);                       \
    solver->solve();                                \
    delete model;                                   \
    delete solver;                                  \
  }catch (BoutException &e) {                       \
    output << "Error encountered\n";                \
    output << e.what() << endl;                     \
    MPI_Abort(BoutComm::get(), 1);                  \
  }                                                 \
  BoutFinalise();                                   \
  return 0;                                         \
}

