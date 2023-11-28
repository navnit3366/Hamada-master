/*
 * ThreeDdiffusion.cxx
 *
 *  Created on: Oct 22, 2015
 *      Author: aske
 */

#include "ThreeDdiffusion.hxx"

int ThreeDdiffusion::init(bool restarting){

    // Get the options
    Options *options = Options::getRoot()->getSection("ThreeDdiffusion");
    Options  *meshoptions = Options::getRoot()->getSection("mesh");

      meshoptions->get("Lx",Lx,1.0);
      meshoptions->get("Ly",Ly,1.0);


      /*this assumes equidistant grid*/
      mesh->dx = Lx/(mesh->GlobalNx - 2*mesh->xstart);

      mesh->dy = Ly/(mesh->GlobalNy - 2*mesh->ystart);


      output.write("SIZES: %d, %d, %e\n", mesh->GlobalNy, (mesh->GlobalNy - 2*mesh->ystart), mesh->dy(0,0));

      SAVE_ONCE2(Lx,Ly);


      Options *cytooptions = Options::getRoot()->getSection("cyto");
      OPTION(cytooptions, Dx, 0.0);
      OPTION(cytooptions, Dy, 20.0);
      OPTION(cytooptions, Dz, 0.0);

      SAVE_ONCE3(Dx, Dy, Dz);

      //set mesh
      output<<mesh->ngz<<std::endl;
      for(int x = 0; x < mesh->ngx-1; ++x){
    	  output<<x<<std::endl;
            for(int y = 0; y < mesh->ngy-1; ++y){
                output << "mesh->g33("<<x<<","<<y<<")=   " << mesh->g33(x,y) << std::endl;
              }
          }
      // Tell BOUT++ to solve N
      SOLVE_FOR(N);
      
      SAVE_ONCE(mesh->g_22);


    return 0;
  }

int ThreeDdiffusion::rhs(BoutReal t) {
    mesh->communicate(N); // Communicate guard cells


    ddt(N) = Dy* Laplace_par(N);//+ (Dx + Dz)* Delp2(N)

    return 0;
  }

int main(int argc, char **argv) {                   \
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0)				      \
    return 0;                                       \
  else if (init_err > 0) 			      \
    return init_err;				      \
  try {                                             \
    ThreeDdiffusion *model = new ThreeDdiffusion();           \
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

