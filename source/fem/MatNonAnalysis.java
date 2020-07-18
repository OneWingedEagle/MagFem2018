package fem;

import static java.lang.Math.PI;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import io.Loader;
import main.Main;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatAsym;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;



public class MatNonAnalysis {
	

	private SpMat Ks = null;

	Mat KK = null;
	MatSolver direct_slv = null;





	Model model;

	int nr_itmax0 = 10;
	int nLoads0 = 1;
	double nr_tol0=1e-3;
	
	double nr_tol;
	
	int n_modifNR;
	int nr_itmax;
	int nLoads;
	
	Vect disp, rhs;

	Vect top;
	SpMat K_hat;
	Vect rhs_hat;
	int totalNRIter = 0;
	
	//Vect nr_err;
	//int nr_it[];

	boolean initialized = false;

	

	public Vect solve(Model model, SpMatSolver solver,SpMat Khat1,Vect bhat1,int step) {
		
		this.K_hat=Khat1.deepCopy();
		this.rhs_hat=bhat1.deepCopy();

		
		
		solver.terminate(false);
		this.model = model;

		direct_slv = new MatSolver();


		nLoads=1;
		nr_itmax=5;
		nr_tol=1e-3;
		
		KK = null;

		boolean direct = true;
	


			if(!initialized){
				
				disp = new Vect(model.Ks.nRow);
				if (disp.length > 10000)
					direct = false;
				
	
			initialized=true;
			
			}

			rhs =rhs_hat.deepCopy();

			System.out.println(" Contact analysis....");


			Vect dF = null;

			Vect residual = null;
		
			Vect load = rhs.deepCopy();
			for (int load_iter = 0; load_iter < nLoads; load_iter++) {

				util.pr("load_iter: " + load_iter);


				double factor = (load_iter + 1.) / nLoads;
				rhs = load.times(factor);


					double er = 1;
					double disp_err = 1;

	
					for (int nr_iter = 0; nr_iter < nr_itmax; nr_iter++) {

				
						calcTangStiff();
						
						if(nr_iter==0) dF=rhs.deepCopy();
						else dF=residual.deepCopy();



						Vect du = solveLinear(solver, Ks.deepCopy(), dF);

						disp = disp.add(du);

						model.setU(disp);

						if (disp.norm() > 0)
							disp_err = du.norm() / disp.norm();
						
						residual = this.calcResidual(disp, rhs);
						
						er = residual.norm() / rhs.norm();

						util.pr("nr_iter: " + nr_iter + "           nr_err: " + er + "       disp_err: " + disp_err);

			
						totalNRIter++;
						if (er < nr_tol) {
							break;
						}

					}

					
			}


			
			if(step==model.nTsteps-1){
			util.pr("NR error");


		
			


		}
	



		model.setU(disp);

		return disp;

	}



	private void calcTangStiff() {

		Ks=model.Ks;

	}

	private Vect solveLinear(SpMatSolver solver, SpMat Ks, Vect dF) {

		Vect du = new Vect(dF.length);

		model.Ci = Ks.scale(dF);

		// if(model.Ls==null)
		model.Ls = Ks.ichol();

		if (dF.abs().max() != 0) {
			 if(model.xp==null){
			du = solver.ICCG(Ks, model.Ls, dF, model.errCGmax, model.iterMax);
			 }
			 else{
			// du=solver.ICCG(Ks,model.Ls,
			// dF,model.errCGmax,model.iterMax,model.xp);
			 du=model.solver.err0ICCG(Ks,model.Ls, dF,model.errCGmax*1e-3,model.iterMax,model.xp);

			 }
		} else {
			util.pr("Solution is zero!");
			du = new Vect(Ks.nRow);
		}

		// model.xp=du.deepCopy();

		du.timesVoid(model.Ci);

		return du;
	}

	private Vect calcResidual(Vect u, Vect b) {

		Vect Fint = K_hat.smul(u);

		Vect dF = b.sub(Fint);



		return dF;

	}



	
	public Vect getDeformation(Model model, SpMatSolver solver, int mode,int step) {

		solver.terminate(false);
		System.out.println(" Static analysis....");

		double loadFactor0 = 1000;

		if (model.dim == 3)
			loadFactor0 = 10;
		
		if (model.centrigForce) {
			model.setNodalMass();
			double rpm = 7000;
			double rps = rpm / 30 * PI;
			double omeag2 = rps * rps;
			for (int i = 1; i <= model.numberOfNodes; i++) {
				Vect v = model.node[i].getCoord();
				double m = model.node[i].getNodalMass();
				Vect F = v.times(omeag2 * m/loadFactor0);
				model.node[i].setF(F);
			}

		}
		

		Vect bU1=model.bU.add(model.getbUt(mode));
		
		bU1.timesVoid((step+1.)/model.nTsteps);
		


		boolean axi = (model.struc2D == 2);

		if (axi)
			loadFactor0 *= 2 * PI;
		
		bU1.timesVoid(loadFactor0);


		Vect u=solve(model, solver,model.Ks,bU1,step);
		
		return u;
	}
	

	
	public static void main(String[] args) {

		new Main();
	}

}

