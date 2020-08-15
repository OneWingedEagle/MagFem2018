package fem;

import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

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



	 Vect nr_err;
	int[]  nr_it;

	int[][] u_index;
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
	//SpMat K_hat;
	Vect rhs_hat;
	int totalNRIter = 0;
	
	Vect [][] gp_strains;
	Vect [][] gp_stresses;
	Vect [][] gp_pl_strains;
	boolean [][] yield_states;

	boolean initialized = false;

	

	public Vect solve(Model model, SpMatSolver solver,SpMat Khat1,Vect bhat1,int step) {
		
		//this.K_hat=Khat1.deepCopy();
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
				
				u_index = new int[model.numberOfNodes + 1][model.dim];
				for (int i = 1; i <= model.numberOfNodes; i++)
					for (int k = 0; k < model.dim; k++)
						u_index[i][k] = -1;

				int ix = 0;
				for (int i = 1; i <= model.numberOfUnknownU; i++) {
					int nodeNumb = model.unknownUnumber[i];

					if (model.node[nodeNumb].isDeformable()) {
						for (int k = 0; k < model.dim; k++) {
							if (!model.node[nodeNumb].is_U_known(k)) {
								u_index[nodeNumb][k] = ix;
								ix++;
							}
						}

					}
				}

				gp_stresses=new Vect[1+model.numberOfElements][];
				gp_strains=new Vect[1+model.numberOfElements][];
				gp_pl_strains=new Vect[1+model.numberOfElements][];
				yield_states=new boolean[1+model.numberOfElements][9];
				
				 nr_err = new Vect(nLoads * nr_itmax);
				 nr_it = new int[nLoads * nr_itmax];
				 
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
										
						else
							dF=calcResidual(rhs);
					
						er = dF.norm() / rhs.norm();
						
						nr_err.el[totalNRIter]=er;

						util.pr("nr_iter: " + nr_iter + "           nr_err: " + er + "       disp_err: " + disp_err);

						if (er < nr_tol && nr_iter>0) {
							break;
						}

					

						Vect du = solveLinear(solver, Ks.deepCopy(), dF);

						disp = disp.add(du);

					//	model.setU(disp);

						if (disp.norm() > 0)
							disp_err = du.norm() / disp.norm();
						
					

			
						returnMapping(du);

						
						totalNRIter++;
					}
		

					
			}


			
			if(step==model.nTsteps-1){
			util.pr("NR error");


		
			


		}
	

			nr_err.show("%12.8e");

		model.setU(disp);

		return disp;

	}



	private void calcTangStiff() {

		
		Ks= model.mechMat.setTangStiffMat(model,false, yield_states);

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

	
	private void returnMapping(Vect du) {
		
		model.setU(du);
		
		
		Mat D =new Mat();


		Vect [] gp_delta_strain;
		Vect [] gp_delta_stress;
		
		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].isDeformable())
			model.node[i].Fms=new Vect(model.dim);
	
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			
			double yield=model.region[ir].getYield();
			
		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
			
			int[] vertNumb=model.element[i].getVertNumb();
			
			boolean deformable=false;
			for(int k=0;k<vertNumb.length;k++){
				if(model.node[vertNumb[k]].isDeformable()){
					deformable=true;
					break;
				}
			}
			if(!deformable) continue; 
			
			if(model.dim==3)
				D=model.femCalc.hook3D(model,i);
			else
				D=model.femCalc.hook(model,i);

			gp_delta_strain=model.femCalc.getGpStrainQuad(model,i);
			
			if(gp_stresses[i]==null){
			gp_stresses[i]=new Vect[gp_delta_strain.length];
			gp_strains[i]=new Vect[gp_delta_strain.length];
			gp_pl_strains[i]=new Vect[gp_delta_strain.length];
			
			for(int k=0;k<gp_delta_strain.length;k++){
				gp_stresses[i][k]=new Vect(3*(model.dim-1));
				gp_strains[i][k]=new Vect(3*(model.dim-1));
				gp_pl_strains[i][k]=new Vect(3*(model.dim-1));
			}
			}
			
			for(int k=0;k<gp_delta_strain.length;k++){
				gp_strains[i][k]=gp_strains[i][k].add(gp_delta_strain[k]);
			}
			
			gp_delta_stress=new Vect[gp_delta_strain.length];
			
			for(int k=0;k<gp_delta_stress.length;k++){
				gp_delta_stress[k]=D.mul(gp_delta_strain[k]);
				
				gp_stresses[i][k]=gp_stresses[i][k].add(gp_delta_stress[k]);
				
				double seq=calcMises(gp_stresses[i][k]);
			//	util.pr(seq+"  "+yield);
				if(seq>yield){
					yield_states[i][k]=true;
					
				///	gp_stresses[i][k]=gp_stresses[i][k].times(.999);
				}
				else yield_states[i][k]=false;
			
				
			}
		}
	}
		

	}

	
	private Vect calcResidual(Vect b) {
		
	
		Vect Fint = new Vect(b.length);


		for(int i=1;i<=model.numberOfElements;i++){
			
			int[] vertNumb=model.element[i].getVertNumb();
			
			boolean deformable=false;
			for(int k=0;k<vertNumb.length;k++){
				if(model.node[i].isDeformable()){
					deformable=true;
					if(model.node[i].Fms==null) model.node[i].Fms=new Vect(model.dim);
					break;
				}
			}
			if(!deformable) continue; 
			
	
			Vect [] gp_stress=gp_stresses[i];
	
			Vect[] nodalForce=model.femCalc.BtSigQuad(model,i,gp_stress);
			for(int j=0;j<model.nElVert;j++){
				int nn=vertNumb[j];		

				
				Vect Fms=model.node[nn].Fms;
				
				for (int k = 0; k < model.dim; k++) {

					if (!model.node[nn].is_U_known(k)) {
						int loc=u_index[nn][k] ;
						Fint.el[loc]+=nodalForce[j].el[k];
						Fms.el[k]+=nodalForce[j].el[k];
					}
					else{
						Fms.el[k]=0;//Fint.el[loc];
					}
				}
				//	model.node[nn].Fms=model.node[nn].Fms.add(nodalForce[j]);
			}
		}
		

	//	model.writeNodalField(model.resultFolder + "\\internal_force.txt", 2);


		
	//	Fint1.show();
	//	Fint.show();
	//	Vect Fint = K_hat.smul(u);
	//	Fint.show();

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
	
	


	public double calcMises(Vect sv){
		double se=0;

		if(model.dim==2){
			double s1,s2,s12,s3;
			s1=sv.el[0];
			s2=sv.el[1];
			s12=sv.el[2];
			s3=0.3*(s1+s2);

			se=pow(s1-s2,2)+pow(s2-s3,2)+pow(s1-s3,2)+6*s12;
		}
		else{
			double s1,s2,s3,s12,s23,s13;
			s1=sv.el[0];
			s2=sv.el[1];
			s3=sv.el[2];
			s12=sv.el[3];
			s23=sv.el[4];
			s13=sv.el[5];
			se=pow(s1-s2,2)+pow(s2-s3,2)+pow(s1-s3,2)+6*s12+6*s23+6*s13;
		}

		se/=2;
		se=sqrt(se);

	return se;
	}

	
	public static void main(String[] args) {

		new Main();
	}

}

