package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.SpMat;
import math.Vect;
import math.util;


public class StaticNonlinearMagSolver{
	int stepNumb;
	boolean usePrev=false;
	int totalNonlinIter;
	boolean relaxedNR=true;

	public StaticNonlinearMagSolver(){	}


	
	
	public Vect solve(Model model,Vect x,boolean echo,int step){

		boolean old=false;
		if(old) return solveOld(model, x, echo, step);

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		

		Vect Ci;
		SpMat L;
		

		int kx=0;
		if(model.Q!=null){
		 kx=model.Q.nCol;
		
		}
		
		
		Vect dA=new Vect(kx+model.numberOfUnknowns);
		
		model.setSolution(x);

		 x.zero();
			
		Vect[] B1,B2;

		model.solver.terminate(false);

		model.magMat.setRHS(model);
		
		double errNR=1;
		double fluxErr=1;
		
		double residual;
		boolean useFluxErr=false;

		int nonLinIter=0;
		model.errNRmax=1e-3; // not working reliably
		double convRatio=1e-3;
		if(model.Q!=null) convRatio=1e-5;
		double updatedConvCrot=1;
		
		double residual0=model.RHS.norm();
			model.solver.resRef=residual0;
		
		
	
		for( nonLinIter=0; errNR>model.errNRmax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+
						"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
				System.out.println(" Flux    error: "+dfe.format(fluxErr));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();
			

			model.setMagMat();
		

		
			if(!model.rotateConnect && model.motor&& model.hasTwoNodeNumb){
					
				if(nonLinIter==0)
					model.magMat.coupleFSMat(model);
				else
					model.magMat.reuseFSMat(model);
			}
		
				Ks=model.Hs.deepCopy();
				//Ks.shownz();

				b=model.RHS.sub(model.HkAk);
				
				residual=b.norm();



			
		if(useFluxErr)
			errNR=fluxErr;
		else
			errNR=residual/residual0;
			
		

			updatedConvCrot=errNR*convRatio;
			
			Ci=Ks.scale(b);
			
			L=Ks.ichol();
			
			
			if(b.abs().sum()>1e-11){
			//	dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);

				dA=model.solver.ICCG(Ks,L, b,updatedConvCrot,iterMax);	
			//dA=model.solver.ICCG(Ks,L, b,model.errCGmax*1e-3,iterMax);	
				//if(b.abs().max()>1e-6)
			//dA=model.solver.err0ICCG(Ks,L,b,errNR*iccgConvRatio,iterMax);	

			}

			if(model.solver.terminate) break;
			

			dA.timesVoid(Ci);
		
			x=x.add(dA);	
				
			B1=model.getAllB();

			model.setSolution(x);	
			
		
			B2=model.getAllB();
			
			fluxErr=model.getFluxErrSquared(B1,B2);

			
	

		}
		

	
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(fluxErr));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+residual0);
		
		System.out.println("=======================================================");
		System.out.println();
		
		return x;

	}

	
	public Vect solveOld(Model model,Vect x,boolean echo,int step){
		

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		double errMax=1e-2;
		boolean fluxErr=true;
		Vect Ci;
		SpMat L;


		int kx=0;
		if(model.Q!=null){
		 kx=model.Q.nCol;
		 x.zero();
		}
		
		Vect dA=new Vect(kx+model.numberOfUnknowns);
			
	
		model.magMat.setRHS(model);


		Vect[] B1,B2;

		model.solver.terminate(false);

		double err=1,err0=1;


		int nonLinIter=0;


		for( nonLinIter=0; err>errMax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(err));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();


			model.setMagMat();

			
			if(!model.rotateConnect && model.motor&& model.hasTwoNodeNumb){
					
				if(nonLinIter==0)
					model.magMat.coupleFSMat(model);
				else
					model.magMat.reuseFSMat(model);
			}
	
				Ks=model.Hs.deepCopy();

				b=model.RHS.sub(model.HkAk);
			




			Ci=Ks.scale(b);
			L=Ks.ichol();

			if(b.abs().max()>1e-11)
				//dA=model.solver.err0ICCG(Ks,L,b,1e-3*model.errCGmax,iterMax);	
					dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);

			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			


			x=x.add(dA);
			

			model.up=x.deepCopy();


			if(fluxErr){
				B1=model.getAllB();

				model.setSolution(x);	
				
				if(model.Q!=null){

					int kp=model.Q.nCol;
					Vect vp=new Vect(kp);
					for(int k=0;k<vp.length;k++)
						vp.el[k]=x.el[x.length-vp.length+k];

				Vect interfaceA=model.Q.mul(vp);
				
				for(int i=1;i<=model.numberOfEdges;i++){
					
					if(model.edgeOnFSIndices[i]>=0){	
						model.edge[i].setA(interfaceA.el[model.edgeOnFSIndices[i]]);
					}
				}

				}


				B2=model.getAllB();
			
				err=model.getDiffMax(B1,B2);//+Math.abs(dA.el[dA.length-1]);
			}
			else{

				if(nonLinIter==1){
					err0=model.RHS.norm();
					if(err0<errMax) err0=1;
				}
				err=errMax/1e-6*dA.norm()/err0;
				model.setSolution(x);	

			}


		}
		

		model.HpAp=model.HkAk.deepCopy();

		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(err));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		

		return x;
	}




}
