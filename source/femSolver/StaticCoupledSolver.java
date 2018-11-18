package femSolver;

import fem.Model;
import math.SpMat;
import math.Vect;
import math.util;


public class StaticCoupledSolver{
	int stepNumb;
	boolean usePrev=false;

	public StaticCoupledSolver(){	}

	public void solveCoupled(Model model, Vect x){
		double t1=System.currentTimeMillis();
		
		NolninearMagSolver solver =new NolninearMagSolver();
		
		boolean twoloops=false;
		model.setStiffMat();
		model.solver.terminate(false);
		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknownEdges);

		Vect[] B1 = null,B2=null;

		double err=1,errMax=1e-3;
		int cascadeIter=0;

		for(cascadeIter=0 ;err>errMax && cascadeIter<100;cascadeIter++)

		{

			if(model.solver.terminate) break;
			System.out.println();

			System.out.println(" >>>>>>  coscade Iter: "+cascadeIter+ "  error: "+err+"  Bmax: "+model.Bmax);
			System.out.println();


			if(!twoloops){

				model.setMagMat();


				Ci=model.Hs.scale(model.RHS);

				L=model.Hs.ichol();

				double t3=System.currentTimeMillis();

				dA=model.solver.err0ICCG(model.Hs,L,model.RHS,1e-8,model.iterMax);	

				double t4=System.currentTimeMillis();

				System.out.println();
				System.out.println(" solving time: "+(t4-t3)/1000+ " seconds ");


				if(model.solver.terminate ) break;

				dA.timesVoid(Ci);

				x=x.add(dA);

			}


			if(twoloops){

				B1=model.getAllB();
				x=solver.solve(model,x,true,cascadeIter);
				B2=model.getAllB();
				err=model.getDiffMax(B1,B2);
			}
			else{
				B1=model.getAllB();
				model.setSolution(x);	
				model.setB();	
				B2=model.getAllB();
			}

			err=model.getDiffMax(B1,B2);

			model.setReluctForce();

			model.setMSForce();
			model.setDeformation();

		}


		System.out.println("    cascade Iteration: "+cascadeIter+ "  error: "+err);

		double t2=System.currentTimeMillis();
		util.pr("[][][][][][][][][][]>>>> cascade time(sec): "+(t2-t1)/1000);
	}




}
