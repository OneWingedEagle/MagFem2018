package femSolver;

import static java.lang.Math.PI;

import fem.Model;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=false;

	public StaticLinearMagSolver(){	}

	public Vect solve(Model model, int step ){
		
		this.stepNumb=step;

	
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);

	model.magMat.setRHS(model);


if(step==0)
	model.setMagMat();


	if(!model.rotateConnect && model.motor&& model.hasTwoNodeNumb){
	model.magMat.coupleFSMat(model);

	}
			

	//=== known values go to right hand side 

	model.RHS=model.RHS.sub(model.HkAk);

	
	SpMat  Ks=model.Hs.deepCopy();


	Vect Ci=Ks.scale(model.RHS);

		L=Ks.ichol();

		if(model.RHS.abs().max()>1e-8){

			if(!usePrev || model.xp==null){
				x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax);
			}
			else{
				//x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,model.xp);
				x=model.solver.err0ICCG(Ks,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	

			}
		}

		else
			x=new Vect(x.length);

		model.xp=x.deepCopy();


		x.timesVoid(Ci);
		

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

		model.setSolution(x);	
		

		model.setB();	

		System.out.println("Bmax ( linear analysis): "+model.Bmax);

		return x;



}




}
