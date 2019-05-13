package femSolver;

import static java.lang.Math.PI;

import fem.Model;
import io.Writer;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=false;

	public StaticLinearMagSolver(){	}

	public Vect solve(Model model,int step,Vect x_init){
		
		this.stepNumb=step;
		
		if(x_init.length==0)
			x_init=new Vect(model.numberOfUnknowns);


	
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);

	model.magMat.setRHS(model);


if(step==0){
	model.setMagMat();

}


	if(!model.rotateConnect && model.motor&& model.hasTwoNodeNumb){
	model.magMat.coupleFSMat(model);
	
	if(x_init.length<model.RHS.length){
	int kp=model.Q.nCol;
	x_init.extend(new Vect(kp));
	}
	}
	
	//=== known values go to right hand side 

	model.RHS=model.RHS.sub(model.HkAk);


	SpMat  Ks=model.Hs.deepCopy();

	if(1>10){
	//	Ks.plot();
	Mat M=Ks.matForm();
	//util.show(M.size());
	for(int i=0;i<M.nRow;i++) M.el[i][i]*=.5;
	M=M.add(M.transp());
	Writer wr=new Writer();
	wr.writeArray(M.el, "C:\\Users\\Hassan\\Desktop\\Km_hassan3.txt");
	wr.writeArray(model.RHS.el, "C:\\Users\\Hassan\\Desktop\\Fe_hassan3.txt");
	}
	


	Vect Ci=Ks.scale(model.RHS);

	x_init.timesVoid(Ci.inv());
		
	L=Ks.ichol();

		if(model.RHS.abs().max()>1e-8){

			if(!usePrev || model.xp==null){
				x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,x_init);
				//x=model.solver.CG(Ks, model.RHS,model.errCGmax,model.iterMax,x_init);
			}
			else{
				x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,model.xp);
				//x=model.solver.err0ICCG(Ks,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	

			}
		}

		else
			x=new Vect(x_init.length);

		model.xp=x.deepCopy();


		x.timesVoid(Ci);
		

		model.setSolution(x);	
		
		


		System.out.println("Bmax ( linear analysis): "+model.Bmax);
		


		return x;



}




}
