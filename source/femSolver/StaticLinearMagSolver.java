package femSolver;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;

import javax.rmi.CORBA.Util;

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
			

	//=== known values go to right hand side 


	model.RHS=model.RHS.sub(model.HkAk);

	
	SpMat Ks=null;
	
	Mat Q=null;
	int kk=0;
	if(model.edgeOnFSIndices!=null)	
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edgeOnFSIndices[i]>=0) kk++;
		}
	//	kk=model.commonNodes[0][0].length;


	int kp=kk;
	

	if(!model.hasTwoNodeNumb){
	 Ks=model.Hs.deepCopy();
	}
	else{
		
	int method=1;



	if(method==0){
		Q=new Mat(kk,kp);
	for(int p=0;p<kp;p++){
		Q.el[p][p]=1;
		
		Q.el[p+kk/2][p]=1;
	}
	}
 
	else{
		kp=kk/2;
		
		if(kp%2==0) kp++;
		
		Q=new Mat(kk,kp);

	
	double tt1=model.alpha1;
	double tt2=model.alpha2;
	double span=tt2-tt1;

	int kph=(kp-1)/2;
	
	int jx=0;
	
	for(int i=1;i<=model.numberOfEdges;i++){
		if(model.edgeOnFSIndices[i]>=0){	
			jx=model.edgeOnFSIndices[i];
			Vect v=model.edge[i].node[0].getCoord();
		double tt=util.getAng(v)-tt1;

			for(int p=0;p<=kph;p++){
				if(p==0){
					Q.el[jx][p]=.5;
					//Q.el[jx+kk][p]=.5;
				}else{
			Q.el[jx][p]=Math.cos(2*PI*(p)*tt/span);
			//Q.el[jx+kk][p]=Math.cos(2*PI*(p)*tt/span);
			
			Q.el[jx][p+kph]=Math.sin(2*PI*(p)*tt/span);
			//Q.el[jx+kk][p+kph]=Math.sin(2*PI*(p)*tt/span);
				
			}
		

			}
		
	}
	}
	}
	//Q.show();
//	util.plot(Q.getColVect(2));

	 Ks=new SpMat(model.numberOfUnknowns+kp);
	 //Ks=model.Hs.deepCopy();

	for(int i=0;i<model.numberOfUnknownEdges;i++){
		Ks.row[i]=new SpVect(model.numberOfUnknowns+kp,model.Hs.row[i].nzLength);
		for(int k=0;k<model.Hs.row[i].nzLength;k++){
		Ks.row[i].el[k]=model.Hs.row[i].el[k];
		Ks.row[i].index[k]=model.Hs.row[i].index[k];

		}
	}

	Mat FQ=new Mat(kp,model.numberOfUnknownEdges);
	for(int k=0;k<kp;k++){
		Vect v=model.Fs.amul(Q.getColVect(k));
		FQ.setRow(v, k);
	}
	

	SpMat Bs=new SpMat(FQ);
	model.Rs.show("%8.3e");
	Mat BtB=Q.transp().mul(model.Rs.mul(Q));
//BtB.show("%8.3e");
	for(int i=0;i<kp;i++){
		Ks.row[i+model.numberOfUnknowns]=new SpVect(model.numberOfUnknowns+kp,Bs.row[i].nzLength+1+i);
		for(int k=0;k<Bs.row[i].nzLength;k++){
			Ks.row[i+model.numberOfUnknowns].el[k]=Bs.row[i].el[k];
			Ks.row[i+model.numberOfUnknowns].index[k]=Bs.row[i].index[k];

			}
		for(int j=0;j<=i;j++){
		Ks.row[i+model.numberOfUnknowns].el[Bs.row[i].nzLength+j]=BtB.el[i][j];

		Ks.row[i+model.numberOfUnknowns].index[Bs.row[i].nzLength+j]=j+model.numberOfUnknowns;
		}
		
	}


	Vect b1=model.RHS.deepCopy();

	model.RHS=new Vect(b1.length+kp);
	for(int k=0;k<b1.length;k++)
		model.RHS.el[k]=b1.el[k];
	
}	
	Ks.show("%6.3e");

	//Ks.shownz();

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
		

		if(Q!=null){

			Vect vp=new Vect(kp);
			for(int k=0;k<vp.length;k++)
				vp.el[k]=x.el[x.length-vp.length+k];

		Vect y1=Q.mul(vp);



		int ix=0;

		
		for(int i=1;i<=model.numberOfEdges;i++){
			
			if(model.edgeOnFSIndices[i]>=0){	
				model.edge[i].setA(y1.el[ix]);
			//	model.edge[model.edge[i].map].setA(model.edge[i].A);
				ix++;
			}
		}

		}


		model.setSolution(x);	
		

		model.setB();	

		System.out.println("Bmax ( linear analysis): "+model.Bmax);

		return x;



}




}
