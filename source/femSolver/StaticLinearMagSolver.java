package femSolver;

import static java.lang.Math.PI;
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
	if(model.commonNodes!=null)	
		kk=model.commonNodes[0][0].length;
	

	int kp=kk-2;
	
	if(!model.hasTwoNodeNumb){
	 Ks=model.Hs.deepCopy();
	}
	else{
	
	

	
	Q=new Mat(2*kk,kp);
	for(int p=0;p<kp;p++){
		Q.el[p+1][p]=1;
		
		Q.el[p+kk+1][p]=1;
	}
	
/*
  	Q.el[1][0]=1;
	 Q.el[2][1]=1;
	 Q.el[5][0]=-1;
	 Q.el[6][1]=-1;*/
		 
	for(int p=0;p<kp;p++){
		//Q.el[jx][p]=Math.cos(2*PI*jx*p/kk)/kk;
			//if(jx==1 || jx==2 ){
		//Q.el[p][p]=1;
		//Q.el[p+kk][p]=1;
	}

	int jx=0;
	for(int i=1;i<=0*model.numberOfEdges;i++){
		if(model.edge[i].common){	
		//	Vect v=model.edge[i].node[0].getCoord();
//v.hshow();
	//	double tt=util.getAng(v);
			for(int p=0;p<kp;p++){
			//Q.el[jx][p]=Math.cos(2*PI*jx*p/kk)/kk;
			//	if(jx==1 || jx==2 ){
			Q.el[jx][p]=Math.cos(2*PI*jx*p/kk);
			Q.el[jx+kk][p]=Q.el[jx][p];
					
			//	Q.el[jx][jx-1]=1;
			//  Q.el[jx+kk][jx-1]=1;
			//	//}
			//}
			jx++;
		}
	}
	}
	//Q.show();
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

	Mat BtB=Q.transp().mul(model.Rs.mul(Q));
	for(int k=0;k<kp;k++)
		BtB.el[k][k]*=2;


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
	
	
	
//Vect	b1=new Vect().linspace(1,0., model.numberOfEdges);
//	b1.el[1]=2;
//	b1.el[2]=1;
//model.Fs.show();

//model.edge[7].setKnownA(.04555);
//model.edge[11].setKnownA(.04555);
//model.edge[20].setKnownA(.04555);
//model.edge[22].setKnownA(.04555);

	//Ks=model.Hs.deepCopy();
/*	Vect b2=new Vect(8);
	b2.el[1]= 0.04555;
	b2.el[2]= 0.04555;
	b2.el[5]= 0.04555;
	b2.el[6]= 0.04555;
	
	Vect b3=model.Fs.amul(b2).times(-1);*/

	//model.RHS=b2.deepCopy();
	//model.RHS.show();
	//model.RHS=model.RHS.add(b3);
	model.RHS=new Vect(b1.length+kp);
	for(int k=0;k<b1.length;k++)
		model.RHS.el[k]=b1.el[k];
	
}	


//Vect vd=Ks.diag();

//for(int i=0;i<vd.length;i++){
//	if(Math.abs(model.edge[model.unknownEdgeNumber[i+1]].node[0].getCoord().norm()-.0574)<1e-3)
//	vd.el[i]=0;	
//}
//vd.show();
	
	Ks.show();

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
//x.show();
//vp.show();
		if(Q!=null){

			Vect vp=new Vect(kp);
			for(int k=0;k<vp.length;k++)
				vp.el[k]=x.el[x.length-vp.length+k];
		Vect y1=Q.mul(vp);
	//	y1.show();
		if(kp==2){
		model.edge[7].setKnownA(y1.el[1]);
		model.edge[11].setKnownA(y1.el[2]);
		model.edge[20].setKnownA(y1.el[5]);
		model.edge[22].setKnownA(y1.el[6]);
		}
//y1.show();
		}


		model.setSolution(x);	
		

		model.setB();	

		System.out.println("Bmax ( linear analysis): "+model.Bmax);

		return x;



}




}
