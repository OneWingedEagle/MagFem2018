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
			

	//=== known values go to right hand side 


	model.RHS=model.RHS.sub(model.HkAk);

	
	SpMat Ks=null;
	
/*	Mat Q=null;
	int kk=0;
	if(model.edgeOnFSIndices!=null)	
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edgeOnFSIndices[i]>=0) kk++;
		}*/
	//	kk=model.commonNodes[0][0].length;


	//int kp=kk;
	//int kph=0;
	

	//if(!model.hasTwoNodeNumb){
	 Ks=model.Hs.deepCopy();
	/*}
	else{
		
	int method=1;



	if(method==0){
	
		kp=kk/2-2;
		Q=new Mat(kk,kp);
	for(int p=0;p<kp;p++){
		Q.el[p+1][p]=1;
		
		Q.el[p+1+kk/2][p]=1;
	}
	}
 
	else if(method==1){
		kp=kk/4;
		
			
		if(kp%2==0) kp++;
	
		kp=Math.min(kp,61);
		
		Q=new Mat(kk,kp);


	double tt1=model.alpha1;
	double tt2=model.alpha2;
	double span=tt2-tt1;

	kph=(kp-1)/2;
	
	int jx=0;
	
	for(int i=1;i<=model.numberOfEdges;i++){
		if(model.edgeOnFSIndices[i]>=0){	
					boolean onBoundary=false;
			for(int k=0;k<4;k++){
			if(model.edge[i].node[0].onBound[k]){
				onBoundary=true;
				break;
			}
			}
			
			if(onBoundary) continue;
			
			jx=model.edgeOnFSIndices[i];
			
			Vect v=model.edge[i].node[0].getCoord();
		double tt=util.getAng(v)-tt1;

			for(int p=0;p<=kph;p++){
				if(p==0){
					Q.el[jx][p]=.5;
				}else{
			Q.el[jx][p]=Math.cos(2*PI*(p)*tt/span);
			
			Q.el[jx][p+kph]=Math.sin(2*PI*(p)*tt/span);
				
			}
		

			}
		
		
	}
	}
	}
	else if(method==2){
		kp=kk/4;
	
		Q=new Mat(kk,kp);
		
	double tt1=model.alpha1;
	double tt2=model.alpha2;
	double span=tt2-tt1;

	double dtt=span/(kp-1);
	
	int jx=0;
	
	for(int i=1;i<=model.numberOfEdges;i++){
		if(model.edgeOnFSIndices[i]>=0){	
					boolean onBoundary=false;
			for(int k=0;k<4;k++){
			if(model.edge[i].node[0].onBound[k]){
				onBoundary=true;
				break;
			}
			}
			
			if(onBoundary) continue;
			
			jx=model.edgeOnFSIndices[i];
			
			Vect v=model.edge[i].node[0].getCoord();
		double tt=util.getAng(v)-tt1;
		
		int it=(int)Math.floor(tt/dtt);
	
		Q.el[jx][it]=1-(tt-it*dtt)/dtt;
		Q.el[jx][it+1]=1-Q.el[jx][it];

		
		
	}
	}
	}

	 Ks=new SpMat(model.numberOfUnknowns+kp);
	for(int i=0;i<model.numberOfUnknownEdges;i++){
		Ks.row[i]=new SpVect(model.numberOfUnknowns+kp,model.Hs.row[i].nzLength);
		for(int k=0;k<model.Hs.row[i].nzLength;k++){
		Ks.row[i].el[k]=model.Hs.row[i].el[k];
		Ks.row[i].index[k]=model.Hs.row[i].index[k];

		}
	}
	
	SpMat Bs=new SpMat(kp);
	for(int i=0;i<kp;i++){
		Vect v=model.Fs.amul(Q.getColVect(i));
		SpVect vs=new SpVect(v);
		Bs.row[i]=new SpVect(model.numberOfUnknowns,vs.nzLength);
		for(int k=0;k<vs.nzLength;k++){
		Bs.row[i].el[k]=vs.el[k];
		Bs.row[i].index[k]=vs.index[k];
		}
	}

	//model.Rs.show("%8.3e");
	//model.Rs.shownz();

	SpMat BtB=new SpMat(kp);
	for(int i=0;i<kp;i++){
		Vect v1=model.Rs.amul(Q.getColVect(i));
		Vect v2=Q.transp().mul(v1);
		SpVect vs=new SpVect(v2);
		
		//vs.shownz();
		int nnz=0;
		for(int k=0;k<vs.nzLength;k++)
			if(vs.index[k]<=i){
				nnz++;
			}
		BtB.row[i]=new SpVect(kp,nnz);
		for(int k=0;k<vs.nzLength;k++){
			if(vs.index[k]<=i){
			BtB.row[i].el[k]=vs.el[k];
			BtB.row[i].index[k]=vs.index[k];

			}
		}
		//BtB.row[i].showr();
	}
	//Mat BtB=Q.transp().mul(model.Rs.mul(Q));
	//BtB.shownz();
	//BtB.show("%8.3e");
	
	for(int i=0;i<kp;i++){
		Ks.row[i+model.numberOfUnknowns]=new SpVect(model.numberOfUnknowns+kp,Bs.row[i].nzLength+BtB.row[i].nzLength);
		for(int k=0;k<Bs.row[i].nzLength;k++){
			Ks.row[i+model.numberOfUnknowns].el[k]=Bs.row[i].el[k];
			Ks.row[i+model.numberOfUnknowns].index[k]=Bs.row[i].index[k];

			}
		for(int j=0;j<BtB.row[i].nzLength;j++){
		Ks.row[i+model.numberOfUnknowns].el[Bs.row[i].nzLength+j]=BtB.row[i].el[j];

		Ks.row[i+model.numberOfUnknowns].index[Bs.row[i].nzLength+j]=BtB.row[i].index[j]+model.numberOfUnknowns;
		}
		
	}


	Vect b1=model.RHS.deepCopy();

	model.RHS=new Vect(b1.length+kp);
	for(int k=0;k<b1.length;k++)
		model.RHS.el[k]=b1.el[k];
	
}	*/
	//Ks.shownz();
	
	
//	Ks.show("%6.3e");

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
		

		if(model.Q!=null){

			int kp=model.Q.nCol;
			Vect vp=new Vect(kp);
			for(int k=0;k<vp.length;k++)
				vp.el[k]=x.el[x.length-vp.length+k];

		Vect y1=model.Q.mul(vp);

		//y1.show();


		
		for(int i=1;i<=model.numberOfEdges;i++){
			
			if(model.edgeOnFSIndices[i]>=0){	
				model.edge[i].setA(y1.el[model.edgeOnFSIndices[i]]);
			}
		}

		}


		model.setSolution(x);	
		

		model.setB();	

		System.out.println("Bmax ( linear analysis): "+model.Bmax);

		return x;



}




}
