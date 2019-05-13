package fem;


import static java.lang.Math.PI;

import java.io.File;
import java.text.DecimalFormat;

import Jama.Matrix;
import Jama.QRDecomposition;
import femSolver.ACMagSolver;
import femSolver.StaticMORNonlinearMagSolver;
import main.Main;
import math.Complex;
import math.DFT;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class POD {
	private DecimalFormat formatter=new DecimalFormat("0.00");

	public POD(){}

	
	public void setMagPOD(Model model, Main main){
		
		
		double tStart=System.currentTimeMillis();

		String fluxFolder="";
		
		if(model.saveFlux){
			fluxFolder = System.getProperty("user.dir")+"\\fluxes";
		
		File dfolder = new File(fluxFolder);
		if(dfolder.exists())
			util.deleteDir(dfolder);
		dfolder.mkdir();

	}
		
		//model.writeMesh(fluxFolder+File.separator+"bun.txt");


		model.setMagBC();


		String solutionfile=System.getProperty("user.dir") +"\\solutionsA180deg0.5Rough.txt";
		//String solutionfile=System.getProperty("user.dir") +"\\solutionsA18LockedRotor.txt";
		
		//String deltasolutionfile=System.getProperty("user.dir") +"\\solutionsdeltaA180deg0.5Rough.txt";
	//	String solutionfile=System.getProperty("user.dir") +"\\solutionsA90deg1Rough.txt";
		//solutionfile=System.getProperty("user.dir") +"\\solutions90Cogging.txt";

		int L=model.numberOfUnknowns;

		int D1=model.snapShot;


		util.pr(" Using "+D1+" raw basis.");
		
		int fineStep=1;//model.nInc;

		Vect TB=new Vect(model.nEnd-model.nBegin+1);
		

		Vect dA=new Vect(model.numberOfUnknowns);
		Vect dAp=new Vect(model.numberOfUnknowns);
		MatSolver ms=new MatSolver();
		
		Mat	As=new Mat(model.loader.loadArrays(L,D1, solutionfile));
	//	Mat	dAs=new Mat(model.loader.loadArrays(L,D1, dsolutionfile));

		
		int order=0;
		int numBlocks=1;
		int blockD=D1/numBlocks;
		int blockSteps=D1*model.nInc/numBlocks;	

		int numCompBlocks=(model.nEnd)/blockSteps+1;
	
		int jx=0;
		int ib0=model.nBegin/blockSteps;
		int nBlockbegin=model.nBegin-ib0*blockSteps;
		
		if(model.analysisMode>0){
			model.magMat.setConductMat(model);
			model.Ss.times(model.nInc);
			}

		for(int ib=ib0;ib<numCompBlocks;ib++){

		Mat W=new Mat(L,blockD);
		
		
		for(int j=0;j<blockD;j++){
			int ix=(ib%numBlocks*blockD+j);
			W.setCol(As.getColVect(ix), j);
		}

		W.normalizeColumns();
	
		
		util.pr("numBlocks: "+numBlocks +".  Using "+blockD+" indepenent basis per block.");
		
		Mat C=W.transp().mul(W).times(1.0/blockD);
			
		 Eigen eg2=new Eigen(C);
		 
		 Mat Q=eg2.V;
		 
		 Mat Phi=W.mul(Q);


		//W=indepW(W,1e-2);
		
		// W=indepWC(W,1e-6);


		 Mat PhiT=Phi.transp();

			int nx0=2352;
			int cmp=0;
			nx0=1558; cmp=1; // reactor
			//nx0=15;
			nx0=24697; cmp=0;// motor half
	
			double dts=model.dt/model.nInc;
		double mechAng0=0;
		double mechAng=0;
		double mechAngStep=model.rotSpeed*dts;

		double t0=280*.01/360;
		
		Mat Sr=null;
		
		if(model.analysisMode>0)
			Sr=PhiT.mul(model.Ss.smul(Phi));


	//	for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){
		for(int i=nBlockbegin+ib*blockSteps;i<=Math.min(model.nEnd,(ib+1)*blockSteps-model.nInc);i+=model.nInc){

		mechAng=mechAng0+i*mechAngStep;

		if(mechAng>=90.0){
				mechAng-=90;
			}

		if(model.rotateRotor){
			model.resetMotor();
			model.setRotorPosition(mechAng);
			//model.setRotorIndex(i/10);
			model.setMagBC();
			}
			//
		//model.writeMesh(model.resultFolder+File.separator+"bun.txt");
	
			int ix1=(i/model.nInc)%D1;
			Vect A1=As.getColVect(ix1);

			
			if(jx==0){

			model.setSolution(A1);
	
			model.setB();
	
			model.magMat.setPODReactMat(model,order);
		}

			Mat Mr1=PhiT.mul(model.Hs.smul(Phi));
				
			int ix2=(ix1+1)%D1;

		
			Vect A2=As.getColVect(ix2);
		
			
			mechAng=mechAng0+((i+model.nInc))*mechAngStep;
			
			if(mechAng>=90.0)
				mechAng-=90;
			
			if(model.rotateRotor){
			model.resetMotor();
			model.setRotorPosition(mechAng);
			model.setMagBC();
			}
	
			model.setSolution(A2);
	
			model.setB();

			
			model.magMat.setPODReactMat(model,order);


			Mat Mr2=PhiT.mul(model.Hs.smul(Phi));
	

		//	Vect dArp=new Vect(Mr2.nCol);
			
			double t1=t0+i*dts;
			
			double t2=t0+(i+model.nInc)*dts;
			
			int last=0;
			
			if(i+model.nInc>=model.nEnd) last=1;
			
			for(int j=i;j<model.nEnd+last &&j<i+model.nInc;j+=fineStep){
				
				main.gui.tfX[0].setText(Integer.toString(j)+"/"+(model.nEnd));
				main.gui.tfX[1].setText(formatter.format(model.TrqZ));
				
			double alpha=(j-i)*1.0/	model.nInc;

			Mat Kr=Mr1.times(1-alpha).add(Mr2.times(alpha));
			
			
			if(model.analysisMode>0){
				Kr=Kr.add(Sr);
			}
	
			double t=t0+j*dts;
		
			model.setdJfromCurrentWaves(t,t1,t2,alpha,1);
			
		

			model.magMat.setRHS(model,false);
					
			Vect br=null;
			if(model.analysisMode==0){
				br=PhiT.mul(model.RHS);
			}
			else{
				br=PhiT.mul(model.RHS.add(model.Ss.smul(dAp)));
			}

	
			Vect dAr;

			if(br.norm()>1e-6){

				dAr=ms.gaussel(Kr, br);
			}
			else
				dAr=new Vect(br.length);	

			
			dAp=dA.deepCopy();
		    dA=Phi.mul(dAr);
		     
		     Vect Ap=A1.times(1-alpha).add(A2.times(alpha));

		    Vect A=Ap.add(dA); 
		    
		//model.setPODSolution(A,alpha);
			mechAng=mechAng0+j*mechAngStep;
			if(mechAng>=90.0)
				mechAng-=90;
		    
			if(model.rotateRotor){
				model.resetMotor();
				model.setRotorPosition(mechAng);
				model.setMagBC();
				}
	

			 model.setSolution(A);
				
				model.setB();
				
		

				int nx=Math.min(nx0,model.numberOfNodes);
				
				nx=10226;

				model.resetReluctForce();
				model.setReluctForce();
				model.setTorque(0,model.rm,1);

				TB.el[jx++]=model.TrqZ;
				if(model.saveFlux)
					if(model.saveFlux){
						String fluxFile = fluxFolder+"\\flux"+j+".txt";
					
						model.writeB(fluxFile);
					}


				if(model.solver.terminate) break;

			}
			
			if(model.solver.terminate) break;

			}
		
		}

		util.plot(TB);
		
		TB.show();
		
		double tEnd=System.currentTimeMillis();
		
		util.pr("-------------------");
		util.pr("Elaspsed Time (sec):");
		util.pr((tEnd-tStart)/1000.0);
		
/*		 Complex[] Y=DFT.dft(TB.el);
		 int N=Y.length;
		  for(int k=0;k<Y.length;k++){
				 util.pr(Y[k].norm()/N);
		}*/
	}
		
	
	public void setVibPOD(Model model, Main main){

	

		util.pr("mech integ mode: "+model.timeIntegMode);
		
		//=======================================
		int[] reg8nodes=new int[1];
		model.mapnr=new int[1];
		if(model.dim==3 && model.transfer2DTo3D){
			String m2dBun= System.getProperty("user.dir") + "\\mot4th2DFine.txt";

			model.m2d=new Model(m2dBun);

			model.mapnr=new int[model.m2d.numberOfNodes+1];
			reg8nodes=model.m2d.getRegNodes(8);
			for(int u=0;u<reg8nodes.length;u++)
				model.mapnr[reg8nodes[u]]=u;
			model.forceLamin=new Vect[1][model.m2d.region[8].getNumbElements()];

			model.m2d.coordCode=1;
		}
		//=======================================


		int N=model.nTsteps;

		Vect T=new Vect(N);
		//	model.writeMesh( model.resultFolder+"\\bun.txt");

		model.setMechBC();

	//	String solutionfile=System.getProperty("user.dir") +"\\eigsReact10\\eigVects.txt";
		//String solutionfile=System.getProperty("user.dir") +"\\eigs40Half\\eigVects.txt";
		String solutionfile=System.getProperty("user.dir") +"\\solutions100motj10.txt";
	//	solutionfile=System.getProperty("user.dir") +"\\solutions.txt";

		int ix=0;

		int mode=1;
		int L=model.numberOfUnknownUcomp;

		int D=100;

		util.pr(" Using "+D+" basis.");



		Vect UU=new Vect(model.nTsteps);



		Mat	W=new Mat(model.loader.loadArrays(L, D, solutionfile));

		W.normalizeColumns();

		model.setStiffMat(true);
		
		
		Mat C=W.transp().mul(W).times(1.0/D);
	
		
		 Eigen eg=new Eigen(C);
		
		 
		 Mat Q=eg.V;
	
		 Mat Phi=W.mul(Q);

		 Mat PhiT=Phi.transp();
		 

		 
		Mat Kr=PhiT.mul(model.Ks.smul(Phi));
		Mat Mr=PhiT.mul(model.Ms.smul(Phi));

		double a1=model.rayAlpha;
		double a2=model.rayBeta;


		Mat Cr=Mr.times(a1).add(Kr.times(a2));


	
			double dt=model.dt;
			double beta=.25;
			double gama=.5;
			double b1=1./beta/Math.pow(dt,2);

			double b2=-1./beta/dt;

			double b3=1-.5/beta;
			double b4=gama*dt*b1;
			double b5=1+gama*dt*b2;
			double b6=dt*(1+gama*b3-gama);

			 Mat Kr2=Kr.add(Mr.times(b1).add(Cr.times(b4)));

			Mat Kr1=Kr.deepCopy();
			
			Kr1.lu();
			Kr2.lu();


		Vect u;


		MatSolver ms=new MatSolver();
		
	

		for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){

						
			String file=model.forceFile[ix];

			main.gui.tfX[1].setText(Integer.toString(i)+"/"+(model.nEnd));
			
			if(model.dim==3 && model.transfer2DTo3D){


						String[] files=new String[model.forceLamin.length];
						for(int k=0;k<files.length;k++){

							int m=i%1800;
							//files[k]=System.getProperty("user.dir") + "\\forcesMotMSz"+k+"\\force"+m+".txt";
							files[k]=System.getProperty("user.dir") + "\\forcesMotMSzerostress\\force"+m+".txt";

							model.m2d.loadNodalField(files[k], 1);
							for(int j=0;j<reg8nodes.length;j++){
								model.forceLamin[k][j]=model.m2d.node[reg8nodes[j]].F;
							}
						}

							model.transfer2DTo3D( file,false,true);
							
							


						}
			else{
			model.loadNodalField(file,1);
			}
		
			
			Vect bU1=PhiT.mul(model.bU.add(model.getbUt(mode)));

			
			Vect bp=new Vect();
			Vect ur;

			if(ix<2 ){


				Vect br=bU1.deepCopy();
				
				ur=ms.solvelu(Kr1, br);
	
				
				if(ix==1)
					model.ud=ur.sub(model.up).times(1.0/dt);	
				
				model.up=ur.deepCopy();
				
				model.ud=new Vect(bU1.length);
				model.udd=new Vect(bU1.length);
			}

			else{
			bp=Mr.mul(model.up.times(b1).add(model.ud.times(-b2)).add(model.udd.times(-b3)))
			.add(Cr.mul(model.up.times(b4).add(model.ud.times(-b5)).add(model.udd.times(-b6))));

		
			
			

			Vect br=bU1.add(bp);
			
			ur=ms.solvelu(Kr2, br);

	

				Vect ud1=model.ud.deepCopy();
				Vect udd1=model.udd.deepCopy();
				Vect up1;
				

				if(model.up==null) up1=new Vect(bU1.length);
				else up1=model.up.deepCopy();
			
					model.ud=ur.sub(up1).times(b4).add(ud1.times(b5)).add(udd1.times(b6));	
					model.udd=ur.sub(up1).times(b1).add(ud1.times(b2)).add(udd1.times(b3));
	

				model.up=ur.deepCopy();
			}

				//u=ur.deepCopy();
					
		u=Phi.mul(ur);


				model.setU(u);
				
	

				int nx0=2352;
				int cmp=0;
				nx0=1558; cmp=1; // reactor
				//nx0=15;
				nx0=24697; cmp=0;// motor half
				

				int nx=Math.min(nx0,model.numberOfNodes);

				UU.el[ix]=model.node[nx].getU(cmp);


				ix++;


				if(model.saveForce)
					model.writeNodalField( model.resultFolder+"\\force"+i+".txt",1);

				if(model.saveDisp)
					model.writeNodalField( model.resultFolder+"\\disp"+i+".txt",-1);

				if(model.solver.terminate) break;



			}
		
		
		UU.timesVoid(1e9);
		util.plot(UU);
		
		UU.show();
		
		}
	
	private Mat indepW(Mat W1,double eps1){

		QRDecomposition qr=new QRDecomposition(new Matrix(W1.el));
		Mat Q=new Mat(qr.getQ().getArray());
		Mat R=new Mat(qr.getR().getArray());
		R.normalizeColumns();

		int kx=0;
		for(int i=0;i<R.nRow;i++){
			if(R.rowVect(i).norm()>eps1) kx++;
		}
		Mat W=new Mat(W1.nRow,kx);
		
		 kx=0;
		for(int i=0;i<R.nCol;i++){
			if(R.rowVect(i).norm()>eps1){
				W.setCol(W1.getColVect(i),kx);
				kx++;
			}
		}
		
		return W;

	}
	
	private Mat indepWC(Mat W,double eps2){
		
		int D=W.nCol;
		
	
		Mat C=W.transp().mul(W).times(1.0/D);
		
		 Eigen eg2=new Eigen(C);
		 
			int kx=0;
			for(int i=0;i<eg2.lam.length;i++){
				if(eg2.lam.el[i]>eps2) kx++;
			}
			
			 
			 Mat Q=eg2.V;
			 
			Mat W2=new Mat(W.nRow,kx);
			kx=0;
			for(int i=0;i<eg2.lam.length;i++){
				if(eg2.lam.el[i]>eps2) {
					W2.setCol(W.getColVect(i),kx);
					kx++;
				}
			}
		 
		 return W2;
	}
	
	public void setMagPOD_AC(Model model, Main main){
		
		
		double tStart=System.currentTimeMillis();

		String fluxFolder="";
		
		if(model.saveFlux){
			fluxFolder = System.getProperty("user.dir")+"\\fluxes";
		
		File dfolder = new File(fluxFolder);
		if(dfolder.exists())
			util.deleteDir(dfolder);
		dfolder.mkdir();

	}
		
		//model.writeMesh(fluxFolder+File.separator+"bun.txt");
		

		String solutionfile=fluxFolder +"\\snaps.txt";

		model.setMagBC();
		
		int L=model.numberOfUnknowns;

		int D1=model.snapShot;
		
		//Mat	As=new Mat(model.loader.loadArrays(L,D1, solutionfile));
		Mat	As=new Mat(L,D1);


		ACMagSolver solver=new ACMagSolver();
		model.dt=1e-2;
		Vect x;
		for(int i=0;i<D1;i++){
			model.setJ(i*model.dt);
			x=solver.solve(model, i);
			util.pr(x.length+"  "+L);
			As.setCol(x, i);
			model.dt/=10;
			
		}





		util.pr(" Using "+D1+" raw basis.");
		


		

		Vect A=new Vect(model.numberOfUnknowns);
		Vect Ap=new Vect(model.numberOfUnknowns);
		MatSolver ms=new MatSolver();
		


	
		int jx=0;

		
		if(model.analysisMode>0){
			model.magMat.setConductMat(model);
			model.Ss.times(model.nInc);
			}


		Mat W=new Mat(L,D1);
		
		
		for(int j=0;j<D1;j++){
			W.setCol(As.getColVect(j), j);
		}

		W.normalizeColumns();
	
		
		
		Mat C=W.transp().mul(W).times(1.0/D1);
			
		 Eigen eg2=new Eigen(C);
		 
		 Mat Q=eg2.V;
		 
		 Mat Phi=W.mul(Q);


		//W=indepW(W,1e-2);
		
		// W=indepWC(W,1e-6);


		 Mat PhiT=Phi.transp();

	
		double dts=model.dt/model.nInc;


		Mat Mr=PhiT.mul(model.Hs.smul(Phi));


		
		Mat Sr=null;
		
		if(model.analysisMode>0)
			Sr=PhiT.mul(model.Ss.smul(Phi));
		
		if(model.analysisMode>0){
			Mr=Mr.add(Sr);
		}


		for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){


			
			double t1=i*dts;
			
			
			int last=0;
			
			if(i+model.nInc>=model.nEnd) last=1;
			
				
			main.gui.tfX[0].setText(Integer.toString(i)+"/"+(model.nEnd));
			main.gui.tfX[1].setText(formatter.format(model.TrqZ));
				
		
	
	
		//	model.setdJfromCurrentWaves(t1);
			model.setJ(t1);
		

			model.magMat.setRHS(model,false);
					
			Vect br=null;
			if(model.analysisMode==0){
				br=PhiT.mul(model.RHS);
			}
			else{
				br=PhiT.mul(model.RHS.add(model.Ss.smul(Ap)));
			}

	
			Vect Ar;

			if(br.norm()>1e-6){

				Ar=ms.gaussel(Mr, br);
			}
			else
				Ar=new Vect(br.length);	

			
			Ap=A.deepCopy();
		    A=Phi.mul(Ar);
		 

		
		    model.setSolution(A);
				
			model.setB();
				
		


				if(model.saveFlux)
					if(model.saveFlux){
						String fluxFile = fluxFolder+"\\flux"+i+".txt";
					
						model.writeB(fluxFile);
					}


				if(model.solver.terminate) break;

		
			
			if(model.solver.terminate) break;

			}
		
		

		
		double tEnd=System.currentTimeMillis();
		
		util.pr("-------------------");
		util.pr("Elaspsed Time (sec):");
		util.pr((tEnd-tStart)/1000.0);
		
/*		 Complex[] Y=DFT.dft(TB.el);
		 int N=Y.length;
		  for(int k=0;k<Y.length;k++){
				 util.pr(Y[k].norm()/N);
		}*/
	}
	
	
	}


















