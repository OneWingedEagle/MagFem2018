package fem;


import static java.lang.Math.PI;

import java.io.File;
import java.text.DecimalFormat;

import Jama.Matrix;
import Jama.QRDecomposition;
import femSolver.StaticMORNonlinearMagSolver;
import femSolver.StaticNonlinearMagSolver;
import femSolver.TransientMORNonlinearMagSolver;
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
public class BlockMOR {
	private DecimalFormat formatter=new DecimalFormat("0.00");

	public BlockMOR(){}

public void solve(Model model, Main main){
		
		double tStart=System.currentTimeMillis();

		String fluxFolder="";
		
		if(model.saveFlux){
			fluxFolder = System.getProperty("user.dir")+"\\fluxes";
		
		File dfolder = new File(fluxFolder);
		if(dfolder.exists())
			util.deleteDir(dfolder);
		dfolder.mkdir();

	}
		
		model.writeMesh(fluxFolder+File.separator+"bun.txt");

		model.setMagBC();

		//String solutionfile=System.getProperty("user.dir") +"\\solutionsA180deg0.5Rough.txt";
		String solutionfile=System.getProperty("user.dir") +"\\solutionsA18LockedRotor.txt";
		//String dsolutionfile=System.getProperty("user.dir") +"\\solutionsdA180deg0.5Rough.txt";
		//String solutionfile=System.getProperty("user.dir") +"\\solutionsA90deg1Rough.txt";
		//String solutionfile=System.getProperty("user.dir") +"\\solutionsA180deg0.5RoughMath.txt";

		int L=model.numberOfUnknowns;

		int D1=model.snapShot;

		util.pr(" Using "+D1+" raw basis.");
		
		int fineStep=1;;//model.nInc;
		
		double dts=model.dt/model.nInc;

		Vect TB=new Vect(model.nTsteps*model.nInc/fineStep);
		
		StaticMORNonlinearMagSolver solver=new StaticMORNonlinearMagSolver();
		TransientMORNonlinearMagSolver transientSolver=new TransientMORNonlinearMagSolver();

		Mat	As=new Mat(model.loader.loadArrays(L,D1, solutionfile));
		//Mat	dAs=new Mat(model.loader.loadArrays(L,D1, dsolutionfile));
		
		int numBlocks=3;
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
	
		
		util.pr("NumBlocks: "+numBlocks +".  Using "+blockD+" indepenent basis per block.");
		
		Mat C=W.transp().mul(W).times(1.0/blockD);
			
		 Eigen eg2=new Eigen(C);
		 
		 Mat Q=eg2.V;
		 
		 Mat Phi=W.mul(Q);

		Vect A=new Vect(As.nRow);

			int nx0=2352;
			int cmp=0;
			nx0=1558; cmp=1; // reactor
			//nx0=15;
			nx0=24697; cmp=0;// motor half
	
		double mechAng0=0;
		double mechAng=0;
		double mechAngStep=model.rotSpeed*dts;
		double t0=280*.01/360;

	
		//for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){
		for(int i=nBlockbegin+ib*blockSteps;i<=Math.min(model.nEnd,(ib+1)*blockSteps-model.nInc);i+=model.nInc){

			int last=0;
			
			if(i+model.nInc>=model.nEnd) last=1;
			
			int ix1=(i/model.nInc)%D1;
			
			Vect A1=As.getColVect(ix1);
			
			int ix2=(ix1+1)%D1;
			
			Vect A2=As.getColVect(ix2);

			for(int j=i;j<model.nEnd+last &&j<i+model.nInc;j+=fineStep){
				
				main.gui.tfX[0].setText(Integer.toString(j)+"/"+(model.nEnd));
				main.gui.tfX[1].setText(formatter.format(model.TrqZ));

	
			double t=t0+j*dts;
			model.setJ(t);	
			
			mechAng=mechAng0+j*mechAngStep;

			if(mechAng>=90.0)
				mechAng-=90;
					
			if(model.rotateRotor){
			model.resetMotor();
			model.setRotorPosition(mechAng);
			//model.setRotorIndex(i/10);
			model.setMagBC();
			}
			//
		//	model.writeMesh(fluxFolder+"\\bun"+i+".txt");


		
			
			double alpha=(j-i)*1.0/	model.nInc;


		
			A=A1.times(1-alpha).add(A2.times(alpha));
			model.setSolution(A);
			
			model.setB();

	
			
			if(model.analysisMode==0)
			A=solver.solve(model, Phi, A, true, j);
			else{
				model.saveAp();
				A=transientSolver.solve(model, Phi, A, true, j);
			
			}
			
			
		    model.setSolution(A);
				
			model.setB();
				
		

				int nx=Math.min(nx0,model.numberOfNodes);
				
				nx=10226;

				//TB.el[jx++]=model.element[nx].getB().el[0];
				model.resetReluctForce();
				model.setReluctForce();
				model.setTorque(0,model.rm,1);
			//	util.pr("torque >>>>>>>"+model.TrqZ);
		    	TB.el[jx++]=model.TrqZ;


				if(model.saveFlux)
					if(model.saveFlux){
						String fluxFile = fluxFolder+"\\flux"+j+".txt";
					
						model.writeB(fluxFile);
						
						String bunFile = fluxFolder+"\\bun"+j+".txt";
						
						model.writeMesh(bunFile);
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
	
	
	}


















