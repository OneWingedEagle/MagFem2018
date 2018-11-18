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

		

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;

		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);

		model.magMat.setRHS(model);
		
		double errNR=1,resNR=0,resNR0=1;
		double errFlux=1;

		int nonLinIter=0;
		model.errNRmax=1e-6;

		for( nonLinIter=0; errNR>model.errNRmax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+
						"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
				System.out.println(" Flux    error: "+dfe.format(errFlux));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();

			
			model.setMagMat();
		

				Ks=model.Hs.deepCopy();

				b=model.RHS.sub(model.HkAk);
				
	

			if(nonLinIter==0){
				resNR0=b.norm();
				model.solver.resRef=resNR0;
			}
			
			resNR=b.norm();
			
		
			errNR=resNR/resNR0;


			Ci=Ks.scale(b);
			L=Ks.ichol();

			if(b.abs().max()>1e-11)
				//dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);
				
				//if(b.abs().max()>1e-6)
				dA=model.solver.err0ICCG(Ks,L,b,1e-3*model.errCGmax,iterMax);	


			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			util.pr(x.length+"  "+dA.length);

			x=x.add(dA);

			
			
			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
			errFlux=model.getDiffMax(B1,B2);

			
	

		}
		

	
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(errFlux));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+resNR);
		
		System.out.println("=======================================================");
		System.out.println();
		


		return x;

	}

	
/*	public Vect solve(Model model,Vect x,boolean echo,int step){

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		double errMax=1e-2;
		boolean fluxErr=true;
		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);

		double err=1,err0=1;


		int nonLinIter=0;


		for( nonLinIter=0; err>errMax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			

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


			Vect dv=new Vect();

			if(model.analysisMode>0)
			{


				if(model.eddyTimeIntegMode==0 || step==0)
				{


					Ks=model.Hs.addNew(model.Ss);
					Vect vp=model.getUnknownAp();

					Vect vi=model.getUnknownA();

					dv=vp.sub(vi);
					
				

					b=model.b.sub(model.HkAk).add(model.Ss.smul(dv));
				}

				else if(model.eddyTimeIntegMode==1){

//crank
					Vect vp=model.getUnknownAp();
					Vect vi=model.getUnknownA();
			

					b=model.b.add(model.bT).sub(model.HkAk.add(model.HpAp).times(0.5)).add(model.Ss.smul(vp.sub(vi)));

					Ks=model.Hs.timesNew(.5).addNew(model.Ss);

				}

				else  if(model.eddyTimeIntegMode<=-2){


					 Ks=getCircuitHs(model,step);


					int nNeut=model.nNeutral;
					int nUnCur=model.numberOfUnknownCurrents;


					Vect vk=model.getUnknownA();

					Vect vk2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vk.length;j++)
						vk2.el[j]=vk.el[j];

					
					Vect vp=model.getUnknownAp();

					Vect vp2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vp.length;j++)
						vp2.el[j]=vp.el[j];

					dv=vp2.sub(vk2);

					model.b=model.b.add(model.Ss.smul(dv));


					
					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];

						double vprev=model.region[nr].terminalVoltagep;
						if(	this.stpNumb==0)
							vprev=model.region[nr].terminalVoltage;
						
				
						double ip=model.region[nr].currentp;

		
						if(model.eddyTimeIntegMode==-3){
														
							double cf=this.theta;
							
							if(model.HpAp!=null){
								model.b=model.b.sub(model.HpAp.times(1-cf));
									}
			
					
						model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=cf*(model.lastRows[i].dot(vp2)
						+(((1-cf)*model.region[nr].getWireRes()*ip
								+(1-cf)*model.vNeutral
								-(cf*model.region[nr].terminalVoltage+(1-cf)*vprev)))*model.dt			
						-this.coilInduct*ip)/model.height;
											

						model.b=model.b.sub(model.lastRows[i].times((1-cf)*model.region[nr].currentp).vectForm());

						}
						else{
				
							model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=model.lastRows[i].dot(vp2)
							+((model.vNeutral-model.region[nr].terminalVoltage)*model.dt
							-this.coilInduct*ip)/model.height;

						}

					}

					SpMat Q=new SpMat(model.b.length,model.b.length);
					for(int i=0;i<model.b.length;i++){
						if(i<vk.length)
							Q.row[i]=null;
						else
							Q.row[i]=Ks.row[i].deepCopy();
					}

			
					Vect v3=new Vect(model.numberOfUnknowns);
					for(int j=0;j<vk.length;j++){
						v3.el[j]=vk.el[j];
					}

					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];
						v3.el[vk.length+i]=model.region[nr].current;

					}
					if(model.nNeutral>0)
					v3.el[vk.length+nUnCur]=model.vNeutral;
					
			
			
					model.b=model.b.sub(Q.smul(v3));


					if(model.eddyTimeIntegMode==-3){
				
						b=model.b.sub(model.HkAk.times((this.theta)));

					}
					else{

						b=model.b.sub(model.HkAk);
						
					}
			

				}
			}


			else
			{
				Ks=model.Hs.deepCopy();

				b=model.b.sub(model.HkAk);
			}




			Ci=Ks.scale(b);
			L=Ks.ichol();

			if(b.abs().max()>1e-6)
				dA=model.solver.err0ICCG(Ks,L,b,1e-3*model.errMax,iterMax);	

			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			
			

			int nUnCur=model.numberOfUnknownCurrents;
			int nNeut=model.nNeutral;

			// very important 
			for(int k=0;k<nUnCur+nNeut;k++){
			dA.el[dA.length-1-k]=dA.el[dA.length-1-k];
			}
			
		//	dA.show();

			x=x.add(dA);
			
			//x=x.times(.5).add(x.add(dA).times(.5));

			model.up=x.deepCopy();


			if(model.eddyTimeIntegMode<=-2){

				
				Vect vk=model.getUnknownA();

				Vect vk2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vk.length;j++)
					vk2.el[j]=vk.el[j];

				Vect vp=model.getUnknownAp();

				Vect vp2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vp.length;j++)
					vp2.el[j]=vp.el[j];

				dv=vp2.sub(vk2);

				if(nNeut>0)
					model.vNeutral=x.el[x.length-nNeut];

				for(int k=0;k<nUnCur;k++){
					int nr=model.unCurRegNumb[k];
					model.region[nr].inducedVoltage=model.lastRows[k].dot(dv)/model.dt*model.height;

					model.region[nr].current=x.el[x.length-nUnCur-nNeut+k];

				}
			}





			if(fluxErr){
				B1=model.getAllB();

				model.setSolution(x);	

				model.setB();	

				B2=model.getAllB();
			
				err=model.getDiffMax(B1,B2);//+Math.abs(dA.el[dA.length-1]);
			}
			else{

				if(nonLinIter==1){
					err0=model.b.norm();
					if(err0<errMax) err0=1;
				}
				err=errMax/1e-6*dA.norm()/err0;
				model.setSolution(x);	
				model.setB();	

			}


		}
		
		int nUnCur=model.numberOfUnknownCurrents;


		for(int k=0;k<nUnCur;k++){
			int nr=model.unCurRegNumb[k];

	 model.region[nr].currentp=model.region[nr].current;
		}


		model.HpAp=model.HkAk.deepCopy();

		System.out.println();

		System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(err));

		System.out.println();

		return x;
	}


*/

}
