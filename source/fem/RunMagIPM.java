package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import components.DynamicPlotter;
import components.GUI;
import main.Main;
import math.Mat;
import math.MatSolver;
import math.Vect;
import math.util;
import meshFactory.MeshManipulator;

public class RunMagIPM {

	
	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){
		
		double tStart=System.currentTimeMillis();
		
		if(model.POD>=0){
			if(model.POD==0){
				BlockMOR mor=new BlockMOR();
				mor.solve(model,main);	
			}
			else if(model.POD==1){
				POD pod=new POD();
				pod.setMagPOD(model,main);	
			}
	
			return;
		}

			String folder=model.resultFolder;


			int nTsteps=model.nTsteps;
		
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

			Vect T=new Vect(nTsteps);
			
			
			Vect time=new Vect(nTsteps);

			
			Vect Vn=new Vect(nTsteps);
			double[][] Is=new double[5][nTsteps];

			double[][] Vf=new double[4][nTsteps];
			double[][] Vs=new double[4][nTsteps];
	

	//*******************************
	 DynamicPlotter plotter=new DynamicPlotter();
			 
		boolean plot=false;
			 
			 if(plot)
			 {
				 plotter.setPlotWinodws();
				 boolean dynamicGrapph=true;
				 if(dynamicGrapph)
					 plotter.setDynamicVisibility();
				 
			 }
	  //*******************************

			

			//double kme;

		//	kme=0.25;
										
			if(model.loadPrevMag){
				
				model.setMagBC();
				
				String initfile = model.resultFolder+"\\initxMag.txt";
				model.loader.loadPrevMag(model,initfile);

				}
			

			// log mode: 0: instantaneous torque, 1: average torque vs. beta, 2: emfs

			int logMode=0;


			double dt=model.dt;

			double elAng, mechAng, elAngStep,mechAngStep;


			mechAngStep=model.rotSpeed*dt;
			
			elAngStep=model.freq*360*dt;

			
		double elAng0=0;
			
			if(!model.circuit) elAng0=280;
			//if(model.circuit) elAng0=0;

			//elAng0=0;
			
			double mechAng0=0;


			if(model.loadFlux) {model.magAnalysis=false;}

			
			String dispFolder="";
			String fluxFolder="";

		

			if(model.saveForce){
			
					dispFolder = model.resultFolder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}
			
			if(model.saveFlux){
				
				fluxFolder = model.resultFolder+"\\fluxes";
			File dfolder = new File(fluxFolder);
			if(dfolder.exists())
				util.deleteDir(dfolder);
			dfolder.mkdir();

		}
			
			if(nTsteps>1)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			boolean writeFiles=model.writeFiles;
			
			
			Vect x=new Vect();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");

			//	if(model.analysisMode>=0) magAnalysis=true;


			double omega=2*PI*model.freq;
			


			Mat W=new Mat();

			if(model.POD==-1){
				
				model.setMagBC();
				W=new Mat(model.numberOfUnknowns,model.snapShot);
				

			}

				mechAng=mechAng0;
				elAng=elAng0;

				int ix=-1;

				for(int i=nBegin;i<=nEnd;i+=inc){
					
					ix++;  
				
					double beta=0;

					elAng=elAng0+i*elAngStep;

						
						mechAng=mechAng0+i*mechAngStep/model.nInc;
						

						if(mechAng>=90.0){
								mechAng-=90;
							}


					main.gui.tfX[0].setText((i)+"/"+nEnd);


				model.setBeta(beta);
				double t0=0;
				if(omega>0)
				t0=(elAng0)/180*PI/omega;
				
					model.setJ(t0+i*model.dt/model.nInc);	
					mechAng=i*model.dt*model.rotSpeed;
					if(abs(mechAng)>0 && model.rotateRotor && model.nRotReg>0){
						model.resetMotor();
						model.setRotorPosition(mechAng);	
				
						//model.setRotorIndex(i/dd);
					}
				

					if(model.magAnalysis){
				
							model.setMagBC();
							
							if(i==nBegin){
							if(main.console)
								io.Console.redirectOutput(main.gui.paramArea);
							model.reportData();
							if(main.console)
								io.Console.redirectOutput(main.gui.progressArea);
							}
							
	
							if(!model.nonLin || i==nBegin){
								
							
									if(!model.loadPrevMag){
										model.saveAp();		

										// x=new Vect(model.numberOfUnknowns);

										x=model.solveMagLin(i-model.nBegin,x);	
									

									}
									else
									{
										Vect vk=model.getUnknownA();

								
										 x=new Vect(model.numberOfUnknowns);
										for(int j=0;j<vk.length;j++){
											x.el[j]=vk.el[j];
										}
										for(int j=0;j<model.numberOfUnknownCurrents;j++){
											int nr=model.unCurRegNumb[j];
											x.el[vk.length+j]=model.region[nr].current;

										}
										if(model.nNeutral>0)
										x.el[vk.length+model.numberOfUnknownCurrents]=model.vNeutral;
									}

								}
								else
									x=this.xp;
								
							if(model.coupled){

								model.solveCoupled(x);
							}
							
							else if(model.nonLin ){

									
									if(model.analysisMode>0 && i!=nBegin )	{	
								
										model.saveAp();				
									}
						
									
									x=model.solveNonLinear(x,true,i-nBegin);
									
								}

						
						if(model.POD==-1 && ix<model.snapShot)
						W.setCol(x, ix);
						
							this.xp=x.deepCopy();

							if(model.saveFlux){
								String fluxFile = fluxFolder+"\\flux"+i+".txt";
							
								model.writeB(fluxFile);

								
								//String vpotFile = fluxFolder+"\\A"+i+".txt";

								//model.writeA(vpotFile);
								
								String meshFile = fluxFolder+"\\bun"+i+".txt";
								
								
								model.writeMesh(meshFile);
		
							}


							util.pr(" >>>>>>>>>> Bmax >>>>>>>"+model.Bmax);

						model.resetReluctForce();
						model.setReluctForce();

						
					//	model.setMSForce();

						
				}

					else if(model.loadFlux)

					{	

						String fluxFile = model.fluxFolderIn+"\\flux"+i+".txt";
						

							model.loadFlux(fluxFile,0);

							model.setReluctForce();
							
/*							int ir=8;
							for(int ie=model.region[ir].getFirstEl();ie<=model.region[ir].getLastEl();ie++){
								model.element[ie].setDeformable(true);
								//model.element[ie].getStress().hshow();
							}
*/
							//model.setMSForce();
							
						
							//String stressFile = model.resultFolder+"\\stressMS"+i+".txt";
							
							//model.writeStress(stressFile);
							
						

					
					}
					else if(model.loadPotentioal) {

						String vPot = model.resultFolder+"\\A"+i+".txt";
						model.loadPotential(vPot);
						model.setB();

					}
					
			
					
					folder=model.resultFolder;
			
					//if(model.dim==2)
						model.setTorque(.0,model.rm,1);
						//else
							//model.setTorque(model.r1,1,1);

					
					
					util.pr("torque >>>>>>>"+model.TrqZ);
			
				T.el[ix]=model.TrqZ;

		
				if(ix>0)
					time.el[ix]=dt+time.el[ix-1];

				if(model.motor){
					calculateVoltages( model, ix, T,  Vn,Is, Vf, Vs, elAng, elAng0, mechAng, logMode);
				}
							

					if(plot){
					
						plotter.updateData(ix, T, Vn, Is, Vf, Vs);
						
					}						
					
				//	writeFiles=true;


					if( writeFiles && i==nEnd){

						model.writeMesh( folder+"\\bun"+i+".txt");
						String fluxFile =  folder+"\\flux"+i+".txt";
						model.writeB(fluxFile);
						
						String forceFile =folder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);

						//String vpotFile = folder+"\\A"+i+".txt";
						//for(int w=1;w<=model.numberOfEdges;w++)
						//	if(abs(model.edge[w].node[0].getR()-0.0547)<1e-5)
							//	util.pr(util.getAng(model.edge[w].node[0].getCoord())+"\t"+model.edge[w].A);

						//model.writeNodalA2D(vpotFile);
					
					}

					if(model.saveForce){
					//	String dispFile =dispFolder+"\\force"+mangs+".txt";
						String forceFile =dispFolder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);
					}


					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));
				//	main.gui.tfX[2].setText(this.formatter.format(model.spwmLevel));
				//	main.gui.tfX[2].setText(this.formatter.format(dang/PI*180));

					if(model.solver.terminate){

						break;
					}

			if(model.solver.terminate) break;

			}
		
			boolean writeInit=true;
			
			if(writeInit){
				writeForInitial(model);
			}
			
			if(model.POD==-1){
				String solutions=model.resultFolder +"\\solutions.txt";
				model.writer.writeArray(W.el, solutions);
			}

			util.plot(time,T);
			

	
			T.show();
			
			double tEnd=System.currentTimeMillis();
			
			util.pr("-------------------");
			util.pr("Elaspsed Time (sec):");
			util.pr((tEnd-tStart)/1000.0);
				
				Vect errs=new Vect(model.solver.totalIter);
				
				util.pr(model.solver.totalIter);
				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);
				
				//util.plot(errs);
			//	time.show();

				//util.pr("-----");
				//util.pr(T.sum()/T.length);
		
		if(logMode==2)					
		writeCurrents(model, Is, Vs);
	}

	
	private void calculateVoltages(Model model,int ix,Vect T, Vect Vn,double[][] Is,
			double[][] Vf,double[][] Vs,double elAng,double elAng0,double mechAng,int logMode){

		
		double kme;

		kme=0.25;
		
	double sumI=0,sumVs=0,sumVf=0;
		for(int j=0;j<1*model.numberOfUnknownCurrents;j++){
		
			int nr=model.unCurRegNumb[j];
			Is[j+1][ix]=1*model.region[nr].current;
			sumI+=model.region[nr].current;
			Vs[j+1][ix]=1*model.region[nr].terminalVoltage;
			sumVs+=model.region[nr].terminalVoltage;
			Vf[j+1][ix]=1*model.region[nr].inducedVoltage;
			sumVf+=model.region[nr].inducedVoltage;					
				
		}

	/*	Vs[ix][1]=1*model.region[9].terminalVoltage-model.region[11].terminalVoltage;
		Vs[ix][2]=1*model.region[11].terminalVoltage-model.region[13].terminalVoltage;
		Vs[ix][3]=1*model.region[13].terminalVoltage-model.region[9].terminalVoltage;*/


		int n1=9,n2=11,n3=13;
		
		if(model.numberOfRegions==8)
		{
			n1=1;n2=3;n3=5;
		}
		else
			if(model.numberOfRegions==4){
				n1=2;n2=2;n3=2;
			}
			
		Is[1][ix]=model.region[n1].current;
		Is[2][ix]=model.region[n2].current;
		Is[3][ix]=model.region[n3].current;

		Is[Is.length-1][ix]=sumI;

	//	Vn.el[ix]=model.vNeutral;
		
		Is[0][ix]=ix;
		Vs[0][ix]=ix;
		Vf[0][ix]=ix;
		

		double tetx=(elAng-elAng0)*kme-mechAng;
		
		
		double[] Ap=new double[1];
		double[] An=new double[1];

		if(logMode==2)
		{
			Ap=new double[1+model.numberOfElements];
			An=new double[1+model.numberOfElements];
		}

		Vn.el[ix]=Math.sin(tetx/180*PI/kme*2);

	
		if(logMode==2 && !model.circuit){
			for(int ie=1;ie<=model.numberOfElements;ie++){
				Ap[ie]=An[ie];
				An[ie]=model.getElement2DA(ie);
				
			}
			
			int[][] rg1={{9,10,11,12},{13,14,15,16},{17,18,19,20}};
			//	int[][] rg2={{10,11,12,13},{14,15,16,17},{18,19,20,21}};
			int[][] rg2={{15,16,17,18},{19,20,21,22},{23,24,25,26}};
			int[][] rg3={{9,10},{11,12},{13,14}};
			int[][] rg4={{1,2},{3,4},{5,6}};
			int[][] rg5={{2,0},{2,0},{2,0}};
			int[][] rg=new int[1][1];
			if(model.numberOfRegions==23) rg=rg1;
			else 	if(model.numberOfRegions==30) rg=rg2;
			else 	if(model.numberOfRegions==8) rg=rg4;
			else if(model.numberOfRegions==4) rg=rg5;
			else rg=rg3;
			double emfx=0;
			double[] emf=new double[3];
			double sr=0,se=0;

			for(int ip=0;ip<3;ip++)
				for(int kk=0;kk<rg[0].length;kk++){
					int ir=rg[ip][kk];
					if(ir==0) continue;
					double nLoop=model.region[ir].nloop;

						
					emfx=0;
					sr=0;

					for(int n=model.region[ir].getFirstEl();n<=model.region[ir].getLastEl();n++){						
						
						 se=model.getElementArea(n);

						 sr+=se;			
						 emfx+=(An[n]-Ap[n])*se/model.dt;

				}


					if(kk<rg[0].length/2)
						emf[ip]+=emfx*nLoop*model.height/sr;
					else
						emf[ip]-=emfx*nLoop*model.height/sr;
					
					

					
						}
					
			Vf[1][ix]=emf[0];
			Vf[2][ix]=emf[1];
			Vf[3][ix]=emf[2];


			
			
			Vs[1][ix]=Vf[1][ix]+model.region[rg[0][0]].getWireRes()*model.region[rg[0][0]].current;
			Vs[2][ix]=Vf[2][ix]+model.region[rg[1][0]].getWireRes()*model.region[rg[1][0]].current;

			Vs[3][ix]=Vf[3][ix]+model.region[rg[2][0]].getWireRes()*model.region[rg[2][0]].current;

			
			
	
		}

	
	}
		
	private void writeForInitial(Model model){
	
		
		double tet=model.tet,tetp=model.tetp,tetpp=model.tetpp;
	
		String initfile = model.resultFolder+"\\initxMag.txt";
		double[][] arr=new double[model.numberOfUnknowns+1+1][2];
		
		Vect vk=model.getUnknownA();

		Vect vp=model.getUnknownAp();

		for(int j=0;j<vp.length;j++){
			arr[j][0]=vp.el[j];
			arr[j][1]=vk.el[j];
		}
		

		for(int i=0;i<model.numberOfUnknownCurrents;i++){

			int nr=model.unCurRegNumb[i];
			
			arr[i+vp.length][0]=model.region[nr].currentp;
			arr[i+vp.length][1]=model.region[nr].current;


		}

		arr[arr.length-2][0]=model.vNeutral;
		arr[arr.length-2][1]=tet;
		
		arr[arr.length-1][0]=tetp;
		arr[arr.length-1][1]=tetpp;



		model.writer.writeArray(arr,initfile);
	
		}
	
	private void writeCurrents(Model model,double[][] Is,double[][] Vs){
		
		int nTsteps=model.nTsteps;
		
		String def = model.resultFolder + "//uu"+model.defMode+".txt";
		String Trq = model.resultFolder + "//torque.txt";
		
		String si="I";
		String sv="V";
		String s="";
		if(model.circuit) s=si;
		else
			s=sv;

		String emf1 = model.resultFolder + "//emf//"+s+"a.txt";
		String emf2 = model.resultFolder + "//emf//"+s+"b.txt";
		String emf3 = model.resultFolder + "//emf//"+s+"c.txt";
		

		try {

			PrintWriter pwu=new PrintWriter(new BufferedWriter(new FileWriter(def)));
			PrintWriter pwt=new PrintWriter(new BufferedWriter(new FileWriter(Trq)));
			PrintWriter pwem1=new PrintWriter(new BufferedWriter(new FileWriter(emf1)));
			PrintWriter pwem2=new PrintWriter(new BufferedWriter(new FileWriter(emf2)));
			PrintWriter pwem3=new PrintWriter(new BufferedWriter(new FileWriter(emf3)));

		

		

			pwem1.println("t    f(t)");
			pwem1.println("begin");
			
			pwem2.println("t    f(t)");
			pwem2.println("begin");
		
			pwem3.println("t    f(t)");
			pwem3.println("begin");
		
		if(!model.circuit){
			Vs[1][0]=Vs[1][nTsteps-1];
			Vs[2][0]=Vs[2][nTsteps-1];
			Vs[3][0]=Vs[3][nTsteps-1];
		}
			

			double ft1=0,ft2=0,ft3=0;
		for(int i=0;i<nTsteps;i++){
		double t=Vs[0][i]*model.dt;
		if(model.circuit){
			 ft1=Is[1][i];
			 ft2=Is[2][i];
			 ft3=Is[3][i];
		}
		else{
		 ft1=Vs[1][i];
		 ft2=Vs[2][i];
		 ft3=Vs[3][i];
		}
		


		pwem1.println(t+"  "+ft1);
		
		pwem2.println(t+"  "+ft2);
		
		pwem3.println(t+"  "+ft3);
		}

		pwem1.println("end");
		
		pwem2.println("end");
		
		pwem3.println("end");
		

pwu.close();
pwt.close();
pwem1.close();
pwem2.close();		
pwem3.close();


} catch (IOException exception) {
// TODO Auto-generated catch-block stub.
exception.printStackTrace();
}
	
	}
	

}
