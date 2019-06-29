package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.io.File;

import java.text.DecimalFormat;

import main.Main;
import math.Vect;
import math.util;

public class RunMag {

	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){
		
		if(model.POD==1){

				POD pod=new POD();
				if(model.AC)
				pod.setMagPOD_AC(model,main);	
			
	
			return;
		}
		
			String folder=model.resultFolder;


			int nTsteps=model.nTsteps;
		
			Vect T=new Vect(nTsteps);
			int ix=0;
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

			
			if(model.loadPrevMag){
				
				model.setMagBC();
				
				String initfile =folder+"\\initxMag.txt";
				model.loader.loadPrevMag(model,initfile);

				}
			



			if(model.loadFlux) {model.magAnalysis=false;}

			
			String dispFolder="";
			String fluxFolder="";

			if(model.saveDisp){
					dispFolder = folder+"\\disps"+model.defMode;
			
			
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveForce){
			
			
					dispFolder = folder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveFlux){
					fluxFolder = folder+"\\fluxes";
				
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

			
			double tav=0;
			
			int nSamples=0;


				for(int i=nBegin;i<=nEnd;i+=inc){
					
				double t0=0;
		
				
				main.gui.tfX[2].setText((i)+"/"+nEnd);
					model.setJ(t0+i*model.dt+i);	
					
					if(model.motor && model.rotateRotor && model.nRotReg>0){
						double mechAng=i*model.dt*model.rotSpeed;

						model.resetMotor();
						model.setRotorPosition(mechAng);	
				
					//	model.setRotorIndex(0);
					}

					if(!model.loadFlux){
						
						
							model.setMagBC();
							
							if(i==nBegin){
							//if(main.console)
							//	io.Console.redirectOutput(main.gui.paramArea);
								model.reportData();
						//	if(main.console)
							//	io.Console.redirectOutput(main.gui.progressArea);
							}
							
	
							if(!model.nonLin || i==nBegin/* ||model.Q!=null*/){
								
							
									if(!model.loadPrevMag){
										model.saveAp();		
										x=model.solveMagLin(i-nBegin,x);	
									

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
							else
								if(model.nonLin ){

									
									if(model.analysisMode>0 && i!=nBegin )	{	
								
										model.saveAp();				
									}
						
									x=model.solveNonLinear(x,true,i-nBegin);
								}

						if(model.analysisMode>0)
							model.setJe();
				
						if(x!=null)
							this.xp=x.deepCopy();

							if(model.saveFlux){
								String fluxFile = fluxFolder+"\\flux"+i+".txt";
								if(i==nBegin)
									model.writeMesh(fluxFolder+"\\bun"+i+".txt");
								if(model.motor && model.rotateRotor && model.nRotReg>0){
									String meshFile = fluxFolder+"\\bun"+i+".txt";	
									
									model.writeMesh(meshFile);
			
								}
							
								model.writeB(fluxFile);
								
		
							}


						model.resetReluctForce();
						model.setReluctForce();
						
						//model.setMSForce();
						
						
				}

					else if(model.loadFlux)

					{	
						
					
						String fluxFile = model.fluxFolderIn+"\\flux"+i+".txt";

					

							model.loadFlux(fluxFile,0);

							model.setReluctForce();

						//	model.setMSForce();
	
					}
					else if(model.loadPotentioal) {

						String vPot = folder+"\\flux\\vPot"+i+".txt";
						model.loadPotential(vPot);
						model.setB();

					}

					
			
					if(model.dim==2)
						model.setTorque(0,model.rm,1);
						else
							model.setTorque(model.r1,1,1);

					
					
					util.pr("torque >>>>>>>"+model.TrqZ);
					
					T.el[ix++]=model.TrqZ;
					
	
					
					
								//************************************************
				//	writeFiles=true;


					if( i==nEnd){

						model.writeMesh( folder+"\\bun"+i+".txt");
						String fluxFile =  folder+"\\flux"+i+".txt";
						model.writeB(fluxFile);
						model.writeJ0(folder+"\\J"+i+".txt");
					//	model.writeJe(folder+"\\Je"+i+".txt");

					}

					if(model.saveForce){
					//	String dispFile =dispFolder+"\\force"+mangs+".txt";
						String forceFile =dispFolder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);
					}


					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));

					if(model.solver.terminate){

						break;
					}

			if(model.solver.terminate) break;


			}
			
				util.plot(T);
				T.show();
				
				Vect errs=new Vect(model.solver.totalIter);
				
			//	util.pr(model.solver.totalIter);
				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);
				
				util.plot(errs);
			//	time.show();

				

				boolean writeInit=false;
				
				if(writeInit){
			
			
				String initfile = folder+"\\initxMag.txt";
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
				arr[arr.length-2][1]=0;
				
				arr[arr.length-1][0]=0;
				arr[arr.length-1][1]=0;



				model.writer.writeArray(arr,initfile);
			
				}




	}

}
