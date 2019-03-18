package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import main.Main;
import math.Mat;
import math.MatSolver;
import math.Vect;
import math.util;

public class RunMagGear {

	
	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){


			
			int nTsteps=model.nTsteps;
		
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

			Vect T=new Vect(nTsteps);





			if(model.loadFlux) {model.magAnalysis=false;}

			
			String dispFolder="";
			String fluxFolder="";

	

			if(model.saveForce){
			
					dispFolder = System.getProperty("user.dir")+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveFlux){
			
					fluxFolder = System.getProperty("user.dir")+"\\fluxes";

				File dfolder = new File(fluxFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			
			if(nTsteps>1)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			
			
			Vect x=new Vect();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");

			//	if(model.analysisMode>=0) magAnalysis=true;


			double omega=2*PI*model.freq;

			double elAng, mechAng, elAngStep,mechAngStep;

			int nSamples=0;

			
			mechAngStep=model.rotSpeed*model.dt;
	
	
			double mechAng0=0;

			
			double ms=0e6;
			double ms2=1e6;
			
			int ir=1;
			
			if(ms>0){
			model.region[ir].setM(new Vect(1,0));
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				model.element[i].setM(model.getElementCenter(i).normalized().times(ms));
			}
		
			 ir=2;
				model.region[ir].setM(new Vect(1,0));
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				model.element[i].setM(model.getElementCenter(i).normalized().times(-ms));
			}
			}
			
			if(ms2>0){
			 ir=9;
				model.region[ir].setM(new Vect(1,0));
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
					model.element[i].setM(model.getElementCenter(i).normalized().times(ms2));
				}
				 ir=10;
					model.region[ir].setM(new Vect(1,0));
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
						model.element[i].setM(model.getElementCenter(i).normalized().times(-ms2));
					}
			
			}	
					
					
			for(int i=model.nBegin;i<=model.nEnd;i+=model.nInc){
				
				main.gui.tfX[0].setText((i+1)+"/"+nEnd);
				
				mechAng=mechAng0+i*mechAngStep;
				
				if(abs(mechAng)>0 && model.rotateRotor && model.nRotReg>0){

					model.resetMotor();
					model.setRotorPosition(mechAng);	
				//	model.setRotorIndex(i);
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


									x=model.solveMagLin(i,x);	

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

				
						this.xp=x.deepCopy();

						if(model.saveFlux){
							if(i==model.nBegin){
								model.writeMesh(fluxFolder+"\\bun"+model.nBegin+".txt");
							}
							String fluxFile = fluxFolder+"\\flux"+i+".txt";
						
							model.writeB(fluxFile);
	
						}


						util.pr(" >>>>>>>>>> Bmax >>>>>>>"+model.Bmax);

					model.resetReluctForce();
					model.setReluctForce();
				//	model.setMSForce();
					
			}
			
				

				if(model.solver.terminate) break;

			}
		

	
	}




	}

	
