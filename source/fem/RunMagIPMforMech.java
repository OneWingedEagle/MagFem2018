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
import meshFactory.MeshManipulator;

public class RunMagIPMforMech {


	private DecimalFormat formatter=new DecimalFormat("0.00");

	boolean extra=false;
	
	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){


	//	Model stf=new Model(System.getProperty("user.dir") + "\\statFrame2D.txt");

		model.m2d=new Model(System.getProperty("user.dir") + "\\statFrame2D.txt");


		model.m2d.loadStress("C:\\Works\\proj8\\resultsShrink3DLamin\\stressLamins\\stress2D"+24+".txt");
		
		//model.m2d.loadStress(System.getProperty("user.dir") + "\\resultsShrink3DLamin\\stressLamins\\stress2D"+24+".txt");

		String m3dBun;
		Model	m3d=new Model();
			
	
		
			String folder=model.resultFolder;

		//	model.timeIntegMode=1;
			
			if(model.transfer2DTo3D)
				{
					
					
					m3dBun= System.getProperty("user.dir") + "\\motor8th.txt";
				//	String m3dData= System.getProperty("user.dir") + "\\dataMot4th2DfineJ.txt";

					m3d=new Model(m3dBun);
				//	m3d.loadData(m3dData);
					m3d.coordCode=model.coordCode;
					
					String folder1 = System.getProperty("user.dir")+"\\forces3DMag";

					File dfolder = new File(folder1);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			
			
					
				}

			int nTsteps=model.nTsteps;
		
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

			Vect T=new Vect(nTsteps);
			Vect T2=new Vect(nTsteps);
			Vect time=new Vect(nTsteps);
			Vect Tdg=new Vect(2*nTsteps);
			Vect hdg=new Vect(2*nTsteps);
			
			Vect Vn=new Vect(nTsteps);
			double[][] Is=new double[5][nTsteps];

			double[][] Vf=new double[4][nTsteps];
			double[][] Vs=new double[4][nTsteps];
	
			//******************************
			Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
			int width = (int)screenSize.getWidth();
			int height = (int)screenSize.getHeight();
		
			 Plot2DPanel[] plot = new Plot2DPanel[4];
			 JFrame[] frame = new JFrame[4];
			 
		boolean plotter=false;
			 
			 if(plotter){
				 
			 double[][] XY=new double[1][2];
			
			 int hh=40;
			 
			 for(int i=0;i<plot.length;i++){
				plot[i]=new Plot2DPanel();
			 
			 }

			 plot[0].addLinePlot("Ia", Color.red, XY);
			 plot[0].addLinePlot("Ib", Color.green, XY);
			 plot[0].addLinePlot("Ic", Color.blue, XY);
			 plot[0].addLinePlot("sumI", Color.black, XY);
				plot[0].setFixedBounds(0,0,50);
				plot[0].setFixedBounds(1,-10,10);

			   frame[0] = new JFrame("Currents");
			   frame[0].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
			  frame[0].setSize(width/3,height/3);
			  frame[0].setContentPane(plot[0]);
			  frame[0].setLocation(width/3,height/3-hh);
			
			  
				plot[1] = new Plot2DPanel();
	
				plot[1].addLinePlot("Va", Color.red, XY);
				plot[1].addLinePlot("Vb", Color.green, XY);
				plot[1].addLinePlot("Vc", Color.blue, XY);
				plot[1].addLinePlot("Vfa", Color.cyan, XY);
				plot[1].addLinePlot("Vfb", Color.pink, XY);
				plot[1].addLinePlot("Vfc", Color.gray, XY);
					plot[1].setFixedBounds(0,0,50);
					plot[1].setFixedBounds(1,-100,100);

				 frame[1] = new JFrame("Voltages");
				  frame[1].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
					 frame[1].setSize(width/3,height/3);
				 frame[1].setContentPane(plot[1]);
				 frame[1].setLocation(2*width/3,height/3-hh);

		
				  
				   plot[2] = new Plot2DPanel();
			
					 plot[2].addLinePlot("Torque", Color.red, XY);
					
						plot[2].setFixedBounds(0,0,50);
						plot[2].setFixedBounds(1,-50,50);
						

					   frame[2] = new JFrame("Torque");
					   frame[2].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
					  frame[2].setSize(width/3,height/3);
					  frame[2].setContentPane(plot[2]);
					  frame[2].setLocation(width/3,2*height/3-hh);

			  
					   plot[3] = new Plot2DPanel();
						
						 plot[3].addLinePlot("Vg", Color.red, XY);
						
							plot[3].setFixedBounds(0,0,50);
							plot[3].setFixedBounds(1,-1.2,1.2);
							
						  frame[3] = new JFrame("Vg");
						  frame[3].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
						 frame[3].setSize(width/3,height/3);
						 frame[3].setContentPane(plot[3]);
						 frame[3].setLocation(2*width/3,2*height/3-hh);

			 }
						  
						 boolean dynamicGrapph=plotter && true;
						  
						  if(dynamicGrapph){
						  
						  frame[0].setVisible(true);
						 frame[1].setVisible(true);
						  frame[2].setVisible(true);

						 frame[3].setVisible(true);
						  }

			  //*******************************
			  
			Vect de=new Vect(nTsteps);
			

			double kme;

			kme=0.25;
										
			if(model.loadPrevMag){
				
				model.setMagBC();
				
				String initfile = System.getProperty("user.dir")+"\\initxMag.txt";
				model.loader.loadPrevMag(model,initfile);

				}
			

			// log mode: 0: instantaneous torque, 1: average torque vs. beta, 2: emfs

			int logMode=0;
			
	
			double tet=model.tet,tetp=model.tetp,tetpp=model.tetpp;
			double up=0,upp=0,uu=0,vv=0,aa=0,bb=0,cc=0,ff=0,ffp=0,udgp=0,vdgp=0,dt=0,dtp;
			double mp=.075,ks=1,cs=.05;

		
			uu=.0;
			up=uu;
			mp=.1;ks=5;cs=.0;
			//double dt=model.dt;
			
			
			MatSolver msolver=new MatSolver();
			
		
			double dx=0,dxp=0,dang=0;

			double elAng, mechAng, elAngStep,mechAngStep;

			int nSamples=0;

			
			mechAngStep=model.rotSpeed*model.dt;
			
			elAngStep=mechAngStep/kme;
			
			
		double elAng0=0*260;
			
			if(model.circuit) elAng0=280;
			//if(model.circuit) elAng0=0;

			//elAng0=0;
			
			double mechAng0=0;
	
					
			double[] Ap=new double[1];
			double[] An=new double[1];

			if(logMode==2)
			{
				Ap=new double[1+model.numberOfElements];
				An=new double[1+model.numberOfElements];
			}


			if(model.loadFlux) {model.magAnalysis=false;}

			
			String dispFolder="";
			String fluxFolder="";

			if(model.saveDisp){
				if(model.fullMotor)
					dispFolder = System.getProperty("user.dir")+"\\disp6AFull-0-90mode"+model.defMode;
				else
					dispFolder = System.getProperty("user.dir")+"\\disp6A4th-0-90mode"+model.defMode;
			
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveForce){
			
					dispFolder = System.getProperty("user.dir")+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveFlux){
				if(model.fullMotor)
					fluxFolder = System.getProperty("user.dir")+"\\flux6AFull-0-90";
				else
					fluxFolder = System.getProperty("user.dir")+"\\flux6A4th-0-90";
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
			

			int nM=8;

			double tav=0;
			
			//String bun=System.getProperty("user.dir") + "\\plungers\\plunger2Dv.txt";
			String bun=System.getProperty("user.dir") + "\\plungers\\EM.txt";
			if(model.dim==3)
				bun=System.getProperty("user.dir") + "\\plungers\\EMSlice.txt";
			
			Model model0=new Model(bun);
			
		for(int m=0;m<=1*nM;m++){
			


				if(logMode!=1)
					if(m!=0*nM) continue;

				//main.gui.tfX[2].setText(Integer.toString(m)+"/"+(nM)+":"+tav);

				tav=0;
				nSamples=0;
				mechAng=mechAng0;
				elAng=elAng0;

				int iy=0;
				int ix=-1;
				
				MeshManipulator mf=new MeshManipulator();
				
				double dtfact=1;
				

				for(int i=nBegin;i<=nEnd;i+=inc){
					
		

				
					double dt0=5e-2/5;
				//	dt0=1;
				
				
					ix++;
		
			
				//	model.dt=dtfact*dt0;
									
					double umax=.0095;
					
					umax=100;
					
					if(uu>umax){
						double ut=uu;
						uu=up;
						up=ut;
						vv=-vv;
						
				
					}
					
					/*if(umax-uu<2e-3) dtfact=.2;
					else*/
						dtfact=1;
				
					//else dtfact=Math.exp(-1e4*abs(uu-up));

					dtp=dt;

					 dt=model.dt;
					double dt2=dt*dt;
					double dt3=dt2*dt;
					double dt4=dt2*dt2;
					
				 int pp=200;
						 
					 
				if(pp==1){

					 //uu=ix*.0004;
					 if(uu>.009) {
						 uu=.009; 
						// break;
					 }
					 
					 mf.meshPlunger(uu,true);

					try {
						Thread.sleep(200);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					Model modelp=model.deepCopy();
					modelp.femCalc=new Calculator(modelp);
					
					 bun=System.getProperty("user.dir") + "\\plungers\\plungerX.txt";
					model.loadMesh(bun);
					
					String 	dataFile= System.getProperty("user.dir") + "\\dataPlunger2D.txt";
					model.loadData(dataFile);
				
				
					model.dt=dtfact*dt0;
					
					 dt=model.dt;
					 dt2=dt*dt;
					 dt3=dt2*dt;
					dt4=dt2*dt2;
					
				
					int[] ncr={1};
					int nt=0;
					for(int j=0;j<ncr.length;j++)
						nt+=model.getRegNodes(ncr[j]).length;
					
					int[] ee=new int[nt];
					boolean[] ec=new boolean[model.numberOfEdges+1];
					int kx=0;
					for(int j=0;j<ncr.length;j++)
					for(int k=model.region[ncr[j]].getFirstEl();k<=model.region[ncr[j]].getLastEl();k++){
						for(int jj=0;jj<model.nElEdge;jj++){
							int en=model.element[k].getEdgeNumb(jj);
							if(!ec[en]){
								ec[en]=true;
								ee[kx++]=en;
							}
						}
						
					}
					
					for(int k:ee)
					{
						
						Vect P=model.edge[k].node[0].getCoord();

						double[] AnAp=modelp.getApAnAt(P);
						
						model.edge[k].Ap=AnAp[0];
						model.edge[k].A=AnAp[1];

						
					}
				}
				 else if(pp==2){

				
					/*Model modelp=model.deepCopy();
					modelp.femCalc=new Calculator(modelp);
					*/
					model.resetReluctForce();
					
				/*	for(int k=1;k<model.numberOfEdges;k++){
						model.edge[k].A=0;
						model.edge[k].Ap=0;
					}
	
						*/
					
				
					if(abs(uu)>1e4 ){
				
					//	util.pr(uu);
						main.gui.tfX[1].setText(this.formatter.format(uu*1e6));
								
						
						main.gui.tfX[2].setText(this.formatter.format(dtfact));

						
						int n1=1561;
						int n2=1834;
						int n3=1210;
						int n4=2146;
						
						int cmp=1;
						if(model.dim==3){
							cmp=2;
							 n1=1601;
							 n2=1881;
							 n3=1241;
							 n4=2201;
							
						}
						

					double y1=model0.node[n1].getCoord(cmp);
					double y2=model0.node[n2].getCoord(cmp);
					double y3=model0.node[n3].getCoord(cmp);
					double y4=model0.node[n4].getCoord(cmp);
					
				
					
					for(int k=1;k<model.numberOfNodes;k++){
						double y=model0.node[k].getCoord(cmp);
						if(y>=y1 && y<=y2) 
							model.node[k].setCoord(cmp,y-uu);
						else if(y<y1 && y>y3)  model.node[k].setCoord(cmp,y-(y-y3)*uu/(y1-y3));
						else if(y>y2 && y<y4)  model.node[k].setCoord(cmp,y-(y-y4)*uu/(y2-y4));
					}
						
					}


					
				
					model.dt=dtfact*dt0;
					
					 dt=model.dt;
					 dt2=dt*dt;
					 dt3=dt2*dt;
					dt4=dt2*dt2;
				
					int[] ncr={};
					int nt=0;
					
					for(int j=0;j<ncr.length;j++)
						nt+=model.getRegNodes(ncr[j]).length;
					
					int[] ee=new int[nt];

					
					boolean[] ec=new boolean[model.numberOfEdges+1];
					
					
					int kx=0;
					for(int j=0;j<ncr.length;j++)
					for(int k=model.region[ncr[j]].getFirstEl();k<=model.region[ncr[j]].getLastEl();k++){
						for(int jj=0;jj<model.nElEdge;jj++){
							int en=model.element[k].getEdgeNumb(jj);

							if(!ec[en]){
								ec[en]=true;
								ee[kx++]=en;
								
							}
						}
						
					}

					
					
			/*		for(int k:ee)
					{
						
						int nn=model.edge[k].endNodeNumber[0];
						Vect P=model.node[nn].getCoord();

						double[] AnAp=modelp.getApAnAt(P);
						
			
						model.edge[k].Ap=AnAp[0];
						model.edge[k].A=AnAp[1];

						
					}*/
				}
		
					  
				
					nSamples++;
					double beta=m*.5/nM*PI;

					elAng=elAng0+i*elAngStep;

					double tets=(mechAng0+i*mechAngStep)*PI/180;
					
						dxp=dx;

						dx=dxp+1e-3*(model.TrqZ-1*9.5*(1-Math.exp(-ix*.05)));
					
			

						tetpp=tetp;
						tetp=tet;
						double tetd=0;
						
						double dTr=(model.TrqZ-10)+20*(tets-tet);
						
						if( !model.loadPrevMag) 
							tet=tets;			
					
						
		
						
						tet=tets;
					
						tetd=180*tet/PI;
						
						mechAng=mechAng0+i*mechAngStep;
						
				model.spwmLevel=.9;
					main.gui.tfX[0].setText((i+1)+"/"+nEnd);

					//elAngRad=elAng*PI/180;	

				model.setBeta(beta);
					
				double t0=0;
				if(omega>0)
				t0=(elAng0)/180*PI/omega;
				
				
					model.setJ(t0+(nSamples-1)*model.dt);	

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

										// x=new Vect(model.numberOfUnknowns);

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

					/*	if(model.analysisMode>0)
							model.setJe();
*/
						
							this.xp=x.deepCopy();

							if(model.saveFlux){
								String fluxFile = fluxFolder+"\\flux"+i+".txt";
							
								model.writeB(fluxFile);
		
							}


							util.pr(" >>>>>>>>>> Bmax >>>>>>>"+model.Bmax);

						model.resetReluctForce();
						model.setReluctForce();
						model.setMSForce();

						
				}

					else if(model.loadFlux)

					{	
					
					//	String fluxFile = System.getProperty("user.dir")+"\\flux-0-1800CalfineOK\\flux"+i+".txt";

						String fluxFile = model.fluxFolderIn+"\\flux"+i+".txt";
						

							model.loadFlux(fluxFile,0);

							model.setReluctForce();
							
							int ir=8;
							for(int ie=model.region[ir].getFirstEl();ie<=model.region[ir].getLastEl();ie++){
								model.element[ie].setDeformable(true);
								//model.element[ie].getStress().hshow();
							}

							model.setMSForce();
							
						
							String stressFile = model.resultFolder+"\\stressMS"+i+".txt";
							
							model.writeStress(stressFile);
							
							if(model.transfer2DTo3D){

							transfer2DTo3D( model, m3d,false, true);
							

							m3d.setTorque(0,1,1);
							
							String forceFile =System.getProperty("user.dir")+"\\forces3DMag\\force"+i+".txt";
							m3d.writeNodalField(forceFile,model.forceCalcMode);
							}
						

					
					}
					else if(model.loadPotentioal) {

						String vPot = System.getProperty("user.dir")+"\\flux\\vPot"+i+".txt";
						model.loadPotential(vPot);
						model.setB();

					}
					
			
					
					folder=System.getProperty("user.dir") + "\\results";
			
					if(model.dim==2)
						model.setTorque(0,model.rm,1);
						else
							model.setTorque(model.r1,1,1);

					
					
					util.pr("torque >>>>>>>"+model.TrqZ);
			
				T.el[ix]=model.TrqZ;
				
				int cmp=model.dim-1;
				int ne=425;
				if(cmp==1) ne=412;
				
		/*		if(model.circuit)
				T.el[ix]=model.region[model.unCurRegNumb[0]].current;
				else
				T.el[ix]=model.element[ne].getB().el[cmp];*/
				
				int[] nn=model.getRegNodes(4);
				//int[] nn={1009,1010,1011,1012,1013,1014,1015};
				Vect F=new Vect(model.dim);
				
				boolean spring=true;
				ffp=ff;

				if(spring) ff=-10*mp*Math.cos(ix*dt0*10);
				else
				{
				if(model.axiSym){
				for(int u=0;u<nn.length;u++)
					if(model.node[nn[u]].F!=null)
					F=F.add(model.node[nn[u]].F.times(2*PI*abs(model.node[nn[u]].getCoord(0))));
			
				ff=-F.el[1];

				}
				else{
					for(int u=0;u<nn.length;u++)
					if(model.node[nn[u]].F!=null)
						F=F.add(model.node[nn[u]].F);
					
					ff=-F.el[2]*360;
				}
					
				
				}



				int mdf=4;

				
				T2.el[ix]=uu;
				
				if(mdf==0){

					
					double kdt=1;
					if(dtp>0)
					 kdt=dt/dtp;
					//implicit euler
					double numax=dt2*ff+mp*((1+kdt)*uu-kdt*up)+cs*dt*uu;
					double denum=mp+cs*dt+ks*dt2;
					upp=up;
					up=uu;
				//	if(ix<-2) uu=.01;
					//else
					uu=numax/denum;
			
				}
				else if(mdf==-1)
				{
					//GN1
				
					double bet1;
					
					bet1=.5;
					
					double ubar=uu+(1-bet1)*dt*vv;
					double vbar=vv+dt*aa;
					
					
					double numax=ff-(ks*ubar);
					double denum=mp+cs*bet1*dt;

					 aa=numax/denum;
			
					uu=ubar+bet1*dt*aa;
					vv=vbar+bet1*dt*aa;
						
					}
				else if(mdf==1)
				{
					//S22
				
					double thet1,thet2;
					
					thet1=.5;
					thet2=.5;
					
					double ubar=uu+thet1*dt*vv;
					
					double fbar=0;

					fbar=(ffp+ff)/2;
					
					
					double numax=fbar-(cs*vv+ks*ubar);
					double denum=mp+cs*thet1*dt+ks*thet2*dt2/2;

					 aa=numax/denum;
						if(ix<-2) {uu=.01;
						vv=0;
						
						}else{
					uu=uu+vv*dt+dt2*aa/2;
					vv=vv+dt*aa;
						}
					}
				else if(mdf==2){
					//velocity Verlet
					
					double a=(ff-cs*vv-ks*uu)/mp;
					
					uu=uu+vv*dt+a*dt2/2;
					
					double vvh=vv+a*dt/2;
					
					a=(ff-cs*vvh-ks*uu)/mp;
					
					vv=vvh+a*dt/2;
				}
				else if(mdf==3){
					//central diff
					
					double numax=dt2*ffp+mp*(2*uu-up)+0.5*cs*dt*up-ks*dt2*uu;
					double denum=mp+0.5*cs*dt;
					up=uu;
					if(ix<-2) uu=.01;
					else
					uu=numax/denum;
				
				}
				else if(mdf==4){
					//GN22
										
					double bet1,bet2;
					
					bet1=.5;
					bet2=.5;
					
					double ubar=uu+dt*vv+(1-bet2)*dt2/2*aa;
					
					double vbar=vv+(1-bet1)*dt*aa;
					
					double numax=ff-(cs*vbar+ks*ubar);
					double denum=mp+cs*bet1*dt+ks*bet2*dt2/2;

					 aa=numax/denum;
			
						up=uu;

					uu=ubar+bet2*dt2*aa/2;
					vv=vbar+bet1*dt*aa;
						
			
				}

				else if(mdf==5){
					//GN32
									
					double bet1,bet2,bet3;
					
					bet1=.5;
					bet2=.5;
					bet3=.5;
					
					double ubar=uu+dt*vv+aa*dt2/2+(1-bet3)*dt3/6*bb;
					
					double vbar=vv+dt*aa+(1-bet2)*dt2/2*bb;
					
					double abar=aa+(1-bet1)*dt*bb;
					
					double numax=ff-(mp*abar+cs*vbar+ks*ubar);
					double denum=bet1*mp*dt+cs*bet2*dt2/2+ks*bet3*dt3/6;

					 bb=numax/denum;
					
					uu=ubar+bet3*dt3*bb/6;
					vv=vbar+bet2*dt2*bb/2;
					aa=abar+bet1*dt*bb;
				
				}
				
				else if(mdf==6){
					//GN44
									
					double bet1,bet2,bet3,bet4;
					
					bet1=.5;
					bet2=.5;
					bet3=.5;
					bet4=.5;
					
					double ubar=uu+dt*vv+aa*dt2/2+dt3/6*bb+(1-bet4)*dt4/24*cc;
					
					double vbar=vv+dt*aa+dt2/2*bb+(1-bet3)*dt3/6*cc;
					
					double abar=aa+dt*bb+(1-bet2)*dt2/2*cc;
					
					double bbar=bb+(1-bet1)*dt*cc;
					
					double numax=ff-(mp*abar+cs*vbar+ks*ubar);
					double denum=bet2*mp*dt2/2+cs*bet3*dt3/6+ks*bet4*dt4/24;

					 cc=numax/denum;
					
					uu=ubar+bet4*dt4*cc/24;
					vv=vbar+bet3*dt3*cc/6;
					aa=abar+bet2*dt2*cc/2;
					bb=bbar+bet1*dt*cc;
			
				
				}
				else if (mdf==7){
					//DG
					
					if(ix==0) util.pr(uu+"  <<<<<<<<<<<<<<<<<");
					
					Mat A=new Mat(4,4);
					A.el[0][0]=mp;
					A.el[0][1]=mp;
					A.el[0][2]=cs+ks*dt;
					A.el[0][3]=cs+ks*dt/2;
					
					A.el[1][0]=dt;					
					A.el[1][1]=dt/2;
					A.el[1][2]=-1;
					A.el[1][3]=-1;
					
					A.el[2][0]=0;
					A.el[2][1]=mp/2;
					A.el[2][2]=ks*dt/2;
					A.el[2][3]=ks*dt/3+cs/2;
					
					A.el[3][0]=dt/2;
					A.el[3][1]=dt/3;
					A.el[3][2]=0;
					A.el[3][3]=-.5;
					
					
					Vect b=new Vect(4);


					b.el[0]=(ffp+ff)*dt/2+mp*vv+cs*uu;
							
					b.el[1]=-uu;
					b.el[2]=(ffp/6+ff/3)*dt;	
					
					Vect y=msolver.gaussel(A, b);
					
					T2.el[ix]=(udgp+uu)/2;

					vdgp=y.el[0];
					udgp=y.el[2];
					vv=vdgp+y.el[1];
					uu=udgp+y.el[3];
					
					hdg.el[iy]=ix*dt;
					Tdg.el[iy++]=udgp;
					hdg.el[iy]=(ix+1)*dt;

					Tdg.el[iy++]=uu;

					
				}
				else if (mdf==8){
					//DG from Li paper

					
					Mat A=new Mat(4,4);
					A.el[0][0]=ks/2;
					A.el[0][1]=ks/2;
					A.el[0][2]=-ks*dt/3;
					A.el[0][3]=-ks*dt/6;
					
					A.el[1][0]=-ks/2;					
					A.el[1][1]=ks/2;
					A.el[1][2]=-ks*dt/6;
					A.el[1][3]=-ks*dt/3;
					
					A.el[2][0]=ks*dt/3;
					A.el[2][1]=ks*dt/6;
					A.el[2][2]=mp/2+cs*dt/3;
					A.el[2][3]=mp/2+cs*dt/6;
					
					A.el[3][0]=ks*dt/6;
					A.el[3][1]=ks*dt/3;
					A.el[3][2]=-mp/2+cs*dt/6;;
					A.el[3][3]=mp/2+cs*dt/3;;
					
					
		Vect b=new Vect(4);
					
		
		b.el[0]=ks*uu;
		b.el[2]=(ffp+ff)*dt/2-(ffp/6+ff/3)*dt+mp*vv;
		b.el[3]=(ffp/6+ff/3)*dt;
				
					
					Vect y=msolver.gaussel(A, b);
					
					hdg.el[2*ix]=ix*dt;
					Tdg.el[2*ix]=uu;
					udgp=y.el[0];
					uu=y.el[1];
					vdgp=y.el[2];				
					vv=y.el[3];
				
				
					Tdg.el[2*ix+1]=udgp;
					hdg.el[2*ix+1]=ix*dt;

					T2.el[ix]=(udgp+uu)/2;
					

					
				}
	
				else if(mdf==9)
				{
					
					// dg for first order equation
					double c1=cs;
					double k1=ks;
					
					Mat B=new Mat(2,2);
					B.el[0][0]=c1+k1*dt;
					B.el[0][1]=c1+k1*dt/2;
					B.el[1][0]=k1*dt/2;
					B.el[1][1]=c1/2+k1*dt/3;
					
					Vect h=new Vect(2);
					
					h.el[0]=(ffp+ff)*dt/2+c1*uu;
					h.el[1]=(ffp/6+ff/3)*dt;			
	
					Vect w=msolver.gaussel(B, h);

					udgp=w.el[0];
					uu=udgp+w.el[1];


					Tdg.el[2*ix]=udgp;
					Tdg.el[2*ix+1]=uu;


					hdg.el[2*ix]=ix*dt;
					hdg.el[2*ix+1]=(ix+1)*dt;

				
			
				}


		
					//T2.el[ix]=.01-uu;
			/*	if(mdf<7)
					T2.el[ix]=uu;*/

		
				if(ix>0)
					time.el[ix]=dt+time.el[ix-1];

				if(model.motor){
		
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
		
					//Vn.el[ix]=10*model.spwmLevel;
					double tetx=(elAng-elAng0)*kme-mechAng;
					
			/*		if(tetx<0) tetx+=360;
					
					 if(tetx>=360)
						tetx=tetx-Math.round(tetx/(360))*360;*/
 
					Vn.el[ix]=Math.sin(tetx/180*PI/kme*2);
		
				
					
					tav+=model.TrqZ;

					
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
					//************************************************
							

					if(plotter){
					
						 double w1=-1e10,w2=1e-10;
						 
						double[][] XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Is[1][j];
							
					
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[0].changePlotData(0, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Is[2][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[0].changePlotData(1, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Is[3][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[0].changePlotData(2, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Is[4][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[0].changePlotData(3, XY);
						
						plot[0].setFixedBounds(0,0,50*(1+ix/50));

						double mx=1.2*Math.max(w1,w2);
						
									
						//************************************************
						
						
						  w1=-1e10;
						  w2=1e-10;
						 
						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vs[1][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(0, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vs[2][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(1, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vs[3][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(2, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vf[1][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(3, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vf[2][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(4, XY);

						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vf[3][j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[1].changePlotData(5, XY);
						
						plot[1].setFixedBounds(0,0,50*(1+ix/50));


						
						//************************************************
						
						
						  w1=-1e10;
						  w2=1e-10;
						 
						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=T.el[j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[2].changePlotData(0, XY);

						plot[2].setFixedBounds(0,0,50*(1+ix/50));
						 mx=1.2*Math.max(w1,w2);
						//plot[2].setFixedBounds(1,-mx,mx);
						
						//************************************************
						
						
						
						//************************************************
						
						
						  w1=-1e10;
						  w2=1e-10;
						 
						XY=new double[ix][2];
						for(int j=0;j<ix;j++){

							XY[j][0]=j;
							XY[j][1]=Vn.el[j];
							if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
						}
						plot[3].changePlotData(0, XY);
						
						plot[3].setFixedBounds(0,0,50*(1+ix/50));

						 mx=1.2*Math.max(w1,w2);
						 mx=400;
						//plot[3].setFixedBounds(1,-mx,mx);
						
						//************************************************
						
					}						
					
								//************************************************
				//	writeFiles=true;


					if( writeFiles /*&& i==nEnd*/){

						model.writeMesh( folder+"\\bun"+i+".txt");
						String fluxFile =  folder+"\\flux"+i+".txt";
						model.writeB(fluxFile);
						
						String forceFile =folder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);

					
					}

					if(model.saveForce){
					//	String dispFile =dispFolder+"\\force"+mangs+".txt";
						String forceFile =dispFolder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);
					}


					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));
				//	main.gui.tfX[2].setText(this.formatter.format(model.spwmLevel));
					main.gui.tfX[2].setText(this.formatter.format(dang/PI*180));

					if(model.solver.terminate){

						break;
					}

			if(model.solver.terminate) break;

			}
			

		}
			
			boolean writeInit=false;
			
			if(writeInit){
		
		
			String initfile = System.getProperty("user.dir")+"\\initxMag.txt";
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
			
			
			
			/*if(nTsteps>2)
				util.plot(time,T);*/
			
			
			util.plot(time,T);
			
			util.plot(hdg,Tdg);

		
			//	util.plotBunch(Is);

				//util.plot(de);
				T.show();
				time.show();
				
				T2.show();
				Tdg.show();
				hdg.show();
				
						
		
		
				
				if(logMode==2){
					
					String def = System.getProperty("user.dir") + "//uu"+model.defMode+".txt";
					String Trq = System.getProperty("user.dir") + "//torque.txt";
					
					String si="I";
					String sv="V";
					String s="";
					if(model.circuit) s=si;
					else
						s=sv;

					String emf1 = System.getProperty("user.dir") + "//emf//"+s+"a.txt";
					String emf2 = System.getProperty("user.dir") + "//emf//"+s+"b.txt";
					String emf3 = System.getProperty("user.dir") + "//emf//"+s+"c.txt";
					

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

	
public void transfer2DTo3D(Model model, Model m3d, boolean flux3DNeeded,boolean force3DNeeded){
		
		for(int i=1;i<=m3d.numberOfNodes;i++){
		m3d.node[i].F=null;
		m3d.node[i].Fms=null;
		}
		

		int ir=8;
		int p=0;
		int nr=1;
		int nL=6;
		int nP=1;

		int nLh=nL;
	/*	if(m3d.tag==17 ) {nL=6; nLh=nL;nP=4;}
		else
			if(m3d.tag==19 )*/ {nL=6; nLh=nL/2;}


	
		if(model.getSliceAngle()<1.6 && m3d.getSliceAngle()>6) nP=4;

		
		if(flux3DNeeded){
			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++){
					Mat R=util.rotMat2D(t*PI/2);
					for(int j=model.region[ir].getFirstEl();j<=model.region[ir].getLastEl();j++){

						Vect B=model.element[j].getB();
						if(t!=0) B=R.mul(B);
						B=B.v3();
						m3d.element[m3d.region[nr].getFirstEl()+p].setB(B);
						p++;
					}
				}
		}

		if(force3DNeeded){
		
			p=0;
	
			int he=model.nElVert;


			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++){
					
					boolean[] nc=new boolean[model.numberOfNodes+1];
					for(int j=model.region[ir].getFirstEl();j<=model.region[ir].getLastEl();j++){
						int[] vn=model.element[j].getVertNumb();

						int[] vn2=m3d.element[m3d.region[nr].getFirstEl()+p].getVertNumb();
						double kf=abs(m3d.node[vn2[vn2.length-1]].getCoord(2)-m3d.node[vn2[0]].getCoord(2))/2;

						p++;

						for(int kk=0;kk<model.nElVert;kk++){
							if(nc[vn[kk]]) continue;
							
							Vect F=model.node[vn[kk]].F;


							Vect Fms=model.node[vn[kk]].Fms;
							
							
							if(F!=null) 
							{
								F=F.v3().times(kf);
							
								if(m3d.node[vn2[kk]].F!=null){
								
								m3d.node[vn2[kk]].setF(F.add(m3d.node[vn2[kk]].F));
								}
								else
									m3d.node[vn2[kk]].setF(F);
								
								if(m3d.node[vn2[kk+he]].F!=null)
									m3d.node[vn2[kk+he]].setF(F.add(m3d.node[vn2[kk+he]].F));
									else
										m3d.node[vn2[kk+he]].setF(F);

								
						

							}

							if(Fms!=null) 
							{

								Fms=Fms.v3().times(kf);
								
								if(m3d.node[vn2[kk]].Fms!=null){
									
									m3d.node[vn2[kk]].setFms(Fms.add(m3d.node[vn2[kk]].Fms));
									}
									else
										m3d.node[vn2[kk]].setFms(Fms);
									
									if(m3d.node[vn2[kk+he]].Fms!=null)
										m3d.node[vn2[kk+he]].setFms(Fms.add(m3d.node[vn2[kk+he]].Fms));
										else
											m3d.node[vn2[kk+he]].setFms(Fms);

							}

							nc[vn[kk]]=true;
					
						
						}
					}


				}


		}
		

				
	}
		


	public void disp3Dto2D(Model model, Model m2d){

		double h0=0;

		double eps=1e-6;
		int ir=3;
		int ir2=7;
		int nt=0;



		for(int j=model.region[ir].getFirstEl();j<=model.region[ir].getLastEl();j++){
			int[] vn=model.element[j].getVertNumb();
			boolean found=false;
			for(int k=0;k<model.nElVert;k++){
				if(abs(model.node[vn[k]].getCoord(2)-h0)<eps) 
				{
					found=true;
					nt=j;
					break;
				}
			}
			if(found) break;
		}

		int kt=0;
		if(model.getElementCenter(nt).el[2]>h0) kt=1;

		int ix=0;
	
		for(int j=m2d.region[ir2].getFirstEl();j<=m2d.region[ir2].getLastEl();j++){
			int[] vn=model.element[ix+nt].getVertNumb();
					ix++;
			int[] vn2=m2d.element[j].getVertNumb();
			for(int k=0;k<m2d.nElVert;k++)
			{
				if(model.node[vn[k+kt*m2d.nElVert]].isDeformable()){
					m2d.node[vn2[k]].setDeformable(true);
					m2d.node[vn2[k]].setU(model.node[vn[k+kt*m2d.nElVert]].getU().v2());
					
				
				}
			}

		}
	}
}
