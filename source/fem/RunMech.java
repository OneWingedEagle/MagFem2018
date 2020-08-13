package fem;


import static java.lang.Math.PI;
import static java.lang.Math.abs;

import main.Main;
import math.Eigen;
import math.Mat;
import math.SpVect;
import math.Vect;
import math.util;

public class RunMech {

	public int[] nn=new int[1];
	boolean terminate=false;
	//SpMatSolver solver=new SpMatSolver();

	public RunMech(){}

	public static void main(String[] args){
		new Main();
	}


	public void runMech(Model model,Main main) {


		if(model.loadDisp){
			runLoadDisp(model,main);
			return;


		}
		
/*		
		for(int n=1;n<=model.numberOfNodes;n++){
			if(model.node[n].getCoord(2)>.000999){
				double r=model.node[n].getCoord().v2().norm();
				if(Math.abs(r-.005)<1e6){
					if(model.node[n].getCoord(1)<.0001) r/=2;
				//if(r<.0011 || r>.00999 || model.node[n].getCoord(0)<.0001 || model.node[n].getCoord(1)<.0001) r/=2;
				util.pr(n+"\t"+"0	.0001  "+(-r*100000));
				}
			}
		}
*/
		//*************************
		//*************************
		//*************************
		//*************************

		/*			model.region[2].thermal=true;
		model.region[2].deltaT=-100;

		model.region[2].thermalCoef=1.1e-5;

		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(model.region[ir].thermal)
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					model.element[i].setHasThermal(true);
					model.element[i].setDeltaT(model.region[ir].deltaT);

				}

		}	*/


		//*************************
		//*************************
		//*************************
		//*************************

		/*		model.setThermalForce();

		for(int i=1;i<=model.numberOfNodes;i++){
			if(model.node[i].Fms!=null)

			model.node[i].F=model.node[i].Fms.deepCopy();
		}
		 */


		//=======================================
		boolean layeredData=false;
		int[] reg8nodes=new int[1];

		if(layeredData){
		model.mapnr=new int[1];
		if(model.dim==3 && model.transfer2DTo3D){
			String m2dBun= System.getProperty("user.dir") + "\\mot4th2DFine.txt";

			model.m2d=new Model(m2dBun);

			model.mapnr=new int[model.m2d.numberOfNodes+1];
			reg8nodes=model.m2d.getRegNodes(8);
			for(int u=0;u<reg8nodes.length;u++)
				model.mapnr[reg8nodes[u]]=u;
			model.forceLamin=new Vect[6][model.m2d.region[8].getNumbElements()];

			model.m2d.coordCode=1;
		}
		//=======================================
		}


		double tstart=System.currentTimeMillis();

		if(model.modal)	model.modalAnalysis(20, 1e-5,true);

		util.pr("mech integ mode: "+model.timeIntegMode);



		if(model.timeIntegMode==-1){
			ModalDecomp mdcp=new ModalDecomp();
			mdcp.setVibModalDecomp(model, main);

		}
		else if(model.timeIntegMode==-2){
			HarmonicSuperposition hsp=new HarmonicSuperposition();
			hsp.setVibFreqDomain(model, main);
		}
		else if(model.timeIntegMode==-3){
			POD pod=new POD();
			pod.setVibPOD(model, main);
		}

		else{

			Vect u=new Vect();

			int N=model.nTsteps;

			Vect T=new Vect(N);
			model.writeMesh( model.eddyFolder+"\\bun.txt");

			int ix=0;

			//	model.setNodalMass();

			if(model.loadPrevMech){

				String initfile = System.getProperty("user.dir")+"\\initxMech.txt";
				model.loader.loadPrevMech(model,initfile);

			}

			Vect UU=new Vect(model.nTsteps);

			boolean snap=false;

			int D=10;
			Mat W=new Mat();

			if(snap)
				W=new Mat(model.numberOfUnknownUcomp,D);

			if(!model.modal)
				for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){



					String file=model.forceFile[0];

					main.gui.tfX[1].setText(Integer.toString(i)+"/"+(model.nEnd));

					if(model.dim==3 && model.transfer2DTo3D){

						if(layeredData){

							String[] files=new String[model.forceLamin.length];

							for(int k=0;k<files.length;k++){

								int m=i%1800;

								files[k]="C:\\Works\\proj8\\forcesMotMSz"+k+"\\force"+m+".txt";
								//files[k]=System.getProperty("user.dir") + "\\forcesMotMSz"+k+"\\force"+m+".txt";
								model.m2d.loadNodalField(files[k], 1);
								for(int j=0;j<reg8nodes.length;j++){
									model.forceLamin[k][j]=model.m2d.node[reg8nodes[j]].F;
								}
							}
							
							model.transfer2DTo3DLayered( file,false,true);
						}
						else{
							model.m2d.loadNodalField(file, 1);
							model.transfer2DTo3D( file,false,true);
						}



					}
					else{



						model.loadNodalField(file,1);

				if(false){
					for (int p = 1; p <= model.numberOfNodes; p++) {
							Vect v = model.node[p].getCoord();
							if(v.el[1]>.299){
								Vect F = new Vect(0,-1,0);						
								
								model.node[p].setF(F);
							}else if(v.el[1]<.0001){
							//	Vect F = new Vect(0,1,0);						
								
							//	model.node[p].setF(F);
							}
						}
					}
/*						double factor=(ix+1.)*1/model.nTsteps;
						for (int p = 1; p <= model.numberOfNodes; p++) {
						//	Vect v = model.node[p].getCoord();
							Vect F = model.node[p].F;
							if(F!=null){
								if(Math.abs(F.el[1])<1e-6)
								{
									model.node[p].setF(F.times(factor));
								}else{
						//	model.node[p].setF(F.times(1e3*factor));
							}
							}
							
						}*/
						




					}


					//	model.writeNodalField( model.resultFolder+"\\force"+i+".txt",1);


					System.gc();
					
					if(model.timeIntegMode==0){
						u=model.setDeformation(ix);


					}
					else
						u=model.setVibration(ix);


					if(snap && ix<D)
						W.setCol(u, ix);

					int nx0=2352;
					nx0=1558;int cmp=1; //react
					//nx0=15;
					//	nx0=24697; cmp=0;// motor half


					int nx=Math.min(nx0,model.numberOfNodes);

					//	UU.el[ix]=model.node[nx].getU(cmp);


					ix++;

		


					if(model.saveForce|| ix==model.nTsteps)
						model.writeNodalField( model.eddyFolder+"\\force_out"+i+".txt",1);

					if(model.saveDisp || ix==model.nTsteps)
						model.writeNodalField( model.eddyFolder+"\\disp"+i+".txt",-1);
					
					//if( ix==model.nTsteps-1)
					//	model.writeMesh( model.eddyFolder+"\\bun"+i+".txt",true);


					if(model.solver.terminate) break;

					if(model.saveStress){
						model.setStress();
						model.writeStress( model.eddyFolder+"\\stress"+i+".txt");
					}


				}

			if(snap){
				String solutions=System.getProperty("user.dir") +"\\solutions.txt";
				model.writer.writeArray(W.el, solutions);
			}

			model.writeData(model.resultFolder+"\\data_out.txt");

/*			UU.timesVoid(1e9);
			util.plot(UU);

			UU.show();*/

			/*		for(int i=1;i<=model.numberOfElements;i++)
			model.element[i].setDeformable(true);
			 */

			/*		model.setStress();
		model.writeStress(model.resultFolder+"\\stress.txt");*/
			/*	util.p		double tend=System.currentTimeMillis();
		double compTime=(tend-tstart)/1000;

		util.pr("Computation time: "+compTime+" seconds");
	lot(T);
		T.show();*/




			boolean writeInit=false;

			if(model.up==null) writeInit=false;

			if(writeInit){


				String initfile = System.getProperty("user.dir")+"\\initxMech.txt";
				double[][] arr=new double[model.numberOfUnknownUcomp][4];

				Vect v1=model.up.deepCopy();

				Vect v2=model.upp.deepCopy();
				Vect v3=model.ud;
				Vect v4=model.udd;


				for(int j=0;j<v1.length;j++){
					if(v1!=null)
						arr[j][0]=v1.el[j];

					if(v2!=null)
						arr[j][1]=v2.el[j];

					if(v3!=null)
						arr[j][2]=v3.el[j];

					if(v4!=null)
						arr[j][3]=v4.el[j];
				}


				model.writer.writeArray(arr,initfile);


			}



		}

		double tend=System.currentTimeMillis();
		double compTime=(tend-tstart)/1000;

		util.pr("Computation time: "+compTime+" seconds");

		try {
			Thread.sleep(100);
		} catch (InterruptedException exception) {
			// TODO Auto-generated catch-block stub.
			exception.printStackTrace();
		}


	}

	public void runLoadDisp(Model model,Main main) {




		double tstart=System.currentTimeMillis();


		util.pr("Loading dispacements");


		int N=model.nTsteps;

		Vect T=new Vect(N);
		model.writeMesh( model.eddyFolder+"\\bun.txt");




		int ix=0;



		for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){



			String file=model.forceFile[ix];
			file =System.getProperty("user.dir") + "\\dispMagHalf\\disp"+i+".txt";


			model.loadNodalField(file,-1);

			main.gui.tfX[1].setText(Integer.toString(i)+"/"+(model.nEnd));


			if(model.saveStress){
				model.setStress();
				model.writeStress( model.resultFolder+"\\stress"+i+".txt");
			}



			int nx0=2352;
			nx0=1558;int cmp=1; //react
			//nx0=15;
			nx0=24697; cmp=0;// motor half

			int nx=Math.min(nx0,model.numberOfNodes);

			T.el[ix]=model.node[nx].getU(cmp);

			ix++;


			if(model.solver.terminate) break;


		}


/*		T=T.times(1e9);
		util.plot(T);
		T.show();
*/




		double tend=System.currentTimeMillis();
		double compTime=(tend-tstart)/1000;

		util.pr("Computation time: "+compTime+" seconds");

		try {
			Thread.sleep(10);
		} catch (InterruptedException exception) {
			// TODO Auto-generated catch-block stub.
			exception.printStackTrace();
		}


	}









}

