package fem;


import static java.lang.Math.PI;

import java.io.File;

import main.Main;
import math.Complex;
import math.DFT;

import math.SpMatComp;

import math.Vect;
import math.VectComp;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class HarmonicSuperposition {


	public HarmonicSuperposition(){}

	public void setVibFreqDomain(Model model, Main main){
		
		
				util.pr("mech integ mode: "+model.timeIntegMode);
				
				util.pr("dt:"+model.dt);
				
				int N=model.nTsteps;
				
				Vect T=new Vect(N);
				
				model.setMechBC();
				
				int ix=0;
				
				int mode=1;
				int Nn=model.numberOfUnknownUcomp;
				
				
			
				
				int Lt=model.nTsteps;
		
				double w0=2*PI/model.dt/Lt;
	 
				
				double a1=model.rayAlpha;
				double a2=model.rayBeta;


				int D=Math.min(20,Lt/2+1);
			
				File fftfolder = new File(System.getProperty("user.dir") + "\\alfft");
			
		boolean fftReady=false;		
				
			if(!fftReady){
			
			Vect[] Fs=new Vect[Nn];
			for(int i=0;i<Fs.length;i++)
				Fs[i]=new Vect(model.nTsteps);
			
				for(int i=	model.nBegin;i<=model.nEnd;	i+=model.nInc){

			
					String file=model.forceFile[ix];

					main.gui.tfX[1].setText(Integer.toString(i)+"/"+(	model.nEnd));
					
					if(model.dim==3 && model.transfer2DTo3D){


			model.transfer2DTo3D( file,false,true);



					}
					else
					model.loadNodalField(file,1);
					
					for(int n=1;n<=0*model.numberOfNodes;n++){
						
						Vect v=model.node[n].getCoord();
				
					/*	if(v.el[1]>.3999 &&v.el[0]>.3999){
							int m=(ix%180);
							model.node[n].F=new Vect(100,100,0).times(Math.cos(2*PI*m/180)*(1.0+m/180.0));
						}*/
						
						
						if( v.el[0]>-.4999){
							int m=(ix%180);
							model.node[n].F=new Vect(0,0,-1);//.times(Math.cos(2*PI*m/180)+0*Math.cos(0*PI*m/180));
						}
						
						
	
							}
	
		

					Vect bF=model.getbUt(mode);
										
					for(int j=0;j<Fs.length;j++)
						Fs[j].el[ix]=bF.el[j];
				ix++;

				if(model.solver.terminate) break;

				}

		
			
				Complex[][] FF=new Complex[Nn][D];
				for(int i=0;i<Nn;i++){
					Complex[] X=DFT.fft(Fs[i].el);
					for(int j=0;j<D;j++){
					FF[i][j]=X[j];
					}
				}
				
	
				if(fftfolder.exists()){
					util.deleteDir(fftfolder);
				}
				else
					fftfolder.mkdir();
		
				for(int j=0;j<D;j++){
					Complex[] vc=new Complex[Nn];
					for(int i=0;i<Nn;i++)
						vc[i]=FF[i][j];
				model.writer.writeFFT(vc,fftfolder + "\\fftc"+j+".txt");
				}			
				
				
				FF=null;
				 Fs=null;
				System.gc();
			
	}
				
				 model.setStiffMat(true);

				
				Complex[][] U=new Complex[Nn][D];
		
					for(int i=0;i<Nn;i++)
						for(int j=0;j<D;j++)
							U[i][j]=new Complex();

	

				Complex[] F=new Complex[Nn];
				
				for(int j=0;j<D;j++){
					
					F=model.loader.loadFFT(fftfolder + "\\fftc"+j+".txt");
					
					util.pr("harmonic "+j);
				
					VectComp b=new VectComp(Nn);
					for(int i=0;i<Nn;i++){
						b.el[i]=F[i].deepCopy();
					}
					VectComp xc;
					
					if(b.norm()<1e-8){
						xc=new VectComp(Nn);
						continue;
					}
		
					
							
				
					SpMatComp Ks=new SpMatComp(model.Ks.addNew(model.Ms.timesNew(-Math.pow(j*w0,2))),model.Ms.timesNew(a1).addNew(model.Ks.timesNew(a2)).timesNew(j*w0));
				
			
		
					Ks.setSymHerm(1);
					
					model.Ci=Ks.scale(b);
					
					SpMatComp Ls=Ks.ichol(1.1);
					Ls.setSymHerm(0);

				
				

						xc=model.solver.COICCG(Ks,Ls,b,model.errCGmax,model.iterMax,new VectComp(Nn),1,true);
						//xc=model.solver.COCG(Ks,Ls,b,model.errMax,model.iterMax,new VectComp(Nn),1,true);
			
					
					xc.timesVoid(model.Ci);	
	

					
					for(int i=0;i<Nn;i++){
					U[i][j]=xc.el[i].deepCopy();
				
					}
					
			
				
				}
			

				
				 ix=0;
					
				 int[][] umap=new int[model.numberOfNodes+1][model.dim];
				 
				for(int i=1;i<=model.numberOfUnknownU;i++){
										
					if(model.node[i].isDeformable()){
						for(int k=0;k<model.dim;k++){
							if(!model.node[ i].is_U_known(k)){
								umap[i][k]=ix;
								ix++;
							}
						}

					}
				}
				
				int cmp=0;

				int nx0=2352;
				nx0=1558; cmp=1; //react
					//nx0=24697;  cmp=0;// motor half
					//nx0=26;
			
					
					int nx=Math.min(nx0,model.numberOfNodes);
				
				int nodeNumb=nx;
				
				int pp=umap[nodeNumb][cmp];
				
			//	pp=0;
				
				Complex[] Utmp=new Complex[Lt];
				for(int j=0;j<Lt;j++)
					if(j<D) 
						Utmp[j]=U[pp][j].deepCopy();
					else
					Utmp[j]=new Complex();

				
				if(D>1)
				for(int j=0;j<D-1;j++)
					Utmp[Lt-1-j]=U[pp][j+1].conj();

				
			

				Complex[] Y=DFT.ifft(Utmp);
				for(int i=0;i<Lt;i++)
					T.el[i]=Y[i].re*1e9;
	
				util.plot(T);
				T.show();
						
				
				if(model.saveDisp){
	/*				if(D>1)
						for(int j=0;j<D-1;j++)
							for(int i=0;i<Nn;i++)
								U[i][Lt-1-j]=U[i][j+1].conj();*/
							
					
						for(int i=0;i<Nn;i++){
						
							Utmp=new Complex[Lt];
							for(int j=0;j<Lt;j++)
								if(j<D) 
									Utmp[j]=U[i][j].deepCopy();
								else
								Utmp[j]=new Complex();

							
							if(D>1)
							for(int j=0;j<D-1;j++)
								Utmp[Lt-1-j]=U[i][j+1].conj();
							
							U[i]=DFT.ifft(Utmp);
						}

						for(int j=0;j<Lt;j++){
							Vect u=new Vect(Nn);
							for(int i=0;i<Nn;i++)
								u.el[i]=U[i][j].re;
							model.setU(u);
							model.writeNodalField( model.resultFolder+"\\disp"+j*model.nInc+".txt",-1);

						}

				
				}
	
				


	}
	
	



	



}


