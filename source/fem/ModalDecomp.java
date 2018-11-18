package fem;


import static java.lang.Math.PI;

import main.Main;
import math.Complex;
import math.DFT;
import math.Mat;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class ModalDecomp {


	public ModalDecomp(){}

	public void setVibModalDecomp(Model model, Main main){


		util.pr("mech integ mode: "+model.timeIntegMode);

		util.pr("dt:"+model.dt);


		int N=model.nTsteps;

		Vect T=new Vect(N);
		//	model.writeMesh( model.resultFolder+"\\bun.txt");

		model.setMechBC();

		//String eigfile=System.getProperty("user.dir") +"\\eigsReact10\\eigVects.txt";
		//String eigfile=System.getProperty("user.dir") +"\\eigsRingHalf100\\eigVects.txt";
		String eigfile=System.getProperty("user.dir") +"\\eigs40Half\\eigVects.txt";
		
		int ix=0;

		int mode=1;
		int L=model.numberOfUnknownUcomp;

		util.pr(model.unknownUnumber[1]+" NNN  ");


		int D=10;
		
		util.pr(" Using "+D+" lowest modes.");


		int Lt=model.nTsteps;



		double w0=2*PI/model.dt/Lt;
		

	//	util.pr(1./model.dt/Lt);

		Mat Q=new Mat(L,D);

		if(model.modal){
			model.modalAnalysis(D, 1e-5,false);
			Q=model.eigVects;
		}

		else{
			Q=new Mat(model.loader.loadArrays(L, D, eigfile));
		}

		model.setStiffMat(true);




		Vect vv=Q.transp().mul(model.Ms.smul(Q)).diagonal();


		for(int i=0;i<Q.nCol;i++){
			Vect q=Q.getColVect(i).times(1.0/Math.sqrt(vv.el[i]));
			Q.setCol(q, i);

		}
		
		
		Mat Qt=Q.transp();

		Vect omega2=Q.transp().mul(model.Ks.smul(Q)).diagonal();

		
		omega2.sqrt().times(.5/PI).show();

		Vect[] Fqs=new Vect[omega2.length];
		for(int i=0;i<Fqs.length;i++)
			Fqs[i]=new Vect(Lt);
		
		util.pr(model.nTsteps);
		
/*		RunMech rm=new RunMech();
		
		Model m2d=new Model();
		if(model.dim==3 && model.transfer2DTo3D){
		String m2dBun= System.getProperty("user.dir") + "\\mot4th2DFine.txt";

		 m2d=new Model(m2dBun);
		m2d.coordCode=1;
		}*/

		for(int i=	model.nBegin;i<=model.nEnd;	i+=model.nInc){


			String file=model.forceFile[ix];

			main.gui.tfX[1].setText(Integer.toString(i)+"/"+(	model.nEnd));
			
			if(model.dim==3 && model.transfer2DTo3D){

				//m2d.loadNodalField(file,1);
				model.transfer2DTo3D(file, false,true);
				
/*				for(int n=1;n<=model.numberOfNodes;n++)
					if(model.node[n].F!=null)
						model.node[n].F=new Vect(1,0,0);
				if(i<85){
					for(int n=1;n<=model.numberOfNodes;n++)
						if(model.node[n].F!=null)
							model.node[n].F=new Vect(1,0,0);
					}
						else{
								for(int n=1;n<=model.numberOfNodes;n++)
									if(model.node[n].F!=null)
										model.node[n].F=new Vect(0,0,0);
						}*/


			}
			else
				model.loadNodalField(file,1);
				for(int n=1;n<=0*model.numberOfNodes;n++){
					
					Vect v=model.node[n].getCoord();
			
					
					if( v.el[0]>-.4999){
						int m=(ix%180);
						model.node[n].F=new Vect(0,0,-10).times(Math.cos(8*PI*m/180)+Math.cos(0*PI*m/180));
					}
	
						}

		

			Vect bF=model.getbUt(mode);
			

			Vect Fq=Qt.mul(bF);

			for(int j=0;j<Fqs.length;j++)
				Fqs[j].el[ix]=Fq.el[j];
			

			ix++;


			if(model.solver.terminate) break;

		}


		Complex[][] FqF=new Complex[D][Lt];
		for(int i=0;i<D;i++){
			Complex[] X=DFT.fft(Fqs[i].el);
			for(int j=0;j<Lt;j++){
				FqF[i][j]=X[j];
					
			}
		}
		
		
		
		double a1=model.rayAlpha;
		double a2=model.rayBeta;


		Complex[][] UqF=new Complex[D][Lt];
		for(int i=0;i<D;i++)
			for(int j=0;j<Lt;j++)
				UqF[i][j]=new Complex();
		
		int Lx=80;

		for(int i=0;i<D;i++){
			for(int j=0;j<Lx;j++){
				UqF[i][j]=FqF[i][j].times(new Complex(omega2.el[i]-Math.pow(j*w0,2),j*w0*(a1+a2*omega2.el[i])).inv());
				
			}
		}

		int N1=Q.nRow;
		
		Complex[][] U=new Complex[N1][Lt];

	//	int pp=0;
		

		
		for(int i=0;i<N1;i++){
		//	if(i!=pp) continue;
			for(int j=0;j<Lx;j++){
				
				
				Complex s=new Complex();
		
				for(int n=0;n<D;n++)					
					s=s.add(UqF[n][j].times(Q.el[i][n]));
			
				
				U[i][j]=s.deepCopy();
				
				
				if(j>0 )
					U[i][Lt-j]=s.conj();

				
			}
		}
		
		
		//================================
		
		for(int i=0;i<N1;i++){
				for(int j=0;j<Lt;j++)
					if(U[i][j]==null)
						U[i][j]=new Complex();
		}
		//================================
		
		

		Complex[] X=new Complex[Lt];

	
		 ix=0;
	
		 int[][] umap=new int[model.numberOfNodes][model.dim];
		 
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
			nx0=1558; cmp=1; /// reactor
		nx0=24697; cmp=0;// motor half
			
			int nx=Math.min(nx0,model.numberOfNodes);

		int nodeNumb=nx;
		
		int pp=umap[nodeNumb][cmp];

		//pp=0;
		
		for(int i=0;i<Lt;i++){
			if(U[pp][i]==null)
				X[i]=new Complex();
			else{
				X[i]=U[pp][i].deepCopy();
				
			}
	
		}
		
		
	

		Complex[] Y=DFT.ifft(X);
		for(int i=0;i<Lt;i++)
			T.el[i]=Y[i].re*1e9;

		util.plot(T);
		T.show();
		
		if(model.saveDisp){
			
			for(int i=0;i<N1;i++){
				U[i]=DFT.ifft(U[i]);
			}

			for(int j=0;j<Lt;j++){
				Vect u=new Vect(N1);
				for(int i=0;i<N1;i++)
					u.el[i]=U[i][j].re;
				model.setU(u);
				model.writeNodalField( model.resultFolder+"\\disp"+j*model.nInc+".txt",-1);

			}

	
	}




	}









}


