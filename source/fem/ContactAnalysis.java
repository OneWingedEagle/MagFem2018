package fem;


import static java.lang.Math.PI;

import main.Main;
import math.Complex;
import math.DFT;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatAsym;
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
public class ContactAnalysis {

	private SpMat Ks=null;
	private SpMat Ksb=null;

	private SpMatAsym node_node=null;
	private SpMatAsym Gc=null;
	private SpMatAsym Gct=null;
	private SpMatAsym Gcf=null;
	private SpMatAsym Gcft=null;

	private SpMat Kc=null;
	private SpMat Kcf=null;
	private Vect lamN;
	private Vect lamT;;

	private Vect aug_N;
	private Vect aug_T;

	private Vect Fc;
	private Vect Fcf;
	
	public int numContacts;
	public double[] penFactor;
	public double[] fric_coef;
	public Node[][] slaveNodes;
	public Edge[][] masterEdges;
	
	public Element[][] edgeElems;
	//int[] slaveNodes;
	//int[][] masterEdges;
	int[][] normalIndex;
	int[][] u_index;
	
	boolean[] rmv=null;
	boolean[] contacting=null;

	Vect[][] normals;

	
	Model model;
	double pf=1e8;
	double pft=1e8;

	boolean check_double=false;


	public Vect solve( Model model,SpMatSolver solver,int mode){


		this.model=model;
		
		int[][] ne=new int[model.numberOfNodes+1][20];
		int[] nz=new int[model.numberOfNodes+1];


		 for(int i=1;i<=model.numberOfElements;i++){
			 int[] vn=model.element[i].getVertNumb();
			 for(int j=0;j<vn.length;j++){
				 ne[vn[j]][nz[vn[j]]++]=i;
					 
				 }
			 }
		
		edgeElems=new Element[numContacts][];
		 for(int contId=0;contId<numContacts;contId++){
			 edgeElems[contId]=new Element[masterEdges[contId].length];
		 
				for(int k=0;k<masterEdges[contId].length;k++){
					
					Node node1=masterEdges[contId][k].node[0];
					Node node2=masterEdges[contId][k].node[1];
					int ie=0;
					 for(int j=0;j<nz[node1.id];j++){
						 for(int p=0;p<nz[node2.id];p++){
						
							 if(ne[node1.id][j]==ne[node2.id][p]){
								 ie=ne[node1.id][j];
								 break;
							 }
							 if(ie>0) break;
						 }
				}
					 
					 if(ie>0) 
						 edgeElems[contId][k]=model.element[ie];
					 

			 }
		 }
		 
		
		contacting=new boolean[model.numberOfNodes+1];
		 rmv=new boolean[model.numberOfNodes+1];

		u_index=new int[model.numberOfNodes+1][model.dim];
		for(int i=1;i<=model.numberOfNodes;i++)
			for(int k=0;k<model.dim;k++)
				u_index[i][k]=-1;
			
		int ix=0;
		for(int i=1;i<=model.numberOfUnknownU;i++){
			int nodeNumb=model.unknownUnumber[i];

			if(model.node[nodeNumb].isDeformable()){
				for(int k=0;k<model.dim;k++){
					if(!model.node[ nodeNumb].is_U_known(k)){
						u_index[nodeNumb][k]=ix;
						ix++;
					}
				}

			}
		}
	//	Mat LU=null;

		int skip=0;

		MatSolver matSolver=new MatSolver();

		Vect bU1=model.bU.add(model.getbUt(mode));

		bU1.timesVoid(1e4);


		Vect u=new Vect(model.Ks.nRow);

		lamN=new Vect(model.Ks.nRow);
		lamT=new Vect(model.Ks.nRow);

		aug_N=new Vect(model.Ks.nRow);
		aug_T=new Vect(model.Ks.nRow);


		Fc=new Vect(model.Ks.nRow);
		Fcf=new Vect(model.Ks.nRow);
		solver.terminate(false);
		System.out.println(" Contact analysis....");


		//readContact();


		int nout=2601-2551+1;
		Vect xr=new Vect(nout);
		
		boolean punch=true;
		if(punch)
		for(int k=0;k<xr.length;k++){
			int n=2551+k;
			xr.el[k]=model.node[n].getCoord(0);
		}
		
		

		pf=0;
		

		 for(int contId=0;contId<numContacts;contId++)
		for(int i=0;i<slaveNodes[contId].length;i++){
			int sn=slaveNodes[contId][i].id;

			int index=model.U_unknownIndex[sn]-1;
			if(index<0) continue;
			int xind=2*index;
			int yind=2*index+1;
			double val1=model.Ks.row[xind].el[model.Ks.row[xind].nzLength-1];
			double val2=model.Ks.row[yind].el[model.Ks.row[yind].nzLength-1];

			if(val1>pf) pf=val1;
			if(val2>pf) pf=val2;
		}

		pf*=1;
		//pf/=slaveNodes[0].length;//
		pft=.1*pf;

		util.pr("pf :"+pf);
		util.pr("pft :"+pft);
		
		
		normalIndex=new int[numContacts][];

		normals=new Vect[numContacts][];

		 for(int contId=0;contId<numContacts;contId++){
				int numSn=slaveNodes[contId].length;
				int numMed=masterEdges[contId].length;
				normalIndex[contId]=new int[numSn];

				normals[contId]=new Vect[numMed];

		 }
		int dof=model.Ks.nRow;


		int nnSize=model.numberOfNodes+1;
		node_node=new SpMatAsym(nnSize,nnSize);
		//	node_node.lower=false;


		//node_node.shownzA();

		//node_node.matForm().plot();

		int dim=model.dim;

		Gc=new SpMatAsym(dof,dof);
		Gcf=new SpMatAsym(dof,dof);

		int itmax=1;
		int nr_itmax=10;
		int nLoads=1;

		Vect[] urs=new Vect[itmax];
		for(int i=0;i<itmax;i++)
			urs[i]=new Vect(xr.length);

		int num_augs_run=0;
		Vect err=new Vect(itmax);

		Vect errf=new Vect(itmax);

		Vect nr_err=new Vect(nLoads*itmax*nr_itmax);

		int totalNRIter=0;
		Vect load=bU1.deepCopy();
		for(int load_iter=0; load_iter<nLoads; load_iter++){

			double factor=(load_iter+1.)/nLoads;
			bU1=load.times(factor);

			for(int cont_iter=0; cont_iter<itmax; cont_iter++){

				util.pr("cont_iter: "+cont_iter);

				for(int nr_iter=0; nr_iter<nr_itmax; nr_iter++){	
				

					int numContacting=obtain_node_node();

					util.pr("numContacting:   "+numContacting);

					assembleConstraintMats();
					
					

					if(numContacting!=0){

						//Gc.shownzA();
						//	Gct=new SpMatAsym(Gc.matForm().transp());
						Gct=Gc.transpose(50);

						//Gct.shownzA();

						Kc=new SpMat(dof,dof); // Gct*Gc

						for(int i=0;i<Gct.nRow;i++){
							if(Gct.row[i].nzLength>0){
								SpVect spv=new SpVect(dof,100);

								int kx=0;
								for(int j=0;j<=i;j++){
									if(Gct.row[j].nzLength>0){

										double dot=Gct.row[i].dot(Gct.row[j]);
										if(dot==0) continue;
										//	util.pr(i+" ---- "+j);
										spv.index[kx]=j;
										spv.el[kx++]=dot;
									}
								}
								spv.trim(kx);
								Kc.row[i]=spv.deepCopy();

							}
						}

						Kc.times(pf);


						Gcft=Gcf.transpose(50);
						

						Kcf=new SpMat(dof,dof); // Gct*Gc

						for(int i=0;i<Gcft.nRow;i++){
							if(Gcft.row[i].nzLength>0){
								SpVect spv=new SpVect(dof,100);

								int kx=0;
								for(int j=0;j<=i;j++){
									if(Gcft.row[j].nzLength>0){

										double dot=Gcft.row[i].dot(Gcft.row[j]);
										if(dot==0) continue;
										//	util.pr(i+" ---- "+j);
										spv.index[kx]=j;
										spv.el[kx++]=dot;
									}
								}
								spv.trim(kx);
								Kcf.row[i]=spv.deepCopy();

							}
						}

						Kcf.times(pft);
						//	Kc.shownzA();
						//	Kc.plot();

						//Mat M=model.Ks.matForm();

						//	M=M.add(Kc.matForm());

						//Ks= new SpMat(M);


						//if(Ksb==null){
						//if(nr_iter%skip==0){
							Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));

						//	LU=Ks.matForm();
						//	LU.lu();



						//}/*else{
		// Ks=Ksb.deepCopy();//model.Ks.addGeneral(Kc.addGeneral(Kcf));

		//}
						Ksb=Ks.deepCopy();

						//}
						//	 else
						//	Ks=Ksb.deepCopy();


						// Ks.plot();


					}
					else
						Ks=model.Ks.deepCopy();

					//Ks.plot();

					///	util.wait(10000);
					//model.Ks.matForm(true).show();

					//	model.Ks.showcl();

					boolean solve=true;


					if(solve){


						Vect Fint=Ks.smul(u);


						Vect dF=bU1.sub(Fint);

						
						boolean opened=false;

						if(Kc!=null){


							Fc=Kc.smul(u).times(1);
							Fcf=Kcf.smul(u).times(1);
						
							Vect g=Gc.mul(u);
		
							for(int i=0;i<g.length;i++){
								if(g.el[i]>0){
									//util.pr(i+"  "+g.el[i]);
									//lamN.el[i]=0;
								//	lamT.el[i]=0;
								//	opened=true;
								}
							}
						
							if(check_double)
							for(int i=1;i<=model.numberOfNodes;i++){
								//for(int k=0;k<model.dim;k++){
									int p=u_index[i][0];
									if(p<0) continue;
									if(g.el[p]>0 && !rmv[i]){
										contacting[i]=false;
										lamN.el[p]=0;
										lamT.el[p]=0;
										rmv[i]=true;
										//opened=true;
										}
								
						
							}
							
					

							aug_N=Gct.mul(lamN);;
							aug_T=Gcft.mul(lamT);;

							dF.sub(Fc).sub(Fcf);

							dF=dF.sub(aug_N).sub(aug_T);



						}




						double er=dF.norm()/bU1.norm();
						nr_err.el[totalNRIter]=er;
						
						util.pr("nr_iter: "+nr_iter+"           nr_err: "+er);

						if(!opened && er<1e-3) break;

						totalNRIter++;

						Vect du=null;

						if(skip>0 && nr_iter%skip!=0){
						//	du=matSolver.solvelu(LU,dF);
						}else{

							model.Ci=Ks.scale(dF);


							//if(model.Ls==null)
							model.Ls=Ks.ichol();


							if(dF.abs().max()!=0){
								if(model.xp==null){
									du=solver.ICCG(Ks,model.Ls, dF,model.errCGmax,model.iterMax);
								}
								else{
									//	u=solver.ICCG(model.Ks,model.Ls, bU1,2e-3,model.iterMax,model.xp);
									du=model.solver.err0ICCG(Ks,model.Ls, dF,model.errCGmax*1e-3,model.iterMax,model.xp);	

								}
							}
							else{
								util.pr("Solution is zero!");
								du=new Vect(Ks.nRow);
							}

							model.xp=du.deepCopy();

							du.timesVoid(model.Ci);
						}

						u=u.add(du);

						model.setU(u);


					}


				}

				Vect gap=Gc.mul(u);


				err.el[cont_iter]=gap.abs().max();


				Vect sld=Gcf.mul(u);
				errf.el[cont_iter]=sld.abs().max();

	if(Kc!=null){

			lamN=lamN.add(gap.times(pf));
			lamT=lamT.add(sld.times(pft));


		}
				if(punch)
				for(int k=0;k<xr.length;k++){
					int n=2551+k;

					int index=model.U_unknownIndex[n]-1;
					if(index<0) continue;


					int com_index=u_index[n][0];
					urs[cont_iter].el[k]=-Math.abs(u.el[com_index])*1e6;

				}
				num_augs_run++;

				if(err.el[cont_iter]<1e-8 && (errf.el[cont_iter]<1e-8)) break;
			}

		}


		util.pr("NR error");

		int nn=0;
		for(int k=0;k<nr_err.length;k++)
			if(nr_err.el[k]>0) nn++;

		Vect nr_distilled=new Vect(nn);

		int kx=0;
		for(int k=0;k<nr_err.length;k++)
			if(nr_err.el[k]>0) nr_distilled.el[kx++]=(nr_err.el[k]);

		nr_distilled.show();
		for(int k=0;k<nr_distilled.length;k++)
			 nr_distilled.el[k]=Math.log10(nr_distilled.el[k]);
		util.plot(nr_distilled);
		//u=aug_N.add(aug_T);


		util.pr("Gap[micon] vs aug_iter");
		err.times(1e6).show();
		//util.plot(err);

		util.pr("slide[micon] vs aug_iter");
		errf.times(1e6).show();
		//util.plot(errf);


		//ur.show();
	
		if(punch){
		double[][] data=new double[xr.length][itmax+1];
		for(int k=0;k<xr.length;k++){

			data[k][0]=xr.el[k];

			for(int i=0;i<itmax;i++)
				data[k][i+1]= urs[i].el[k];

		}
	//	util.show(data);
		util.plotBunch(data);
		//util.plot(x,ur);
		}

		model.setU(new Vect(u.length));
		 for(int contId=0;contId<numContacts;contId++)		
			for(int k=0;k<masterEdges[contId].length;k++){
				
				Node node1=masterEdges[contId][k].node[0];
				Node node2=masterEdges[contId][k].node[1];
				int mn1=node1.id;
				int mn2=node2.id;

				Vect normal=normals[contId][k];//new Vect(-v12.el[1],v12.el[0]).normalized();
				if(normal==null) continue;
				//normal.hshow();
				model.node[mn1].u=normal.deepCopy();
				model.node[mn2].u=normal.deepCopy();
				//node2.u=normal.deepCopy();
			}
		model.writeNodalField( model.resultFolder+"\\master_normal.txt",-1);
		model.setU(Fc.add(Fcf).times(-1));
		model.writeNodalField( model.resultFolder+"\\contact_force.txt",-1);

		model.setU(u);

		return u;

	}
	
	
	
	private int obtain_node_node(){
		
		
		int numContacting=0;
		
		int nnSize=model.numberOfNodes+1;
		 for(int contId=0;contId<numContacts;contId++)
		for(int i=0;i<slaveNodes[contId].length;i++){
			
			Node node=slaveNodes[contId][i];
			int sn=node.id;
			
			if(rmv[sn] && check_double){
				rmv[sn]=false;
				continue;
			}
			
			Vect v=model.node[sn].getCoord().add(model.node[sn].u);
			
			for(int k=0;k<masterEdges[contId].length;k++){
				
				Node node1=masterEdges[contId][k].node[0];
				Node node2=masterEdges[contId][k].node[1];
				int mn1=node1.id;
				int mn2=node2.id;
				Vect v1=node1.getCoord().add(node1.u);
				Vect v2=node2.getCoord().add(node2.u);
				Vect v12=v2.sub(v1);
				//v12.hshow();
				Vect v1v=v.sub(v1);
				Vect v2v=v.sub(v2);

				Vect edgeDir=v12.normalized();

				double dot1=v1v.dot(v12);
				double dot2=v2v.dot(v12);
				if(dot1*dot2>0) continue;

				Element elem=edgeElems[contId][k];
				 int[] vn=elem.getVertNumb();
				 for(int j=0;j<vn.length;j++){
					if(vn[j]!=mn1 && vn[j]!=mn2){
						Node node3=model.node[vn[j]];
	
						Vect v3=node3.getCoord().add(node3.u);	
						Vect v13=v3.sub(v1);
						
						Vect cross1=edgeDir.v3().cross(v13.v3());
						Vect cross2=edgeDir.v3().cross(cross1.v3());
						normals[contId][k]=cross2.normalized().v2();
						break;
					}
				 }
						 
					 

				Vect normal=normals[contId][k];//new Vect(-v12.el[1],v12.el[0]).normalized();
				//if(v1.el[1]<.0201) 
				//normals[contId][k].timesVoid(-1);
				//normal.hshow();

				double pen=v1v.dot(normal);

				if(pen>0e-8 && (!contacting[sn] && !check_double)) {

					continue;
				}
				

				double edgeLength=v12.norm();

				if(pen<-2*edgeLength) continue;

				normalIndex[contId][i]=k;

				double beta=0;
				double alpha=0;
				double v1n=v1v.norm();

				if(v1n==0){
					beta=0;
					alpha=1;
				}
				///if(sn>=120) util.pr(sn+"   "+v1n);
				else{
					beta=v1v.dot(edgeDir)/edgeLength;			
					alpha=1-beta;
				}
	

				node_node.row[sn]=new SpVect(nnSize,2);
				node_node.row[sn].index[0]=node1.id;
				node_node.row[sn].index[1]=node2.id;
				node_node.row[sn].el[0]=alpha;
				node_node.row[sn].el[1]=beta;
				
				contacting[sn]=true;

				numContacting++;
				//	break;//
			}

		}
		
		return numContacting;
	}

	
	private void assembleConstraintMats(){
		
		int dim=model.dim;
		int nnSize=model.numberOfNodes+1;
		
		 for(int contId=0;contId<numContacts;contId++)
		for(int i=0;i<slaveNodes[contId].length;i++){
			
			Node node=slaveNodes[contId][i];
			int sn=node.id;

			if(node_node.row[sn].nzLength>0){

				int index=model.U_unknownIndex[sn]-1;
				if(index<0) continue;


				int com_index=u_index[sn][0];
				


				Vect normal=normals[contId][normalIndex[contId][i]];
				
				Edge edge=masterEdges[contId][normalIndex[contId][i]];

				int mn1=edge.node[0].id;
				int mn2=edge.node[1].id;



				int com_index1=u_index[mn1][0];
				int com_index2=u_index[mn2][0];

				double alpha=node_node.row[sn].el[0];
				double beta=node_node.row[sn].el[1];


				Gc.row[com_index]=new SpVect(nnSize,6);
				int kx=0;
				Gc.row[com_index].index[kx]=u_index[sn][0];
				Gc.row[com_index].el[kx++]=normal.el[0];
				Gc.row[com_index].index[kx]=u_index[sn][1];
				Gc.row[com_index].el[kx++]=normal.el[1];

				if(com_index1>=0){
					Gc.row[com_index].index[kx]=u_index[mn1][0];;
					Gc.row[com_index].el[kx++]=-alpha*normal.el[0];
					Gc.row[com_index].index[kx]=u_index[mn1][1];;
					Gc.row[com_index].el[kx++]=-alpha*normal.el[1];
				}
				if(com_index2>=0){
					Gc.row[com_index].index[kx]=u_index[mn2][0];
					Gc.row[com_index].el[kx++]=-beta*normal.el[0];
					Gc.row[com_index].index[kx]=u_index[mn2][1];;
					Gc.row[com_index].el[kx++]=-beta*normal.el[1];
				}
				Gc.row[com_index].sortAndTrim(kx);;




				Vect tang=new Vect(-normal.el[1],normal.el[0]);
				
				if(fric_coef[contId]==0){
					tang.zero();
				}
				Gcf.row[com_index]=new SpVect(nnSize,6);
				kx=0;
				Gcf.row[com_index].index[kx]=u_index[sn][0];
				Gcf.row[com_index].el[kx++]=tang.el[0];
				Gcf.row[com_index].index[kx]=u_index[sn][1];
				Gcf.row[com_index].el[kx++]=tang.el[1];

				if(com_index1>=0){
					Gcf.row[com_index].index[kx]=u_index[mn1][0];
					Gcf.row[com_index].el[kx++]=-alpha*tang.el[0];
					Gcf.row[com_index].index[kx]=u_index[mn1][1];
					Gcf.row[com_index].el[kx++]=-alpha*tang.el[1];
				}
				if(com_index2>=0){
					Gcf.row[com_index].index[kx]=u_index[mn2][0];
					Gcf.row[com_index].el[kx++]=-beta*tang.el[0];
					Gcf.row[com_index].index[kx]=u_index[mn2][1];
					Gcf.row[com_index].el[kx++]=-beta*tang.el[1];
				}

				Gcf.row[com_index].sortAndTrim(kx);;
				
	


				//if(mn1==408) util.pr((com_index1+1)+" ---  "+ alpha+"  "+beta);
				//	if(mn2==408) util.pr((com_index2+1)+" ---  "+ alpha+"  "+beta);
				//util.pr(sn+"   "+index);
				//Gc.row[index].showr();
				//Gc.row[index].shownz();
			}
		}
	}
	
	
/*	private void readContact(){
		int numSn=13;
		int numMed=8;
		slaveNodes=new int[numSn];
		for(int k=0;k<slaveNodes.length;k++)
			slaveNodes[k]=359+k;

		masterEdges=new int[numMed][2];
		for(int k=0;k<masterEdges.length;k++){
			masterEdges[k][0]=406+k;
			masterEdges[k][1]=masterEdges[k][0]+1;
			//util.hshow(masterEdges[k]);
		}

		boolean thick=false;

		if(thick){
			numSn=21;
			numMed=10;
			slaveNodes=new int[numSn];
			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=106+k;

			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=127+k;
				masterEdges[k][1]=masterEdges[k][0]+1;
				//util.hshow(masterEdges[k]);	
			}
		}

		boolean ipm=true;
		if(ipm){
			numSn=73;
			numMed=82;

			// numSn=3;
			// numMed=3;
			slaveNodes=new int[numSn];
						 int kx1=0;
			 slaveNodes[kx1++]=1298;;;
			 slaveNodes[kx1++]=2408;
			 slaveNodes[kx1++]=1272;

			String file="D:\\JavaWorks\\FEM problems\\structural\\2D-contact\\ipm\\cont.txt";
			int[][] data=model.loader.loadIntArray(file, numSn, 1,1);

			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=data[k][0];

			data=model.loader.loadIntArray(file, numMed, 2,numSn+2);

			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=data[k][0];
				masterEdges[k][1]=data[k][1];
				//util.hshow(masterEdges[k]);	
			}	

		}
		boolean punch=false;
		if(punch){
			numSn=22;
			numMed=20;

			// numSn=3;
			// numMed=3;
			slaveNodes=new int[numSn];
						 int kx1=0;
			 slaveNodes[kx1++]=1298;
			 slaveNodes[kx1++]=2408;
			 slaveNodes[kx1++]=1272;

			String file="D:\\JavaWorks\\FEM problems\\structural\\2D-contact\\punch\\cont.txt";
			int[][] data=model.loader.loadIntArray(file, numSn, 1,1);

			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=data[k][0];

			data=model.loader.loadIntArray(file, numMed, 2,numSn+2);

			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=data[k][0];
				masterEdges[k][1]=data[k][1];
				//util.hshow(masterEdges[k]);	
			}	

		}
	}*/
}


