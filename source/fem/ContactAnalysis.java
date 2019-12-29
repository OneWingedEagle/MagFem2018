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
	

	boolean twice_check=false;


	private SpMat Ks=null;
//	private SpMat Ksb=null;

	private SpMatAsym node_node=null;
	private SpMatAsym Gc=null;
	private SpMatAsym Gct=null;
	private SpMatAsym Gcf=null;
	private SpMatAsym Gcft=null;
	private SpMatAsym G_stk=null;
	private SpMatAsym G_stk_prev=null;
	private SpMatAsym G_stkt=null;

	private SpMat Kc=null;
	private SpMat Kcf=null;
	private Vect lamN,gap;
	private Vect lamT,slide;;
	
	//private Vect gapDetec;;

	private Vect gap0;;
	private Vect slide0;;
	
	private Vect aug_N;
	private Vect aug_T;

	private Vect ref_stick;
	
	private Vect Fc;
	private Vect Fcf;
	
	private Vect weights;
	
	public int numContacts;
	public double[] penFactor;
	public double[] fric_coef;
	public Node[][] slaveNodes;
	public Edge[][] masterEdges;
	boolean stick[],landed_stick[];
	
	public Element[][] edgeElems;

	int[][] normalIndex;
	int[][] u_index;
	int[] numContactingNodes;
	int totalnumContactingNodes;
	
	boolean[] rmv=null;
	boolean[] contacting=null;

	Vect[][] normals;
	Vect[][] tangentials;
	
	Model model;
	double pf=1e8;
	double pft=1e8;
	double margin=0.05;
	double lamNupFactor=1;


	public Vect solve( Model model,SpMatSolver solver,int mode){


		twice_check=false;
		
		int itmax=3;
		int nr_itmax=10;
		int nLoads=5;
		int n_modifNR=1;
		
		double fp=1;
		double fr=.01;
		lamNupFactor=1;

		Vect u=new Vect(model.Ks.nRow);
int nmu=1;

Vect mm=new Vect(nmu);
 for(int contId=0;contId<numContacts;contId++)
 fric_coef[contId]=.2;
 
for(int im=0;im<nmu;im++){
	
	u.zero();
	model.setU(u);
	//	 for(int contId=0;contId<numContacts;contId++)
			// fric_coef[contId]=.1*(im);
		 
		double thickness=0.001;
		double load_scale=1./thickness;;
		
		boolean centrif=false;
		if(centrif){
			model.setNodalMass();
			double rpm=2000;
			double rps=rpm/30*PI;
			double omeag2=rps*rps;
			for(int i=1;i<=model.numberOfNodes;i++){
				Vect v=model.node[i].getCoord();
				double m=model.node[i].getNodalMass();
				Vect F=v.times(omeag2*m);
				model.node[i].setF(F);
			}
			
		}
		
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
					 else
						 util.pr("master edge ( "+node1.id+", "+node2.id+" ) belongs to no element.");

			 }
				//for(int k=0;k<masterEdges[contId].length;k++)
				//util.hshow(edgeElems[contId][k].getVertNumb());
		 }
		 
		 landed_stick=new boolean[model.numberOfNodes+1];
		 
		 stick=new boolean[model.numberOfNodes+1];
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

		Vect bU1=model.bU.add(model.getbUt(mode));

		bU1.timesVoid(load_scale);


	
		
		ref_stick=new Vect(model.Ks.nRow);

		lamN=new Vect(model.Ks.nRow);
		lamT=new Vect(model.Ks.nRow);

		aug_N=new Vect(model.Ks.nRow);
		aug_T=new Vect(model.Ks.nRow);
		
		weights=new Vect(model.Ks.nRow);//.ones(model.Ks.nRow);

		gap=new Vect(model.Ks.nRow);
		gap0=new Vect(model.Ks.nRow);
		slide0=new Vect(model.Ks.nRow);
		
		Fc=new Vect(model.Ks.nRow);
		Fcf=new Vect(model.Ks.nRow);
		solver.terminate(false);
		System.out.println(" Contact analysis....");


		

		//readContact();


	//	int nout=2601-2551+1;
		int nout=slaveNodes[0].length;
		Vect xr=new Vect(nout);
		
		boolean plot=true;//model.numberOfNodes==2832;
		if(plot)
		for(int k=0;k<xr.length;k++){
			//int n=2551+k;
			int n=slaveNodes[0][k].id;
			xr.el[k]=model.node[n].getCoord(0);
		}
		

		pf=0;
		


		 for(int contId=0;contId<numContacts;contId++)
		for(int i=0;i<slaveNodes[contId].length;i++){
			int sn=slaveNodes[contId][i].id;

			int index=model.U_unknownIndex[sn]-1;
			
	
			if(index<0) continue;
			int xind=u_index[sn][0];
			int yind=u_index[sn][1];
			double val1=0;
			double val2=0;

			if(xind!=-1)
			 val1=fp*model.Ks.row[xind].el[model.Ks.row[xind].nzLength-1];
			
			if(yind!=-1)
			 val2=fp*model.Ks.row[yind].el[model.Ks.row[yind].nzLength-1];

			if(val1>pf) pf=val1;
			if(val2>pf) pf=val2;
			int p=u_index[sn][0];
			if(p<0) 
				p=u_index[sn][1];
			weights.el[p]=penFactor[contId]*util.max(val1,val2);
		}

	
		///pf/=slaveNodes[0].length;//
		pft=fr*pf;

		util.pr("pf :"+pf);
		util.pr("pft :"+pft);

		pf=1e0;
		pft=fr*pf;
		
		normalIndex=new int[numContacts][];

		normals=new Vect[numContacts][];
		tangentials=new Vect[numContacts][];

		 for(int contId=0;contId<numContacts;contId++){
				int numSn=slaveNodes[contId].length;
				int numMed=masterEdges[contId].length;
				normalIndex[contId]=new int[numSn];
				normals[contId]=new Vect[numMed];
				tangentials[contId]=new Vect[numMed];

		 }


		int nnSize=model.numberOfNodes+1;
		node_node=new SpMatAsym(nnSize,nnSize);


		Vect dF=null;


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

			util.pr("load_iter: "+load_iter);
			
			lamN.zero();
			lamT.zero();

			double factor=(load_iter+1.)/nLoads;
			bU1=load.times(factor);

			for(int cont_iter=0; cont_iter<itmax; cont_iter++){

				util.pr("cont_iter: "+cont_iter);


				for(int nr_iter=0; nr_iter<nr_itmax; nr_iter++){	
				

					assembleContactMatrices();
					//if(num_augs_run>0)
					//Kcf.shownzA();

					
					if(Kc!=null){
					Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));
					}		
					else{
						Ks=model.Ks.deepCopy();
					}

					int ns=0;
					for(int sb=0;sb<n_modifNR;sb++){
						ns++;
					if(sb>0)		
						assembleContactMatrices();
				
						Vect dF1=calcResidual(Ks,u,bU1);

						Vect du1=solveLinear(solver,Ks.deepCopy(),dF1);
						
						double er1=dF1.norm()/bU1.norm();
						
						util.pr("nr_iter_sub: "+sb+"           nr_err_sub: "+er1);

						u=u.add(du1);
						
						model.setU(u);
						
						if(twice_check)
						checkPositiveGap(u);
						
						if(er1<1e-6) break;

						
					}


					if(ns>0){
					util.pr("ns ="+ns);
					assembleContactMatrices();
					
					if(Kc!=null){
					Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));
					}		
					else{
						Ks=model.Ks.deepCopy();
					}
					}
						
					dF=this.calcResidual(Ks, u, bU1);
				

					double er=dF.norm()/bU1.norm();
					nr_err.el[totalNRIter]=er;
					
					util.pr("nr_iter: "+nr_iter+"           nr_err: "+er);
					


					totalNRIter++;

					Vect du=solveLinear(solver,Ks.deepCopy(),dF);


					u=u.add(du);					

					model.setU(u);
					
					if(twice_check)
					checkPositiveGap(u);

					if(er<1e-3){
						break;
					}


				}
				

				gap=Gc.mul(u).add(gap0);


				err.el[cont_iter]=gap.abs().max();


				slide=G_stk.mul(u).add(slide0);
				if(G_stk_prev!=null)
				slide=slide.sub(G_stk_prev.mul(ref_stick).add(slide0));
				
				//util.plot(slide);
			///	new SpVect(slide).shownz();
				
				errf.el[cont_iter]=slide.abs().max();

	if(Kc!=null){
		for(int k=0;k<gap.length;k++){
			gap.el[k]*=weights.el[k];;
			slide.el[k]*=weights.el[k];;
		}
	//	util.plot(slide);
		if(lamN.norm()==0)
			lamN=lamN.add(gap.times(pf));
		else
			lamN=lamN.add(gap.times(lamNupFactor*pf));
			lamT=lamT.add(slide.times(pft));

			///util.pr("_________________________________________________");
			 for(int contId=0;contId<numContacts;contId++)	
			for(int i=0;i<slaveNodes[contId].length;i++){
				Node node=slaveNodes[contId][i];
				int sn=node.id;
				int p=u_index[sn][0];
				if(p<0) 
					p=u_index[sn][1];
			//	util.ph("node "+sn+ "lamT:  ");
			//	util.pr(lamT.el[p],"%4.4e");
				
			}


		}
	
	
		checkStickSlip();
	//xr.show();
				if(plot)
				for(int k=0;k<xr.length;k++){
				//	int n=2551+k;
					int n=slaveNodes[0][k].id;

					int index=model.U_unknownIndex[n]-1;
					if(index<0) continue;


					int com_index=u_index[n][0];
				//	urs[cont_iter].el[k]=-Math.abs(u.el[com_index])*1e6;
					if(com_index==-1) continue;
					urs[cont_iter].el[k]=u.el[com_index]*1e6;
					if(xr.el[k]<0)
						urs[cont_iter].el[k]*=-1;

				}

			//	if(err.el[cont_iter]<1e-6 && (errf.el[cont_iter]<1e-6)) break;
				
				num_augs_run++;
				
				aug_N=Gct.mul(lamN);;
				aug_T=Gcft.mul(lamT);;

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

		nr_distilled.show("%5.4e");
		for(int k=0;k<nr_distilled.length;k++)
			 nr_distilled.el[k]=Math.log10(nr_distilled.el[k]);
		util.plot(nr_distilled);
		//u=aug_N.add(aug_T);


		util.pr("Gap[micon] vs aug_iter");
		err.show("%5.4e");
		//util.plot(err);

		util.pr("slide[micon] vs aug_iter");
		errf.show("%5.4e");
		//util.plot(errf);


		//ur.show();
	
		if(plot){
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

mm.el[im]=u.abs().max();


}
if(nmu>1)
util.plot(mm);

	//	model.setU(new Vect(u.length));
		 for(int contId=0;contId<numContacts;contId++)		
			for(int k=0;k<masterEdges[contId].length;k++){
				
				Node node1=masterEdges[contId][k].node[0];
				Node node2=masterEdges[contId][k].node[1];
				int mn1=node1.id;
				int mn2=node2.id;

				Vect normal=normals[contId][k];//new Vect(-v12.el[1],v12.el[0]).normalized();
				if(normal==null) continue;
		
				model.node[mn1].Fms=normal.deepCopy();
				model.node[mn2].Fms=normal.deepCopy();
				//node2.u=normal.deepCopy();
			}
		model.writeNodalField( model.resultFolder+"\\master_normal.txt",2);
		//model.setU(Fc.add(Fcf.times(1)).times(-1));
		model.setU(aug_N.times(0).add(aug_T.times(1)).times(1));
		model.writeNodalField( model.resultFolder+"\\contact_force.txt",-1);

		model.setU(u);
		
		

		return u;

	}
	
	
	
	private void assembleContactMatrices(){
		
		int dof=model.Ks.nRow;
		
		obtain_node_node();

		assembleConstraintMats();
		
		util.pr("totalnumContactingNodes: "+totalnumContactingNodes);

		if(totalnumContactingNodes!=0){

			Gct=Gc.transpose(50);
			//Gct.shownzA();

			Kc=new SpMat(dof,dof); // Gct*Gc

			for(int i=0;i<Gct.nRow;i++){
				if(Gct.row[i].nzLength>0){
					SpVect spv=new SpVect(dof,100);
					
					SpVect spv1=Gct.row[i].deepCopy();
					for(int k=0;k<spv1.nzLength;k++){
						int ind=spv1.index[k];
						spv1.el[k]*=weights.el[ind];
					}

					int kx=0;
					for(int j=0;j<=i;j++){
						if(Gct.row[j].nzLength>0){
	
							double dot=spv1.dot(Gct.row[j]);
							if(dot==0) continue;
		
							spv.index[kx]=j;
							spv.el[kx++]=dot;
						}
					}
					spv.trim(kx);
					Kc.row[i]=spv.deepCopy();

				}
			}

			Kc.times(pf);


			if(lamN.norm()>0) {
			//	G_stk.times(0);
			///	Gcf.times(0);
			}

		
			Gcft=Gcf.transpose(50);
			G_stkt=G_stk.transpose(50);

		
			Kcf=new SpMat(dof,dof); // Gct*Gc

			for(int i=0;i<G_stkt.nRow;i++){
				if(G_stkt.row[i].nzLength>0){
					SpVect spv=new SpVect(dof,100);
					SpVect spv1=G_stkt.row[i].deepCopy();
					for(int k=0;k<spv1.nzLength;k++){
						int ind=spv1.index[k];
						spv1.el[k]*=weights.el[ind];
					}

					int kx=0;
					for(int j=0;j<=i;j++){
						if(G_stkt.row[j].nzLength>0){

							double dot=spv1.dot(G_stkt.row[j]);
							if(dot==0) continue;
			
							spv.index[kx]=j;
							spv.el[kx++]=dot;
						}
					}
					spv.trim(kx);
					Kcf.row[i]=spv.deepCopy();

				}
			}
			
			Kcf.times(pft);

		//	Kcf.shownzA();


		}
	}
	
	private void obtain_node_node(){
		
		
		///weights.zero();
		gap.zero();
		
		numContactingNodes=new int[numContacts];
		
		totalnumContactingNodes=0;
		//int numContacting=0;
		
		int nnSize=model.numberOfNodes+1;
		 for(int contId=0;contId<numContacts;contId++){
			 double mu=fric_coef[contId];
		for(int i=0;i<slaveNodes[contId].length;i++){
			Node node=slaveNodes[contId][i];
			int sn=node.id;
			
			node_node.row[sn]=new SpVect(nnSize);

			
			if(rmv[sn] && twice_check){
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

				double edgeLength=v12.norm();

				Vect edgeDir=v12.times(1./edgeLength);

				
				boolean node_to_node=false;
				double dot1=v1v.dot(v12);
				double dot2=v2v.dot(v12);
				if(dot1!=0) dot1/=v1v.norm()*edgeLength;
				if(dot2!=0) dot2/=v2v.norm()*edgeLength;
				
				
				if(dot1*dot2>.001) continue;
				else{
					node_to_node=true;
				}

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



				double pen=v1v.dot(normal);
			
			//	util.pr(sn+"    "+pen);
				if(pen>1e-15 &&(!contacting[sn] || !twice_check)) {
/*					stick[sn]=false;
					landed_stick[sn]=false;
					int p=u_index[sn][0];
					if(p<0) 
						p=u_index[sn][1];
					lamN.el[p]=0;
					lamT.el[p]=0;*/
				
					continue;
				}
				


				if(pen<-10*edgeLength) continue;

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
	
			///	util.ph("normal.dot(edgeDir)   "+normal.dot(edgeDir));

				node_node.row[sn]=new SpVect(nnSize,2);
				node_node.row[sn].index[0]=node1.id;
				node_node.row[sn].index[1]=node2.id;
				node_node.row[sn].el[0]=alpha;
				node_node.row[sn].el[1]=beta;
				
				if(mu!=0) {
					if(!stick[sn]&& !contacting[sn]){
					stick[sn]=true;
					landed_stick[sn]=true;
					}
				}
				
				
				contacting[sn]=true;

				
				//util.pr("node "+sn+" contacted.");
				
				
			//	gap.el[u_index[sn][0]]=pen;
				
				Vect delaDisp=v.sub(v1.times(alpha).add(v2.times(beta)));
				
				double proj=delaDisp.dot(normal);
				
			//	slide.el[u_index[sn][0]]=pen;
				
				
				Vect tang=delaDisp.sub(normal.times(proj)).normalized();
				//if(tang.norm()==0) 
					tang=edgeDir.deepCopy();
				
				//Vect tang=delaDisp.normalized();
				//if(sld1.norm2()==0)
				//	sld1=edgeDir.deepCopy();
				// tang=edgeDir.deepCopy();
				//if(sld1.norm2()==0)
				
				tangentials[contId][k]=tang.deepCopy();
				
				int p=u_index[sn][0];
				if(p<0) 
					p=u_index[sn][1];
								
				gap0.el[p]=node.getCoord().sub(node1.getCoord().times(alpha).add(node2.getCoord().times(beta))).dot(normal);
				
			//	gap.el[p]=pen;

				slide0.el[p]=-node.getCoord().sub(node1.getCoord().times(alpha).add(node2.getCoord().times(beta))).dot(tang);

				totalnumContactingNodes++;
				
				numContactingNodes[contId]++;
					break;//
			}

		}
	}
	
		 for(int contId=0;contId<numContacts;contId++)
			 util.pr((contId+1)+", num slave nodes: "+slaveNodes[contId].length+",  numContactingNodes: "+numContactingNodes[contId]);

	}

	
	private void assembleConstraintMats(){
		
		int dof=model.Ks.nRow;
		
		Gc=new SpMatAsym(dof,dof);
		Gcf=new SpMatAsym(dof,dof);
		
		if(G_stk!=null)
		G_stk_prev=G_stk.deepCopy();
		
		G_stk=new SpMatAsym(dof,dof);

		 for(int contId=0;contId<numContacts;contId++)
		for(int i=0;i<slaveNodes[contId].length;i++){
			
			Node node=slaveNodes[contId][i];
			int sn=node.id;

			if(node_node.row[sn].nzLength>0){

				int index=model.U_unknownIndex[sn]-1;
				if(index<0) continue;


				int px=u_index[sn][0];
				int py=u_index[sn][1];
				
				int p=px;
				if(p==-1) p=py;
				
				Gc.row[p]=new SpVect(dof);



				Vect normal=normals[contId][normalIndex[contId][i]];
				
				Edge edge=masterEdges[contId][normalIndex[contId][i]];

				int mn1=edge.node[0].id;
				int mn2=edge.node[1].id;


				int p1x=u_index[mn1][0];
				int p1y=u_index[mn1][1];
		
				int p2x=u_index[mn2][0];
				int p2y=u_index[mn2][1];

				double alpha=node_node.row[sn].el[0];
				double beta=node_node.row[sn].el[1];


				Gc.row[p]=new SpVect(dof,6);
								
				int kx=0;
				if(px!=-1){
				Gc.row[p].index[kx]=px;
				Gc.row[p].el[kx++]=normal.el[0];
				}
				if(p1y!=-1){
				Gc.row[p].index[kx]=py;
				Gc.row[p].el[kx++]=normal.el[1];
				}

				
//util.pr(weights.el[com_index]);

				if(p1x>=0){
					Gc.row[p].index[kx]=p1x;
					Gc.row[p].el[kx++]=-alpha*normal.el[0];
				}
				if(p1y>=0){
					Gc.row[p].index[kx]=p1y;;
					Gc.row[p].el[kx++]=-alpha*normal.el[1];
				}
	
				if(p2x>=0){
					Gc.row[p].index[kx]=p2x;
					Gc.row[p].el[kx++]=-beta*normal.el[0];
				}
				if(p2y>=0){
					Gc.row[p].index[kx]=p2y;;
					Gc.row[p].el[kx++]=-beta*normal.el[1];
				}
				Gc.row[p].sortAndTrim(kx);;
	

				//Vect tang=new Vect(-normal.el[1],normal.el[0]);
				
				Vect tang=tangentials[contId][normalIndex[contId][i]];
	
				if(fric_coef[contId]==0){
					tang.zero();
				}
				Gcf.row[p]=new SpVect(dof,6);
				kx=0;
				if(px!=-1){
					Gcf.row[p].index[kx]=px;
					Gcf.row[p].el[kx++]=tang.el[0];
					}
				if(p1y!=-1){
					Gcf.row[p].index[kx]=py;
					Gcf.row[p].el[kx++]=tang.el[1];
				}

				
				if(p1x>=0){
					Gcf.row[p].index[kx]=p1x;
					Gcf.row[p].el[kx++]=-alpha*tang.el[0];
				}
				if(p1y>=0){
					Gcf.row[p].index[kx]=p1y;;
					Gcf.row[p].el[kx++]=-alpha*tang.el[1];
				}
	
				if(p2x>=0){
					Gcf.row[p].index[kx]=p2x;
					Gcf.row[p].el[kx++]=-beta*tang.el[0];
				}
				if(p2y>=0){
					Gcf.row[p].index[kx]=p2y;;
					Gcf.row[p].el[kx++]=-beta*tang.el[1];
				}
				

				if(stick[sn]){
				G_stk.row[p]=Gcf.row[p].deepCopy();
			//	util.pr("node -_______________> "+sn+" stick "+stick[sn]);
				if(landed_stick[sn]){
					Vect u=model.node[sn].getU();
					if(px>=0)
					ref_stick.el[px]=u.el[0];
					if(py>=0)
						ref_stick.el[py]=u.el[1];
					
					 u=model.node[mn1].getU();
					if(p1x>=0)
					ref_stick.el[p1x]=u.el[0];
					if(p1y>=0)
						ref_stick.el[p1y]=u.el[1];
					
					 u=model.node[mn2].getU();
						if(p2x>=0)
						ref_stick.el[p2x]=u.el[0];
						if(p2y>=0)
							ref_stick.el[p2y]=u.el[1];
					
					
					landed_stick[sn]=false;

				}
				}else{
					G_stk.row[p]=Gcf.row[p].times(0);
				}
				
		
			}
		}

	///G_stk.shownzA();

			
	}
	
	
	private void checkPositiveGap(Vect u){


		gap=Gc.mul(u).add(gap0);
		
	
		for(int i=1;i<=model.numberOfNodes;i++){
			//for(int k=0;k<model.dim;k++){
				int p=u_index[i][0];
				if(p<0)  p=u_index[i][1];;
				if(p<0) continue;
			//	if(gap.el[p]!=0) util.pr(i+"    "+gap.el[p]);
				if(twice_check){
				if(gap.el[p]>0 && !rmv[i]){
				//	util.pr(p+"  "+gap.el[p]);
				//	gap.el[p]=0;
				//	slide.el[p]=0;
					contacting[i]=false;
					lamN.el[p]=0;
					lamT.el[p]=0;
					rmv[i]=true;
					}
				}else{
					
					if(gap.el[p]>0){
			
						lamN.el[p]=0;
						lamT.el[p]=0;
					}
				}
			
	
		}





	}
	

	private void checkStickSlip(){

		 for(int contId=0;contId<numContacts;contId++){
			 double mu=this.fric_coef[contId];
			 if(mu==0) continue;
				for(int i=0;i<slaveNodes[contId].length;i++){
					
					Node node=slaveNodes[contId][i];
					int sn=node.id;

			

						int index=model.U_unknownIndex[sn]-1;
						if(index<0) continue;


						int px=u_index[sn][0];
						int py=u_index[sn][1];
						
						int p=px;
						if(p==-1) p=py;
						if(node_node.row[sn].nzLength>0){
						if(lamN.el[p]<=0){
							double abs_lamT=Math.abs(lamT.el[p]);
							double muFn=mu*Math.abs(lamN.el[p]);
							if(abs_lamT>muFn*(1+margin)){
								if(lamT.el[p]>0)
									lamT.el[p]=muFn;
								else
									lamT.el[p]=-muFn;
								stick[sn]=false;
								landed_stick[sn]=false;
							}else{
								if(!stick[sn]){
								stick[sn]=true;

								landed_stick[sn]=true;
								}
							}
						}else{
							stick[sn]=false;
							landed_stick[sn]=false;
							lamT.el[p]=0;
						}
					}else{
						stick[sn]=false;
						landed_stick[sn]=false;
						lamT.el[p]=0;
					}
				}

				for(int i=0;i<slaveNodes[contId].length;i++){
					Node node=slaveNodes[contId][i];
					int sn=node.id;
					if(contacting[sn]){
						int px=u_index[sn][0];
						int py=u_index[sn][1];
						
						int p=px;
						if(p==-1) p=py;
					//	double abs_lamT=Math.abs(lamT.el[p]);
					//	double muFn=mu*Math.abs(lamN.el[p]);
					//	util.pr("lamT.el[p] "+lamT.el[p] +" lamN.el[p] "+lamN.el[p]);
				/*	util.ph("node "+sn+ "lamT:  ");
					util.pr(lamT.el[p],"%4.4e");
					util.ph("node "+sn+ "lamN:  ");
					util.pr(lamN.el[p],"%4.4e");*/
						util.pr("node "+sn+" stick "+stick[sn]);
					}
					else
						util.pr("node "+sn+" free ");
				}
		 }


	}

	private Vect solveLinear(SpMatSolver solver, SpMat Ks, Vect dF){
		
		Vect du=new Vect(dF.length);

		
		model.Ci=Ks.scale(dF);


			//if(model.Ls==null)
			model.Ls=Ks.ichol();


			if(dF.abs().max()!=0){
				//if(model.xp==null){
					du=solver.ICCG(Ks,model.Ls, dF,model.errCGmax,model.iterMax);
				//}
			//	else{
				//	du=solver.ICCG(Ks,model.Ls, dF,model.errCGmax,model.iterMax,model.xp);
					//du=model.solver.err0ICCG(Ks,model.Ls, dF,model.errCGmax*1e-3,model.iterMax,model.xp);	

				//}
			}
			else{
				util.pr("Solution is zero!");
				du=new Vect(Ks.nRow);
			}

			//model.xp=du.deepCopy();

			du.timesVoid(model.Ci);
		
			return du;
	}
	
	private Vect calcResidual(SpMat Ks,Vect u, Vect b){
		
		Vect Fint=model.Ks.smul(u);

		Vect dF=b.sub(Fint);

		Fc=Kc.smul(u);
		Fcf=Kcf.smul(u);
		
		dF=dF.sub(Fc).sub(Fcf); //
			
		dF=dF.sub(aug_N).sub(aug_T);
		
			return dF;

	}
	
	
	public static void main(String[] args){

		new Main();
	}

}


