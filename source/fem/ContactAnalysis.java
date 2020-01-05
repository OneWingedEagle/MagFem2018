package fem;


import static java.lang.Math.PI;

import main.Main;
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
	boolean twice_check_mnr=false;

	private SpMat Ks=null;

	private SpMatAsym node_node=null;
	private SpMatAsym Gc=null;
	private SpMatAsym Gct=null;
	private SpMatAsym Gcf=null;
	private SpMatAsym Gcft=null;
	private SpMatAsym G_stk=null;
	private SpMatAsym G_stkt=null;
	private SpMatAsym G_stk_pr=null;


	private SpMat Kc=null;
	private SpMat Kcf=null;
	private Vect lamN,gap;
	private Vect lamT,slide;//,slide_prev;

//private Vect gapDetec;;

private Vect gap0;;
//private Vect slide0;;

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
double lamTupFactor=1;
boolean allow_sep=true;
boolean allow_sep_mnr=true;

public int method=0;


public Vect solve( Model model,SpMatSolver solver,int mode){


	if(method==1)
		return solveSingle(model,solver,mode);
	
	MatSolver slv=new MatSolver();


	allow_sep=true;
	twice_check=false;

	allow_sep_mnr=true;
	twice_check_mnr=false;

	int itmax=3;
	int nr_itmax=8;
	int nLoads=1;
	int n_modifNR=0;

	double fp=10;
	double fr=.01;
	lamNupFactor=0;
	lamTupFactor=1;
	

	boolean axi=true;

	double loadFactor0=1000;//*2*PI;//1000/2/PI;//*1.57000;//*23;//4.95;//*2.439;//100*.8;
	if(axi) loadFactor0*=2*PI;
	Vect u=new Vect(model.Ks.nRow);
	int nmu=1;
	Vect mus=new Vect(nmu);
	Vect mm=new Vect(nmu);
	
	if(axi)
	for(int i=1;i<=model.numberOfNodes;i++){
		Vect F=model.node[i].F;
		double r=model.node[i].getCoord(0);

		if(F!=null)
			model.node[i].F=model.node[i].F.times(r);
	
	}
	Mat KK=null;
	boolean direct=false;

	for(int im=0;im<nmu;im++){
		u.zero();
		model.setU(u);
		mus.el[im]=.1+.1*(im);
		for(int contId=0;contId<numContacts;contId++) fric_coef[contId]=mus.el[im];

		double thickness=1;//0.001;
		double load_scale=1./thickness*loadFactor0;;

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
		slide=new Vect(model.Ks.nRow);
		gap0=new Vect(model.Ks.nRow);
		//	slide_prev=new Vect(model.Ks.nRow);

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
				
				weights.el[p]=1;//penFactor[contId]*util.max(val1,val2);
				
	
			}


		///pf/=slaveNodes[0].length;//
		pft=fr*pf;

		util.pr("pf :"+pf);
		util.pr("pft :"+pft);

		//pf=1e0;
		//pft=fr*pf;

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
		Vect aug_disp_err=new Vect(itmax);
		Vect errf=new Vect(itmax);

		Vect nr_err=new Vect(nLoads*itmax*(nr_itmax+n_modifNR));
		int nr_it[]=new int[nLoads*itmax*(nr_itmax+n_modifNR)];

		int totalNRIter=0;
		Vect load=bU1.deepCopy();
		for(int load_iter=0; load_iter<nLoads; load_iter++){

			util.pr("load_iter: "+load_iter);

			//	lamN.zero();
			//	lamT.zero();

			double factor=(load_iter+1.)/nLoads;
			bU1=load.times(factor);

			for(int cont_iter=0; cont_iter<itmax; cont_iter++){

				util.pr("cont_iter: "+cont_iter);
				
				//	lamN.zero();
				//	lamT.zero();
				//	slide.zero();
				//	gap.zero();
					//aug_N.zero();
					//aug_T.zero();
				
				//	u.zero();
				//	model.setU(u);
					
				Vect uaug=u.deepCopy();

		
				for(int nr_iter=0; nr_iter<nr_itmax; nr_iter++){	

					assembleContactMatrices(twice_check);


					if(Kc!=null){
		
						Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));
					}		
					else{
						Ks=model.Ks.deepCopy();
					}

					if(direct&& n_modifNR>0){
						KK=Ks.matForm();
						KK.lu();
					}

					int ns=0;

					double er1=1;//
					
					Vect uc=u.deepCopy();
					
					if(nr_iter>0)
					for(int sb=0;sb<n_modifNR;sb++){
						ns++;
						if(sb>0)		
							assembleContactMatrices(twice_check_mnr);

						Vect dF1=calcResidual(Ks,u,bU1);

						Vect du1=null;
						
						if(direct){
							du1=slv.solvelu(KK, dF1);	
						}else
						du1=solveLinear(solver,Ks.deepCopy(),dF1);

						er1=dF1.norm()/bU1.norm();

						util.pr("nr_iter_sub: "+sb+"           nr_err_sub: "+er1);

						u=u.add(du1);

						model.setU(u);

						if(allow_sep_mnr)
							checkPositiveGap(u);
						
						//nr_err.el[totalNRIter]=er1;
						//nr_it[totalNRIter]=totalNRIter;

						//totalNRIter++;

						if(er1<1e-2&& sb>2){
							break;
						}


					}
					


					if(er1>10){
						u=uc.deepCopy();
					}



					util.pr("ns ="+ns);

					if(n_modifNR>0){
						assembleContactMatrices(twice_check);

						if(Kc!=null && n_modifNR>0){
							Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));
						}		
						else{
							Ks=model.Ks.deepCopy();
						}
					}

					dF=this.calcResidual(Ks, u, bU1);


					double er=dF.norm()/bU1.norm();
					nr_err.el[totalNRIter]=er;
					nr_it[totalNRIter]=totalNRIter;

					util.pr("nr_iter: "+nr_iter+"           nr_err: "+er);



					totalNRIter++;
					

					Vect du=solveLinear(solver,Ks.deepCopy(),dF);



					u=u.add(du);					

					model.setU(u);

					if(allow_sep)
						checkPositiveGap(u);
					

					if(er<1e-3 /*&& nr_iter>4*/){
						break;
					}


				}

				
			double dip_err=u.sub(uaug).norm()/u.norm();
				
				aug_disp_err.el[cont_iter]=dip_err;
				
				if(dip_err<1e-6) break;

				//	gap=Gc.mul(u).add(gap0);


				err.el[cont_iter]=gap.abs().max();


				Vect pgap=gap.deepCopy();
				Vect pslide=slide.deepCopy();//.sub(slide_prev);
				errf.el[cont_iter]=slide.abs().max();



					for(int k=0;k<gap.length;k++){
						pgap.el[k]=weights.el[k]*gap.el[k];;
						pslide.el[k]=weights.el[k]*slide.el[k];;
					}

					if(lamN.norm()==0)
						lamN=lamN.add(pgap.times(pf));
					else
						lamN=lamN.times(lamNupFactor).add(pgap.times(pf));

				
						lamT=lamT.times(lamTupFactor).add(pslide.times(pft));


			checkStickSlip();
			
			lamN.timesVoid(lamNupFactor);
			lamT.timesVoid(lamTupFactor);

			
				//xr.show();
				if(plot)
					for(int k=0;k<xr.length;k++){
						//	int n=2551+k;
						int n=slaveNodes[0][k].id;

						int index=model.U_unknownIndex[n]-1;
						if(index<0) continue;


						int com_index=u_index[n][1];
						//	urs[cont_iter].el[k]=-Math.abs(u.el[com_index])*1e6;
						if(com_index==-1) continue;
						urs[cont_iter].el[k]=u.el[com_index]*1e3;
						if(xr.el[k]<0)
							urs[cont_iter].el[k]*=-1;

					}

				//	if(err.el[cont_iter]<1e-6 && (errf.el[cont_iter]<1e-6)) break;
				
	


				num_augs_run++;

			}
		}

		util.pr("NR error");

		int nn=0;
		for(int k=0;k<nr_err.length;k++)
			if(nr_err.el[k]>0) nn++;

		Vect nr_distilled=new Vect(nn);

		int nr_it_distilled[]=new int[nn];
		int kx=0;
		for(int k=0;k<nr_err.length;k++){
			if(nr_err.el[k]>0){
				nr_it_distilled[kx]=nr_it[k];

				nr_distilled.el[kx++]=(nr_err.el[k]);
			}
		}

		//	nr_distilled.show("%5.4e");
		for(int k=0;k<nr_distilled.length;k++){
			System.out.format(" %d\t%6.4e\n",nr_it_distilled[k], nr_distilled.el[k]);
			nr_distilled.el[k]=Math.log10(nr_distilled.el[k]);
		}
		util.plot(nr_distilled);
		//u=aug_N.add(aug_T);


		
		util.pr("aug_disp changes vs aug_iter");
		aug_disp_err.show("%5.4e");

		
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
		util.plot(mus,mm);
	mm.show("%5.4e");

	//	s6.show("%5.4e");
	//	util.plot(s6);

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
	
	
	for(int i=1;i<=model.numberOfNodes;i++){
		Vect F=model.node[i].Fms;
		if(F!=null)
			model.node[i].Fms=null;
	}
	
	for(int contId=0;contId<numContacts;contId++)		
		for(int k=0;k<slaveNodes[contId].length;k++){

			Node node=slaveNodes[contId][k];
			
			int ind=	normalIndex[contId][k];
			
			if(ind<0) continue;

			Vect normal=normals[contId][ind];//new Vect(-v12.el[1],v12.el[0]).normalized();
			if(normal==null) continue;

			model.node[node.id].Fms=normal.times(-1);
			//node2.u=normal.deepCopy();
		}
	model.writeNodalField( model.resultFolder+"\\slave_normal.txt",2);
	model.setU(Fc.add(Fcf.times(1)).times(-1));
	//model.setU(aug_N.times(1).add(aug_T.times(1)).times(-1));
	model.writeNodalField( model.resultFolder+"\\contact_force.txt",-1);

	model.setU(u);



	return u;

}


public Vect solveSingle( Model model,SpMatSolver solver,int mode){


	MatSolver slv=new MatSolver();
	allow_sep=false;
	twice_check=true;

	allow_sep_mnr=false;
	twice_check_mnr=false;

	int itmax=1;
	int nr_itmax=10;
	int nLoads=1;
	int n_modifNR=0;

	double fp=1;
	double fr=.1;
	lamNupFactor=1;

	double loadFactor0=5;//*23;//4.95;//*2.439;//100*.8;
	Vect u=new Vect(model.Ks.nRow);


	for(int contId=0;contId<numContacts;contId++)

		fric_coef[contId]=0.5;
	int nmu=1;
	Vect mus=new Vect(nmu);
	Vect mm=new Vect(nmu);
	//for(int contId=0;contId<numContacts;contId++)
	// fric_coef[contId]=.1;


	for(int im=0;im<nmu;im++){
		u.zero();
		model.setU(u);
		mus.el[im]=.5+.1*(im);
		//	for(int contId=0;contId<numContacts;contId++) fric_coef[contId]=mus.el[im];

		double thickness=1;//0.001;
		double load_scale=1./thickness*loadFactor0;;

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
		slide=new Vect(model.Ks.nRow);
		gap0=new Vect(model.Ks.nRow);
		//	slide_prev=new Vect(model.Ks.nRow);

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
				weights.el[p]=penFactor[contId];//*util.max(val1,val2);
			}


		///pf/=slaveNodes[0].length;//
		pft=fr*pf;

		util.pr("pf :"+pf);
		util.pr("pft :"+pft);

		//pf=1e0;
		//pft=fr*pf;

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

		Mat KK=null;
		boolean direct=true;

		Vect[] urs=new Vect[itmax];
		for(int i=0;i<itmax;i++)
			urs[i]=new Vect(xr.length);

		int num_augs_run=0;
		Vect err=new Vect(itmax);

		Vect errf=new Vect(itmax);

		Vect nr_err=new Vect(nLoads*itmax*nr_itmax);
		int nr_it[]=new int[nLoads*itmax*nr_itmax];

		int totalNRIter=0;
		Vect load=bU1.deepCopy();
		for(int load_iter=0; load_iter<nLoads; load_iter++){

			util.pr("load_iter: "+load_iter);

			//	lamN.zero();
			//	lamT.zero();

			double factor=(load_iter+1.)/nLoads;
			bU1=load.times(factor);

			for(int cont_iter=0; cont_iter<itmax; cont_iter++){

				util.pr("cont_iter: "+cont_iter);


				for(int nr_iter=0; nr_iter<nr_itmax; nr_iter++){	


					assembleContactMatrices(twice_check);



					if(Kc!=null){

						if(direct){
							Ks=model.Ks.addGeneral(Kc);
							KK=Ks.matForm();
							for(int i=0;i<KK.nRow;i++)
								for(int j=i+1;j<KK.nRow;j++)
									KK.el[i][j]=	KK.el[j][i];

							Mat Kf=Kcf.matForm();
							for(int i=0;i<Kf.nRow;i++)
								for(int j=i+1;j<Kf.nRow;j++){
									if(G_stk.row[i].norm()>0)
										Kf.el[i][j]=	Kf.el[j][i];
									else
										Kf.el[i][j]=	0.00*Kf.el[j][i];
								}
							KK=KK.add(Kf);
							//	KK.show();

							KK.lu();
						}else{
							Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));

						}

					}		
					else{
						Ks=model.Ks.deepCopy();
					}


					int ns=0;

					double er1=1;//
					Vect uc=u.deepCopy();
					for(int sb=0;sb<n_modifNR;sb++){
						ns++;
						if(sb>0)		
							assembleContactMatrices(twice_check_mnr);

						Vect dF1=calcResidual(Ks,u,bU1);

						Vect du1=null;
						if(direct){


							du1=slv.solvelu(KK, dF1);
						}
						else
							du1=solveLinear(solver,Ks.deepCopy(),dF1);

						er1=dF1.norm()/bU1.norm();

						util.pr("nr_iter_sub: "+sb+"           nr_err_sub: "+er1);

						u=u.add(du1);

						model.setU(u);

						if(allow_sep_mnr)
							checkPositiveGap(u);

						if(er1<1e-3&& sb>2){
							break;
						}


					}

					if(er1>10){
						u=uc.deepCopy();
					}

					if(ns>0){
						util.pr("ns ="+ns);
						assembleContactMatrices(twice_check);

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
					nr_it[totalNRIter]=nr_iter;

					util.pr("nr_iter: "+nr_iter+"           nr_err: "+er);



					totalNRIter++;

					Vect du=null;

					if(direct){

						if(nr_itmax==1)
							du=new Vect(dF.length);
						else
							du=slv.solvelu(KK, dF);
					}
					else
						du=solveLinear(solver,Ks.deepCopy(),dF);


					u=u.add(du);					

					util.pr("u_max ===============>  "+u.abs().max());


					//	model.setU(dF);
					//	model.writeNodalField( model.resultFolder+"\\nrdisp\\disp"+nr_iter+".txt",-1);

					model.setU(u);

					//	model.writeNodalField( model.resultFolder+"\\nrdisp\\disp"+nr_iter+".txt",-1);


					if(allow_sep)
						checkPositiveGap(u);

					if(er<1e-3 && nr_iter>0){
						break;
					}

					Vect pgap=gap.deepCopy();
					Vect pslide=slide.deepCopy();//.sub(slide_prev);
					errf.el[cont_iter]=slide.abs().max();


					for(int k=0;k<gap.length;k++){
						pgap.el[k]=weights.el[k]*gap.el[k];;
						pslide.el[k]=weights.el[k]*slide.el[k];;

					}




					if(lamN.norm()==0)
						lamN=lamN.times(0).add(pgap.times(pf));
					else
						lamN=lamN.times(0).add(pgap.times(lamNupFactor*pf));

					//	if(nr_iter>4) lamN=lamN.ones(lamN.length).times(-1e4);
					lamT=lamT.times(0).add(pslide.times(pft));
					//	new SpVect(pslide).shownzA();
					if(nr_iter>2){ 
						//	lamT=lamT.times(-1);
						//margin=1000.1;
						//	fric_coef[0]=100;
						//	lamN=lamT.times(2.5);
					}
					//	if(nr_iter<4)
					new SpVect(lamT).shownzA();
					checkStickSlip();
					new SpVect(lamT).shownzA();

					//	util.pr("======  new SpVect(lamT).shownzA() ======");
					//Gcft.shownzA();
					//new SpVect(lamT).shownzA();
					//new SpVect(pslide).shownzA();
					//	s.show();
					aug_T=Gcft.mul(lamT);
					//new SpVect(aug_T).shownzA();
					//	aug_N=Gct.mul(lamN);

					//	new SpVect(aug_T).shownzA();
					//new SpVect(dF).shownzA();
					util.pr("ooooooooooooooooooooooooooooo");
				}





				err.el[cont_iter]=gap.abs().max();

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

			}

		}


		util.pr("NR error");

		int nn=0;
		for(int k=0;k<nr_err.length;k++)
			if(nr_err.el[k]>0) nn++;

		Vect nr_distilled=new Vect(nn);

		int nr_it_distilled[]=new int[nn];
		int kx=0;
		for(int k=0;k<nr_err.length;k++){
			if(nr_err.el[k]>0){
				nr_it_distilled[kx]=nr_it[k];

				nr_distilled.el[kx++]=(nr_err.el[k]);
			}
		}

		//	nr_distilled.show("%5.4e");
		for(int k=0;k<nr_distilled.length;k++){
			System.out.format(" %d\t%6.4e\n",nr_it_distilled[k], nr_distilled.el[k]);
			nr_distilled.el[k]=Math.log10(nr_distilled.el[k]);
		}
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
		util.plot(mus,mm);
	mm.show("%5.4e");

	util.pr("uMax ============================: "+u.abs().max()+"\n");
	//	s6.show("%5.4e");
	//	util.plot(s6);

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
	Vect FF=Fc.add(Fcf.times(1)).times(-1);
	FF=FF.add(aug_N.times(1).add(aug_T.times(1)).times(-1));
	model.setU(FF);
	model.writeNodalField( model.resultFolder+"\\contact_force.txt",-1);

	model.setU(u);



	return u;

}


private void assembleContactMatrices(boolean twchk){

	int dof=model.Ks.nRow;

	obtain_node_node(twchk);

	assembleConstraintMats();

	util.pr("totalnumContactingNodes: "+totalnumContactingNodes);

	if(totalnumContactingNodes!=0){

		Gct=Gc.transpose(50);

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

		Gcft=Gcf.transpose(50);
		G_stkt=G_stk.transpose(50);

		//G_stk.shownzA();
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




	}
}

private void obtain_node_node(boolean twchk){


	//slide_prev.zero();
	gap.zero();
	gap0.zero();

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


			if(rmv[sn] && twchk){

				rmv[sn]=false;
				continue;
			}


			Vect u=node.u;
			Vect v=node.getCoord().add(u);

			for(int k=0;k<masterEdges[contId].length;k++){

				Node node1=masterEdges[contId][k].node[0];
				Node node2=masterEdges[contId][k].node[1];
				int mn1=node1.id;
				int mn2=node2.id;
				Vect u1=node1.u;
				Vect u2=node2.u;
				Vect v1=node1.getCoord().add(u1);
				Vect v2=node2.getCoord().add(u2);
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


				if(dot1*dot2>.001) {
					contacting[sn]=false;
					continue;
				}
				else if (dot1*dot2>-.001) node_to_node=true;



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



				Vect normal=normals[contId][k];



				double pen=v1v.dot(normal);

				if(pen>0e-14 &&  (!contacting[sn] || !twchk)) {
					contacting[sn]=false;
					landed_stick[sn]=false;
					int p=u_index[sn][0];
					if(p<0) 
						p=u_index[sn][1];
					if(p<0) continue;
					lamN.el[p]=0;
					lamT.el[p]=0;
					
					continue;
				}



				if(pen<-10*edgeLength) {
					contacting[sn]=false;
					continue;
				}

				normalIndex[contId][i]=k;

				double beta=0;
				double alpha=0;


				if(node_to_node){
					if(Math.abs(dot1)<1e-3){
						beta=0;			
						alpha=1-beta;

					}else if(Math.abs(dot2)<1e-3){
						beta=1;			
						alpha=1-beta;
					}
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

				if(mu!=0) {
					if(!stick[sn]&& !contacting[sn]){
						stick[sn]=true;
						landed_stick[sn]=true;
					}
				}


				contacting[sn]=true;


				//util.pr("node "+sn+" contacted.");


				Vect deltaDisp=u.sub(u1.times(alpha).add(u2.times(beta)));
				///	deltaDisp.times(1e9).hshow();
				double proj=deltaDisp.dot(normal);


				Vect disp_tang=deltaDisp.sub(normal.times(proj));

				//	util.pr("-----< "+disp_tang.norm());

				Vect tang=disp_tang.normalized();
				//if(tang.norm()==0) 
				tang=edgeDir.deepCopy();


				tangentials[contId][k]=tang.deepCopy();

				int p=u_index[sn][0];
				if(p<0) 
					p=u_index[sn][1];
				if(p<0) continue;

			//	gap0.el[p]=node.getCoord().sub(node1.getCoord().times(alpha).add(node2.getCoord().times(beta))).dot(normal);

				gap.el[p]=pen;

				//	util.pr(sn+"  gap "+	gap.el[p]);
				//double s0=node.getCoord().sub(node1.getCoord().times(alpha).add(node2.getCoord().times(beta))).dot(tang);
				//	slide_prev.el[p]=slide.el[p];
				//if(!landed_stick[sn])
				slide.el[p]=deltaDisp.dot(tang);
				//	else
				//	slide.el[p]=0;
				//	deltaDisp.hshow();
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
		G_stk_pr=G_stk.deepCopy();

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


	///	gap=Gc.mul(u).add(gap0);


	for(int i=1;i<=model.numberOfNodes;i++){
		//for(int k=0;k<model.dim;k++){
		int p=u_index[i][0];
		if(p<0)  p=u_index[i][1];;
		if(p<0) continue;
		//	if(gap.el[p]!=0) util.pr(i+"    "+gap.el[p]);

		if(gap.el[p]>0 && !rmv[i]){
			//	util.pr(p+"  "+gap.el[p]);

			contacting[i]=false;
			//lamN.el[p]=0;
			//lamT.el[p]=0;
			rmv[i]=true;
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
				if(lamN.el[p]<0){
					double abs_lamT=Math.abs(lamT.el[p]);
					double muFn=mu*Math.abs(lamN.el[p]);
					util.ph("node "+sn+ "  p "+p+ " lamT:  ");
					util.pr(lamT.el[p]," %4.4e");
					util.ph("node "+sn+ " lamN:  ");
					util.pr(lamN.el[p]," %4.4e");
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
							//if(method==1)	lamT.el[p]=0;
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

				util.pr("node "+sn+" stick "+stick[sn]);
				///	util.pr("u "+model.node[sn].getU(p)+" uref "+ref_stick.el[p]);
			}
			else
				util.pr("node "+sn+" free ");
		}
	}

	/// G_stkt.shownzA();

	//	Kcf.shownzA();

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

	Vect pg=gap.times(pf);
	Fc=Gct.mul(pg);
	//Fc=Kc.smul(u);
	//Fcf=Kcf.smul(u);
	Vect ps=slide.times(pft);

	Fcf=G_stkt.mul(ps);

	dF=dF.sub(Fc).sub(Fcf); //
	

	aug_N=Gct.mul(lamN);;
	aug_T=Gcft.mul(lamT);;


	dF=dF.sub(aug_N).sub(aug_T);

	return dF;

}


/*	private Vect calcResidual(SpMat Ks,Vect u, Vect b){

		Vect Fint=model.Ks.smul(u);

		Vect dF=b.sub(Fint);

		//Fc=Kc.smul(u);
		//Fcf=Kcf.smul(u);
	
	//	util.pr("_______________________");
		//Gct.shownzA();

	//	new SpVect(Fc).shownzA();
		Vect ps=slide.times(pft);

		for(int i=0;i<G_stk.nRow;i++)
			if(G_stk.row[i].norm()==0)
				ps.el[i]=0;
	//	new SpVect(slide).shownzA();


		if(G_stk_pr!=null){
		Vect psp=G_stk_pr.mul(ref_stick).times(pft);
		//SpMatAsym M=G_stk_pr.deepCopy();
	//	M.transpose(50);
		//Vect f=M.mul(psp);
	//	Fcf=Fcf.add(f);
	//	ps=ps.sub(psp);
		}
		Fcf=G_stkt.mul(ps);//.sub(Kcf.smul(ref_stick));
	//	G_stkt.shownzA();
	//	new SpVect(Fcf).shownzA();
	//	util.plot(gap);
	//node_node.shownzA();
		//Kc.shownzA();
	//	Fcf=Fcf.add(aug_T);

	//	new SpVect(u).shownzA();

		dF=dF.sub(Fc).sub(Fcf); //

		dF=dF.sub(aug_N).sub(aug_T);

	//	new SpVect(lamT).shownzA();
	//	new SpVect(aug_T).shownzA();
		//new SpVect(dF).shownzA();
	//	util.pr("ooooooooooooooooooooooooooooo");
			return dF;

	}

 */
	public static void main(String[] args){

		new Main();
	}

}


