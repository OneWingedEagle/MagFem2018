package fem;

import static java.lang.Math.abs;

import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class MagMatrix {

	private Calculator calc;

	public  MagMatrix(){}

	public  MagMatrix(Model model){

		this.calc=new Calculator(model);
	}
	



	public void setMagMat(Model model){
		
		setReactMat(model);
		if(model.analysisMode>0)
			setConductMat(model);
		
		if(model.analysisMode>1 && model.dim==3){
			
			setA_Phi_couplingMat( model);
			
			if(!model.AC){
	
				for(int i=0;i<model.numberOfVarNodes;i++)
					model.Hs.row[i+model.numberOfUnknownEdges]=model.Ps.row[i].deepCopy();
			
				model.Qs.times(model.dt);

				for(int i=0;i<model.numberOfVarNodes;i++){
					model.Hs.row[i+model.numberOfUnknownEdges]=model.Hs.row[i+model.numberOfUnknownEdges].augh(	model.Qs.row[i]);
				}
				

			}
		}
		

	}
	
	public void setReactMat(Model model){		


		double eps=1e-10,cPB=model.cpb; 
		boolean nonLinear;
		double[][] H1=new double[model.nElEdge][model.nElEdge];

		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;
	
		int[] nz=new int[model.numberOfUnknowns];

		model.Hs=new SpMat(model.numberOfUnknowns);

	

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}




		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}

		//model.RHS=new Vect(model.numberOfUnknowns);
		model.HkAk= new Vect(model.numberOfUnknowns);



		for(int ir=1;ir<=model.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				nonLinear=model.element[i].isNonlin();

					H1=this.calc.He(model,nBH,nLam,i,nonLinear,false,false,false);
			

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;
							if(nonLinear){
								model.H2[j][k]*=cPB;
							}
						

						}	

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {

							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}
						
	
						
						
						model.HkAk.el[matrixRow]+=H1[j][k]*Ak;

						if(model.edge[columnEdgeNumb].edgeKnown  ){
						
							continue;
						}



						//=======================

						if(nonLinear){

							H1[j][k]+= model.H2[j][k];
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	


							if(abs(H1[j][k])>eps  ){	

								model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];


								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
									model.Hs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);
						
						}


					}			
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


	}
	
	public void setConductMat(Model model){		





		double eps=1e-10,cPB=model.cpb; 
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,ext=6;

		int[] nz=new int[model.numberOfUnknowns];
		
		model.Ss=new SpMat(model.numberOfUnknowns);
			for(int i=0;i<model.numberOfUnknowns;i++){

				model.Ss.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
			}
		




		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}




		for(int ir=1;ir<=model.numberOfRegions;ir++){


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
		

				if(!model.element[i].isConductor()) continue;


				this.calc.He(model,0,0,i,false,true,false,false);




				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];

						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

								model.H3[j][k]*=cPB;

						}	

						//===========================

						if(model.edge[columnEdgeNumb].edgeKnown){

								continue;
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Ss.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	


							if(abs(model.H3[j][k])>eps  ){	

								model.Ss.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Ss.row[matrixRow].el[nz[matrixRow]]=model.H3[j][k];


								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Ss.row[matrixRow].nzLength-1){
									model.Ss.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

						
								model.Ss.row[matrixRow].addToNz(model.H3[j][k],m);
							}
						}


					}			
				}
			}
		


			for(int i=0;i<model.numberOfUnknownEdges;i++){
				model.Ss.row[i].sortAndTrim(nz[i]);
			}


			if(model.eddyTimeIntegMode==0 || model.eddyTimeIntegMode==1){

				model.Ss.times(1.0/model.dt);
				

			}



	}
	

	public void setA_Phi_couplingMat(Model model){
	
		
	model.Ps=getPs(model);


	model.Qs=getQs(model);
	


}
	
	public void setRHS(Model model){
		setRHS(model,true);
	}
	
	public void setRHS(Model model,boolean includeM){		

		model.RHS=new Vect(model.numberOfUnknowns);

		int matrixRow=0, rowEdgeNumb;
		
		double[] Cj=new double[model.nElEdge];
		
		boolean hasM, hasJ;
		
		Vect J=null;


		for(int ir=1;ir<=model.numberOfRegions;ir++){

			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				hasJ=model.element[i].hasJ();

				hasM=model.element[i].hasM();
				
				if(!hasJ && !hasM) continue;

				
				if(hasJ){
					J=model.element[i].getJ();
				}


				this.calc.He(model,0,0,i,false,false,hasJ,hasM);


				if(hasJ ){
					for(int j=0;j<model.nElEdge;j++){

						if(model.dim==2)
							Cj[j]=J.el[2]*model.Cj2d[j];
						else{
							Cj[j]=J.dot(model.Cj[j]);
							
						}


					}

				}


				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					//===========  right-hand side
					if( hasM && includeM  ){					
						model.RHS.el[matrixRow]+=model.C[j];	
					}

					if(hasJ ){

						model.RHS.el[matrixRow]+=Cj[j];	


					}

					
				}
			}
		}


	}


/*


	public void setMagMatEdge(Model model){		

		if(model.analysisMode>0 && model.circuit && model.eddyTimeIntegMode<=-2){

			setMagMatEdgeCircuit(model);
			return;
		}


		boolean fillSs=true;//(model.analysisMode>0 && model.Ss==null);


		double eps=1e-10,cPB=model.cpb; 
		boolean nonLinear,eddy;
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		Vect J=new Vect();
		double[] Cj=new double[model.nElEdge];
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;
		boolean hasJ,hasM;

		int[] nz=new int[model.numberOfUnknowns];
		int[] nzs=new int[model.numberOfUnknowns];
		model.Hs=new SpMat(model.numberOfUnknowns);

		if(fillSs){
			model.Ss=new SpMat(model.numberOfUnknowns);
			for(int i=0;i<model.numberOfUnknowns;i++){

				model.Ss.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
			}
		}

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}




		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}

		model.RHS=new Vect(model.numberOfUnknowns);
		model.HkAk= new Vect(model.numberOfUnknowns);



		for(int ir=1;ir<=model.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				hasJ=model.element[i].hasJ();
				if(hasJ){
					J=model.element[i].getJ();

				}

				hasM=model.element[i].hasM();
				nonLinear=model.element[i].isNonlin();

				if(fillSs &&  model.element[i].isConductor())
					eddy=true;
				else
					eddy=false;



				H1=this.calc.He(model,nBH,nLam,i,nonLinear,eddy,hasJ,hasM);


				if(hasJ ){

					for(int j=0;j<model.nElEdge;j++){

						if(model.dim==2)
							Cj[j]=J.el[2]*model.Cj2d[j];
						else
							Cj[j]=J.dot(model.Cj[j]);


					}

				}


				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					//===========  right-hand side
					if( hasM   ){					
						model.RHS.el[matrixRow]+=model.C[j];	
					}

					if(hasJ ){

						model.RHS.el[matrixRow]+=Cj[j];	


					}

					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;
							if(nonLinear){
								model.H2[j][k]*=cPB;
							}
							if(eddy)
								model.H3[j][k]*=cPB;

						}	

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {

							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}


						if(model.edge[columnEdgeNumb].edgeKnown || model.nonLin ){

							model.HkAk.el[matrixRow]+=H1[j][k]*Ak;


							if(model.edge[columnEdgeNumb].edgeKnown )
								continue;
						}



						//=======================

						if(nonLinear){

							H1[j][k]+= model.H2[j][k];
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	


							if(abs(H1[j][k])>eps  ){	

								model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];


								if(fillSs && model.H3!=null && abs(model.H3[j][k])>eps ){

									model.Ss.row[matrixRow].index[nz[matrixRow]]=columnIndex;

									model.Ss.row[matrixRow].el[nz[matrixRow]]=model.H3[j][k];
									nzs[matrixRow]++;

									if(nzs[matrixRow]==model.Ss.row[matrixRow].nzLength-1){
										model.Ss.row[matrixRow].extend(ext);
									}

								}

								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
									model.Hs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);
							if(fillSs && eddy){
								model.Ss.row[matrixRow].addToNz(model.H3[j][k],m);
							}
						}


					}			
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


		if(fillSs){

			for(int i=0;i<model.numberOfUnknownEdges;i++){
				model.Ss.row[i].sortAndTrim(nzs[i]);

			}


			if(model.eddyTimeIntegMode==0 || model.eddyTimeIntegMode==1){

				model.Ss.times(1.0/model.dt);

			}


		}


		if( model.analysisMode==2 && model.dim==3){
			
			SpMat T;
			if(step==0){
				T=new SpMat(model.numberOfVarNodes, model.numberOfUnknownEdges,0);

			}
			else
			{
				T=getPs(model);

				//T.times(model.dt);


				Vect x=T.smul(model.getUnknownA());
				for( int i=0;i<x.length;i++){
					model.RHS.el[model.numberOfUnknownEdges+i]+=x.el[i];
				}
			//}

			for(int i=0;i<model.numberOfVarNodes;i++)
				model.Hs.row[i+model.numberOfUnknownEdges]=T.row[i].deepCopy();

			T=getQs(model);

			T.times(model.dt);

			for(int i=0;i<model.numberOfVarNodes;i++){
				model.Hs.row[i+model.numberOfUnknownEdges]=model.Hs.row[i+model.numberOfUnknownEdges].augh(T.row[i]);
			}

		}


	}
*/

	public void setMagMatEdgeCircuit(Model model){	



		boolean fillSs=(model.Ss==null);

		//************************
		// for plunger moving mesh
		//fillSs=true;
		//******************

		double eps=1e-20,cPB=model.cpb; 
		boolean nonLinear,eddy;
		Vect J=new Vect();
		double[] Cj=new double[model.nElVert];
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;
		boolean hasJ,hasM,circuit;

		int[] nz=new int[model.numberOfUnknowns];
		int[] nzs=new int[model.numberOfUnknowns];
		int nzStranded=0;
		model.Hs=new SpMat(model.numberOfUnknowns);

		if(fillSs){
			model.Ss=new SpMat(model.numberOfUnknowns);
			for(int i=0;i<model.numberOfUnknowns;i++){

				model.Ss.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
			}
		}

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}


		if(fillSs){

			model.lastRowsAll=new SpVect[model.numberOfCurrents];

			for(int i=0;i<model.numberOfCurrents;i++)
				model.lastRowsAll[i]=new SpVect(model.numberOfUnknowns,model.numberOfUnknownEdges);

			model.lastRows=new SpVect[model.numberOfUnknownCurrents];



		}


		model.RHS=new Vect(model.numberOfUnknowns);

		model.HkAk= new Vect(model.numberOfUnknowns);
		model.HpAp= new Vect(model.numberOfUnknowns);




		for(int ir=1;ir<=model.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			circuit=model.region[ir].circuit;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				hasJ=model.element[i].hasJ();
				if(hasJ)
					J=model.element[i].getJ();

				hasM=model.element[i].hasM();
				nonLinear=model.element[i].isNonlin();

				if(fillSs &&  model.element[i].isConductor())
					eddy=true;
				else
					eddy=false;

				H1=this.calc.He(model,nBH,nLam,i,nonLinear,eddy,hasJ||circuit,hasM);


				if(hasJ ){

					for(int j=0;j<model.nElEdge;j++){

						if(model.dim==2)
							Cj[j]=J.el[2]*model.Cj2d[j];
						else
							Cj[j]=J.dot(model.Cj[j]);
					}
				}

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

					//===========  right-hand side
					if( hasM )	{					
						model.RHS.el[matrixRow]+=model.C[j];	
						//util.pr(model.C[j]);
					}

					if(hasJ )					
						model.RHS.el[matrixRow]+=Cj[j];	



					if( fillSs && model.region[ir].circuit)	{

						double norm=0;
						if(model.dim==2)
							norm=abs(model.Cj2d[j]);
						else
							norm=model.Cj[j].norm();


						if(norm>eps  ){

							int ix=model.region[ir].currentIndex;

							int n=util.search(model.lastRowsAll[ix].index,nzStranded-1,matrixRow);


							double val=0;

							if(model.dim==2){
								val=-model.region[ir].NtS*model.Cj2d[j];
							}

							else if(model.dim==3){
								Vect t=model.node[model.edge[rowEdgeNumb].endNodeNumber[1]].getCoord().sub
										(model.node[model.edge[rowEdgeNumb].endNodeNumber[0]].getCoord());
								t.normalize();

								val=-model.region[ir].NtS*model.Cj[j].dot(t);

							}


							if(n<0){


								model.lastRowsAll[ix].el[nzStranded]+=val;
								model.lastRowsAll[ix].index[nzStranded]=matrixRow;
								nzStranded++;

							}

							else{


								model.lastRowsAll[ix].addToNz(val,n);

							}



						}
					}



					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];

						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;

							if(nonLinear)
								model.H2[j][k]*=cPB;
							if(eddy)
								model.H3[j][k]*=cPB;

						}	

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {
							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}

						if(model.edge[columnEdgeNumb].edgeKnown || model.nonLin ){


							model.HkAk.el[matrixRow]+=H1[j][k]*Ak;


							if(model.edge[columnEdgeNumb].edgeKnown )
								continue;
						}

						if(!model.edge[columnEdgeNumb].edgeKnown && !model.nonLin ){

							model.HpAp.el[matrixRow]+=H1[j][k]*Ak;

						}


						//=======================

						if(nonLinear)
							H1[j][k]+= model.H2[j][k];


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;


						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);

						if(m<0)
						{	

							if(abs(H1[j][k])>eps  ){	

								model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];

								if(!model.region[ir].circuit &&fillSs && abs(model.H3[j][k])>eps ){

									model.Ss.row[matrixRow].index[nz[matrixRow]]=columnIndex;

									model.Ss.row[matrixRow].el[nz[matrixRow]]=model.H3[j][k];
									nzs[matrixRow]++;

									if(nzs[matrixRow]==model.Ss.row[matrixRow].nzLength-1){
										model.Ss.row[matrixRow].extend(ext);
									}						

								}

								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
									model.Hs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);
							if(fillSs && !model.region[ir].circuit )
								model.Ss.row[matrixRow].addToNz(model.H3[j][k],m);
						}


					}			
				}
			}
		}

		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


		if(fillSs){

			for(int i=0;i<model.numberOfUnknownEdges;i++){
				model.Ss.row[i].sortAndTrim(nzs[i]);
			}

			model.Ss.times(1.0/model.dt);

			for(int j=0;j<model.lastRowsAll.length;j++)
				model.lastRowsAll[j]=model.lastRowsAll[j];


			for(int j=0;j<model.lastRows.length;j++){
				int nr1=model.unCurRegNumb[j];

				int kk=model.region[nr1].currentIndex;
				model.lastRows[j]=model.lastRowsAll[kk].deepCopy();
				for(int ir=1;ir<model.numberOfRegions;ir++){
					if(model.region[ir].curMap1!=nr1 || !model.region[ir].circuit) continue;

					int mm=model.region[ir].currentIndex;
					double kf=model.region[ir].currCoef1;
					model.lastRows[j]=model.lastRows[j].addGeneral(model.lastRowsAll[mm].times(kf));
				}

			}
		}

		if( model.analysisMode==2 && model.dim==3){
			SpMat T;
			
				T=getPs(model);
				Vect x=T.smul(model.getUnknownA());
				for( int i=0;i<x.length;i++){
					model.RHS.el[model.numberOfUnknownEdges+i]+=x.el[i];
				}
			
			for(int i=0;i<model.numberOfVarNodes;i++)
				model.Hs.row[i+model.numberOfUnknownEdges]=T.row[i].deepCopy();

			T=getQs(model);

			for(int i=0;i<model.numberOfVarNodes;i++){
				model.Hs.row[i+model.numberOfUnknownEdges]=model.Hs.row[i+model.numberOfUnknownEdges].augh(T.row[i]);
			}

		}


	}


	private  SpMat getQs(Model model){

		double eps=1e-6;
		double[][] He;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[model.numberOfVarNodes];

		SpMat Qs=new SpMat(model.numberOfVarNodes, model.numberOfVarNodes,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;

				He=this.calc.Qe(model,i);
				

				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){



					if(!model.node[vertNumb[j]].isPhiVar() || model.node[vertNumb[j]].isPhiKnown() ) continue;

					matrixRow=model.nodeVarIndex[vertNumb[j]]-1;

					for(int k=0;k<model.nElVert;k++){
						nodeNumber=vertNumb[k];
						columnIndex=model.nodeVarIndex[nodeNumber]-1;								
						if(columnIndex>matrixRow) continue;


						m=util.search(Qs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(He[j][k])>eps ){
								Qs.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								Qs.row[matrixRow].el[nz[matrixRow]++]=He[j][k];
							}

						}

						else{

							Qs.row[matrixRow].addToNz(He[j][k],m);
						}

					}			
				}
			}
		}

		Qs.sortAndTrim(nz);

		return Qs;
	}


	private  SpMat getPs(Model model){

		double ePs=1e-6;
		double[][] He;

		int m,edgeNumber,nodeNumber,columnIndex,matrixRow;
		int[] nz=new int[model.numberOfVarNodes];

		SpMat Ps=new SpMat(model.numberOfVarNodes, model.numberOfUnknownEdges,model.nNodEd);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				
				if(!model.element[i].isConductor()) continue;

				He=this.calc.Pe(model,i);

				int[] vertNumb=model.element[i].getVertNumb();
				int[] edgeNumb=model.element[i].getEdgeNumb();


				for(int j=0;j<model.nElVert;j++){
					nodeNumber=vertNumb[j];
					if(!model.node[nodeNumber].isPhiVar() || model.node[nodeNumber].isPhiKnown()  ) continue;
					matrixRow=model.nodeVarIndex[nodeNumber]-1;

					for(int k=0;k<model.nElEdge;k++){		
						edgeNumber=edgeNumb[k];
						if(model.edge[edgeNumber].edgeKnown) continue;

						columnIndex=model.edgeUnknownIndex[edgeNumber]-1;
						m=util.search(Ps.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	
							if(abs(He[j][k])>ePs ){
								Ps.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								Ps.row[matrixRow].el[nz[matrixRow]++]=He[j][k];
							}

						}

						else{

							Ps.row[matrixRow].addToNz(He[j][k],m);
						}

					}			
				}
			}
		}
		Ps.sortAndTrim(nz);
		
		


		return Ps;
	}

	public void setPODReactMat(Model model,int order){
		


		double eps=1e-10,cPB=model.cpb; 
		boolean nonLinear;
		double[][] H1=new double[model.nElEdge][model.nElEdge];

		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;
	
		int[] nz=new int[model.numberOfUnknowns];

		model.Hs=new SpMat(model.numberOfUnknowns);

	

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}




		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}

		//model.RHS=new Vect(model.numberOfUnknowns);
		//model.HkAk= new Vect(model.numberOfUnknowns);



		for(int ir=1;ir<=model.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				nonLinear=model.element[i].isNonlin();

				H1=this.calc.He(model,nBH,nLam,i,nonLinear,false,false,false);

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;
							if(nonLinear){
								model.H2[j][k]*=cPB;
							}
						

						}	

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {

							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}


						if(model.edge[columnEdgeNumb].edgeKnown  ){
						
							continue;
						}



						//=======================

						if(nonLinear && order==1){

							H1[j][k]+= model.H2[j][k];
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	


							if(abs(H1[j][k])>eps  ){	

								model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];


								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
									model.Hs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);
						
						}


					}			
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


	
	}



}
