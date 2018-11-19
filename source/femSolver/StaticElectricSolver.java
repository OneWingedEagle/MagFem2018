package femSolver;

import static java.lang.Math.abs;

import java.util.Arrays;

import fem.Calculator;
import fem.Model;
import fem.PhiCoil;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticElectricSolver{
	int stepNumb;

	public SpMat phi_matrix;
	public SpMat t_matrix;
	public SpMat conductiveMat;
	public Vect RHS;
	private Calculator calc;
	private int[] phiVarIndex;
	private boolean byCircuit;
	private int numberOfUnknownPhis,numberOfUnknowns,nCurrents;
	double vps_volatge;
	public boolean open_vps;
	
	public StaticElectricSolver(){

		byCircuit=true;	
		vps_volatge=1;
	}

	public void setMatrix(Model model){

		this.calc=new Calculator(model);

		setPhiMat(model);
		setTmat(model);
		AnnexTmat();


	}

	public Vect solve(Model model ){

		Vect x=new Vect(numberOfUnknowns);

		if(numberOfUnknowns==0) return x;

		SpMat L=new SpMat();


		model.solver.terminate(false);

		//RHS.show();
		SpMat Ks=conductiveMat.deepCopy();

		//	Ks.diag().show();

		//Ks.shownz();

		Vect Ci=Ks.scale(RHS);

		if(open_vps)
			L=Ks.ichol(1.);
		else
			L=Ks.ichol();
		
		double errMax=1e-11;
		if(open_vps) errMax=1e-10;


		x=model.solver.ICCG(Ks,L, RHS,errMax,model.iterMax);


		x.timesVoid(Ci);

		for(int i=1;i<=model.numberOfNodes;i++){
			int index=phiVarIndex[i]-1;	

			if(index>=0){
				model.node[i].setPhi(x.el[index]);	
			}
		}

		return x;

	}

	public  void setPhiMat(Model model){

		double eps=1e-8;

		int ext=10;
		Mat Ke;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[numberOfUnknownPhis];

		phi_matrix=new SpMat(numberOfUnknownPhis, numberOfUnknownPhis,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(!model.region[ir].isConductor) continue;

			double conductivity=model.region[ir].getSigma().el[0];

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;

				Ke=this.calc.elemPhiMat(model,i);
				Ke=Ke.times(conductivity);

				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){

					if(!model.node[vertNumb[j]].isPhiVar() || model.node[vertNumb[j]].isPhiKnown() ) continue;

					matrixRow=phiVarIndex[vertNumb[j]]-1;
					for(int k=0;k<model.nElVert;k++){
						nodeNumber=vertNumb[k];
						columnIndex=phiVarIndex[nodeNumber]-1;								

						if(columnIndex==-1 || columnIndex>matrixRow) continue;	


						m=util.search(phi_matrix.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(Ke.el[j][k])>eps ){

								//===========================
								if(nz[matrixRow]==phi_matrix.row[matrixRow].nzLength-1){
									phi_matrix.row[matrixRow].extend(ext);
								}
								//===========================

								phi_matrix.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								phi_matrix.row[matrixRow].el[nz[matrixRow]++]=Ke.el[j][k];
							}

						}

						else{

							phi_matrix.row[matrixRow].addToNz(Ke.el[j][k],m);
						}

					}			
				}
			}
		}
/*		if(byCircuit){
			double Rext=1e-10;
			util.pr(numberOfUnknownPhis+" "+numberOfUnknowns);
			matrixRow=numberOfUnknowns-1;
			nz[matrixRow]=model.phiCoils.length+1;
			for(int ic=0;ic<model.phiCoils.length;ic++){
				phi_matrix.row[matrixRow].index[ic]=matrixRow-model.phiCoils.length+ic;
				phi_matrix.row[matrixRow].el[ic]=-1;
			}
				phi_matrix.row[matrixRow].index[model.phiCoils.length]=matrixRow;
				phi_matrix.row[matrixRow].el[model.phiCoils.length]=-Rext;
		}

*/
		phi_matrix.sortAndTrim(nz);

	}
	
	
	public  void setTmat(Model model){

		if(nCurrents>1) {
			setTmat2(model);
			return;
		}
			t_matrix=new SpMat(nCurrents,numberOfUnknowns);
			double Rext=1e-10;
			int matrixRow=0;
			t_matrix.row[matrixRow]=new SpVect(numberOfUnknowns,model.phiCoils.length+1);

			for(int ic=0;ic<model.phiCoils.length;ic++){
				t_matrix.row[matrixRow].index[ic]=numberOfUnknownPhis-model.phiCoils.length+ic;
				t_matrix.row[matrixRow].el[ic]=-1;
			}
			
			t_matrix.row[matrixRow].index[model.phiCoils.length]=numberOfUnknowns-1;
			t_matrix.row[matrixRow].el[model.phiCoils.length]=-Rext;
	
	}
	public  void setTmat2(Model model){

	
			t_matrix=new SpMat(nCurrents,numberOfUnknowns);
			double Rext=1e-10;
			int matrixRow=0;
			t_matrix.row[matrixRow]=new SpVect(numberOfUnknowns,model.phiCoils.length);

			for(int ic=0;ic<model.phiCoils.length-1;ic++){
				t_matrix.row[matrixRow].index[ic]=numberOfUnknownPhis-model.phiCoils.length+ic;
				t_matrix.row[matrixRow].el[ic]=-1;
			}
			
			
			t_matrix.row[matrixRow].index[model.phiCoils.length-1]=numberOfUnknowns-2;
			t_matrix.row[matrixRow].el[model.phiCoils.length-1]=-Rext;
			
			matrixRow=1;
			t_matrix.row[matrixRow]=new SpVect(numberOfUnknowns,1);
			t_matrix.row[matrixRow].index[0]=numberOfUnknowns-1;
			t_matrix.row[matrixRow].el[0]=-5;
			t_matrix.shownz();
	
	}
	
	private void AnnexTmat(){

		conductiveMat=new SpMat(numberOfUnknowns);
		
		for(int i=0;i<numberOfUnknownPhis;i++){
			conductiveMat.row[i]=new SpVect(numberOfUnknowns,phi_matrix.row[i].nzLength);

			conductiveMat.row[i].index=phi_matrix.row[i].index;
			conductiveMat.row[i].el=phi_matrix.row[i].el;
		}

		for(int i=0;i<t_matrix.nRow;i++){
			int matrixRow=i+numberOfUnknownPhis;
			util.pr(t_matrix.row[i].nzLength);
		conductiveMat.row[matrixRow]=new SpVect(numberOfUnknowns,t_matrix.row[i].nzLength);
		for(int j=0;j<t_matrix.row[i].nzLength;j++){
			conductiveMat.row[matrixRow].index=t_matrix.row[i].index;
			
			conductiveMat.row[matrixRow].el=t_matrix.row[i].el;
		
		}
	}
		//conductiveMat.size();

	}


	public void setRHS(Model model){		


		RHS=new Vect(numberOfUnknowns);

			RHS.el[RHS.length-1]=-vps_volatge;
		
	}


	public void setRHS(Model model,double factor){		



		RHS=new Vect(numberOfUnknowns);

		Vect elemRHS=new Vect(model.nElVert);

		int matrixRow=0;

		boolean isConductive;
		double conductivity;

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(!model.region[ir].isConductor) continue;

			conductivity=model.region[ir].getSigma().el[0];

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				isConductive=model.element[i].isConductor();


				if(!isConductive) continue;

				int[] vertNumb=model.element[i].getVertNumb();


				boolean elemHasPhiVar=false;
				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];

					if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){
						elemHasPhiVar=true;
						break;
					}
				}

				if(!elemHasPhiVar) continue;

				elemRHS=this.calc.elemPhiVect(model,i);

				elemRHS=elemRHS.times(conductivity);


				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];

					if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){
						matrixRow=phiVarIndex[nodeNumber]-1;
						RHS.el[matrixRow]+=-elemRHS.el[k];
					}
				}



			}


		}
		
		if(!open_vps)
			RHS.el[RHS.length-1]=-vps_volatge;
//RHS.show();

	}

	public  void setBoundaryCondition(Model model){




		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;

		}


		for(int ir=1;ir<=model.numberOfRegions;ir++){


			if(!model.region[ir].isConductor) continue;

			if(coilIndices[ir]<0) continue;
	

			PhiCoil coil=model.phiCoils[coilIndices[ir]];
			
			int[] infaceNodes1=new int[model.numberOfNodes];
			boolean[] nc=new boolean[model.numberOfNodes];

			int nx=0;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;


				if(model.element[i].getRegion()!=coil.regNo) continue;

		
				int[] vertNumb=model.element[i].getVertNumb();

				for(int k=0;k<model.nElVert;k++){

					int nodeNumber=vertNumb[k];

					model.node[nodeNumber].setPhiVar(true);

					Vect coord=model.node[nodeNumber].getCoord();

				
					if(coord.el[0]>=coil.faceBox[0][0] &&coord.el[0]<=coil.faceBox[0][1]&&
							coord.el[1]>=coil.faceBox[0][2] &&coord.el[1]<=coil.faceBox[0][3] &&
							(model.dim==2 || coord.el[2]>=coil.faceBox[0][4] &&coord.el[2]<=coil.faceBox[0][5])){
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(1.);
						if(nc[nodeNumber]==false){
							infaceNodes1[nx++]=nodeNumber;
							nc[nodeNumber]=true;
						}

					}

					else if(coord.el[0]>=coil.faceBox[1][0] &&coord.el[0]<=coil.faceBox[1][1] &&
							coord.el[1]>=coil.faceBox[1][2] &&coord.el[1]<=coil.faceBox[1][3] &&
							(model.dim==2 || coord.el[2]>=coil.faceBox[1][4] &&coord.el[2]<=coil.faceBox[1][5])){
						//		else if(model.node[nodeNumber].getCoord(0)<.001){
						if(nc[nodeNumber]==false){
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(1e-10);
					}

					}

				}

			}

			coil.infaceNodes=Arrays.copyOf(infaceNodes1, nx);
			//	util.hshow(coil.infaceNodes);
		}




		setPhiIndices(model);
		
		setCurrentIndices(model);
	}



	public void setPhiIndices(Model model){


		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;
		}

		int ix=1;
		phiVarIndex=new int[model.numberOfNodes+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++){


			if(!model.region[ir].isConductor) continue;

			if(coilIndices[ir]<0) continue;


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;


				int[] vertNumb=model.element[i].getVertNumb();


				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];
					if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown())
					{
						if(phiVarIndex[nodeNumber]==0)
							phiVarIndex[nodeNumber]=ix++;


					}
				}
			}
		}

			for(int ic=0;ic<model.phiCoils.length;ic++){
				PhiCoil coil=model.phiCoils[ic];	

				phiVarIndex[coil.infaceNodes[0]]=ix++;
				model.node[coil.infaceNodes[0]].setPhiKnown(false);
				for(int i=1;i<coil.infaceNodes.length;i++){
					model.node[coil.infaceNodes[i]].setPhiKnown(false);
					phiVarIndex[coil.infaceNodes[i]]=phiVarIndex[coil.infaceNodes[0]];
				}

			
			}
			numberOfUnknownPhis=ix-1;
				

		
	//	util.pr("------->>>     "+numberOfUnknownPhis);
	}

	public void setCurrentIndices(Model model){
		
		nCurrents=1;
		numberOfUnknowns=numberOfUnknownPhis;
		for(int i=0;i<nCurrents;i++)
			numberOfUnknowns++;
		
	}

	
	public void openVPS(Model model){

				//for(int row=0;row<t_matrix.nRow;row++)
					int row=0;
				t_matrix.row[row].el[t_matrix.row[row].nzLength-1]=-1e12;

	}


}





