package femSolver;

import static java.lang.Math.abs;

import java.util.Arrays;

import fem.Calculator;
import fem.Model;
import fem.PhiCoil;
import math.Mat;
import math.SpMat;
import math.Vect;
import math.util;


public class StaticElectricSolver{
	int stepNumb;

	private SpMat phiMat;
	private Vect RHS;
	private Calculator calc;
	private int[] phiVarIndex;
	private boolean byCircuit;
	private int numberOfUnknownPhis;
	double vps_volatge;
	public boolean open_vps;
	
	public StaticElectricSolver(){

		byCircuit=true;	
		vps_volatge=1;
	}

	public void setMatrix(Model model){

		this.calc=new Calculator(model);

		setPhiMat(model);


	}

	public Vect solve(Model model ){

		Vect x=new Vect(numberOfUnknownPhis);

		if(numberOfUnknownPhis==0) return x;

		SpMat L=new SpMat();


		model.solver.terminate(false);

		//RHS.show();
		SpMat Ks=phiMat.deepCopy();

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

		phiMat=new SpMat(numberOfUnknownPhis, numberOfUnknownPhis,model.nNodNod);

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


						m=util.search(phiMat.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(Ke.el[j][k])>eps ){

								//===========================
								if(nz[matrixRow]==phiMat.row[matrixRow].nzLength-1){
									phiMat.row[matrixRow].extend(ext);
								}
								//===========================

								phiMat.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								phiMat.row[matrixRow].el[nz[matrixRow]++]=Ke.el[j][k];
							}

						}

						else{

							phiMat.row[matrixRow].addToNz(Ke.el[j][k],m);
						}

					}			
				}
			}
		}
		if(byCircuit){
			double Rext=1e-10;
			matrixRow=numberOfUnknownPhis-1;
			nz[matrixRow]=model.phiCoils.length+1;
			for(int ic=0;ic<model.phiCoils.length;ic++){
				phiMat.row[matrixRow].index[ic]=matrixRow-model.phiCoils.length+ic;
				phiMat.row[matrixRow].el[ic]=-1;
			}
				phiMat.row[matrixRow].index[model.phiCoils.length]=matrixRow;
				phiMat.row[matrixRow].el[model.phiCoils.length]=-Rext;
		}


		phiMat.sortAndTrim(nz);

		//	util.pr("PHIMAT");
		
	//	Vect v=phiMat.diag();
		
	////	for(int k=v.length-5;k<v.length;k++)
	//		util.pr(v.el[k]);
		//		phiMat.shownz();
	}


	public void setRHS(Model model){		

		Mat Ke;

		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;
		}

		RHS=new Vect(numberOfUnknownPhis);

		Vect elemPhiVec=new Vect(model.nElVert);

		Vect elemRHS=new Vect(model.nElVert);

		int matrixRow=0, columnIndex;

		boolean isConductive;
		double conductivity;

		if(!byCircuit)
			for(int ir=1;ir<=model.numberOfRegions;ir++){

				if(!model.region[ir].isConductor) continue;

				conductivity=model.region[ir].getSigma().el[0];

				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					isConductive=model.element[i].isConductor();


					if(!isConductive) continue;

					int[] vertNumb=model.element[i].getVertNumb();


					boolean elemHasKnownPhi=false;
					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];
						if(model.node[nodeNumber].isPhiVar() && model.node[nodeNumber].isPhiKnown() ){
							elemHasKnownPhi=true;
							break;
						}
					}

					if(!elemHasKnownPhi) continue;

					Ke=this.calc.elemPhiMat(model,i);

					Ke=Ke.times(conductivity);

					elemPhiVec.timesVoid(0);


					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];
						if(model.node[nodeNumber].isPhiVar() && model.node[nodeNumber].isPhiKnown() ){
							elemPhiVec.el[k]=model.node[nodeNumber].getPhi();
						}
					}

					elemRHS=Ke.mul(elemPhiVec);

					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];

						if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){

							matrixRow=phiVarIndex[nodeNumber]-1;
							RHS.el[matrixRow]+=-elemRHS.el[k];
						}
					}



				}


			}


		if(byCircuit){
		//	for(int ic=0;ic<model.phiCoils.length;ic++){
			//	PhiCoil coil=model.phiCoils[ic];	
				
		//}
			RHS.el[RHS.length-1]=-vps_volatge;
		}
	}


	public void setRHS(Model model,double factor){		



		RHS=new Vect(numberOfUnknownPhis);

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
		if(byCircuit){
			for(int ic=0;ic<model.phiCoils.length;ic++){
				PhiCoil coil=model.phiCoils[ic];	

				phiVarIndex[coil.infaceNodes[0]]=ix++;
				model.node[coil.infaceNodes[0]].setPhiKnown(false);
				for(int i=1;i<coil.infaceNodes.length;i++){
					model.node[coil.infaceNodes[i]].setPhiKnown(false);
					phiVarIndex[coil.infaceNodes[i]]=phiVarIndex[coil.infaceNodes[0]];
				}

			
			}
			ix++; // for current
			
		}

		
		numberOfUnknownPhis=ix-1;
	//	util.pr("------->>>     "+numberOfUnknownPhis);
	}


	
	public void openVPS(Model model){

		if(byCircuit){
				int row=numberOfUnknownPhis-1;
				phiMat.row[row].el[model.phiCoils.length]=-1e12;
		}
	}


}





