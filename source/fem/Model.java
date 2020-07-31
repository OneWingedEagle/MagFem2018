package fem;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import femSolver.FEMsolver;
import femSolver.StaticNonlinearMagSolver;
import io.Loader;
import io.Writer;
import materialData.BHCurve;
import materialData.BHSCurve;
import materialData.CurrentWaveForm;
import materialData.LamBCurve;
import materialData.LamBSCurve;
import math.*;
import meshFactory.MeshManipulator;
import static java.lang.Math.*;




public class Model{

	public MagMatAssembler magMat;
	public MechMatrix mechMat;
	public int dim=3,iterMax=10000,nonLinIterMax=30,cpb=1,nRotReg;
	public double cpm=PI/2,Rg=1e5;
	public int numberOfRegions, numberOfNodes,numberOfElements,numberOfEdges,nXYedges,
	nBlocks,nBH,nLam;
	public int nodeNodeConnectionMax=27,nElVert=8,nBoundary=6;
	public int nElEdge=12;
	public double[] spaceBoundary;
	public int[] BCtype;
	public String[] blockMaterial;
	public Vect unifB;
	public Region[] region;
	public PhiCoil[] phiCoils;
	public Node[] node;
	public Element[] element;
	public Edge[] edge;
	public Loader loader=new Loader();
	public Writer writer=new Writer();
	private Combiner combiner=new Combiner();
	public BoundarySet bcSet=new BoundarySet();
	public Calculator femCalc;
	public Force forceCalc;
	public int[] edgeUnknownIndex,unknownEdgeNumber,nodeUnknownIndex,unknownNodeNumber,
	U_unknownIndex,unknownUnumber,A_unknownIndex,unknownAnumber,T_unknownIndex,unknownTnumber;
	public int[] knownEdgeNumber,nodeVarIndex,varNodeNumber,unCurRegNumb;
	public double[] knownEdgeValue;
	public double scaleFactor=1,Bmax1=0,Bmin1=0,Bmax,Bmin,stressMax=0,nodalStressMax=0,
			stressMin,nodalStressMin=0,Jmin,Jmax,Jemin,Jemax,maxDim,minEdgeLength,maxEdgeLength,
			FmsMin,FmsMax=0,FreluctMax=0,FreluctMin=0,uMax=0,AuMax,defScale;
	public int numberOfUnknownEdges,numberOfMappedEdges,numberOfKnownEdges,numberOfVarNodes
	,numberOfKnownPhis,numberOfUnknowns,analysisMode,stressViewCode=-1
	,numberOfUnknownU,numberOfUnknownUcomp,numberOfCurrents,numberOfUnknownCurrents,nNeutral,defMode=1,numberOfUnknownT;
	public boolean deform,thermal,hasJ,hasM,forceLoaded,fluxLoaded,potentialLoaded,
	stressLoaded,forceCalcLoaded,coupled,cpvms,calCurve,nonLin,hasPBC,hasMechPBC=true,rotIndex,centrigForce;
	public int nEdEd=34,nNodNod=27,nEdNod=18,nNodEd=54,nEdgeHasJ;
	static int[][] facetVert={{6,2,5},{7,4,3},{6,7,2},{0,4,1},{4,7,5},{0,1,3}};

	public byte elCode=4;
	public double nu0=1e7/(4*PI);
	public double rotAng,rotSubAng,rotSpeed,meshAngStep,freq=1,dt,errCGmax=1,
			errNRmax=1e-6,errCG_NR=1e-1,errFluxMax=1e-3;
	public int nTsteps,nBegin,nEnd,nInc,nRotorElements,tagx;
	public int coordCode=0,timeIntegMode=0,eddyTimeIntegMode;
	public BHSCurve[] BHS;
	public double[] BHstress;
	public LamBSCurve[] lamBS;
	public LamBCurve[] lamB;
	public BHCurve[] BH;
	public CurrentWaveForm ia,ib,ic,va,vb,vc;
	public String elType="hexahedron";
	public SpMat Hs,Ms,Ks,Ls,Cs,Ss,Ps,Qs,Fs;
	public SpMat Rs,Bs,BtBs;
	public Mat eigVects,bigMat,Q;
	public Vect lams,RHS,bU,bT,HpAp,HkAk;
	public boolean AC,motor,modal,hasTwoNodeNumb,rotateConnect,fullMotor,writeFiles,Tmethod,
	circuit,stranded,wavePC,loadFlux,loadPotentioal,loadPrevMag,loadPrevMech,loadForce,loadDisp,saveFlux,saveForce,saveDisp,saveStress,
	transfer2DTo3D,magAnalysis,mechAnalysis,rotateRotor,axiSym;
	public int forceCalcMode=1,dataType,snapShot=0,POD=-1;
	public double alpha1,alpha2,r1,r2,h1,h2,rm,TrqZ,height=1,gw=1e4,vNeutral,tet,tetp,tetpp,spwmLevel,rayAlpha, rayBeta;
	public int[] PBCpair,threePhaseRegs=new int[3];
	public int[][][] commonNodes;
	public int[] edgeOnFSIndices=null;
	public SpMatSolver solver=new SpMatSolver();
	public StaticNonlinearMagSolver nonlinearSolver=null;
	FEMsolver femSolver=new FEMsolver();;

	public Vect up,upp,xp,Ci,ud,udd,hp;
	public String eddyFolder,filePath,meshFilePath,dataFilePath,eddyFilePath,fluxFilePath,resultFolder,fluxFolderIn;
	public double[][] H2,H3;
	public double[] C,Cj2d;
	public Vect[] Cj;

	public String[][] mechBC;
	public String[][] seepBC=new String[50][20];
	public String[] forceFile;
	public String forceFolder;
	public SpVect[] lastRows,lastRowsAll;
	public Model m2d;
	public Vect[][] forceLamin;
	public int[] mapnr;
	public boolean hasBunif,open_vps;
	public int nCLNstages;
	public Network network;
	
	public ContactAnalysis contact;
	
	public int struc2D=0;// 0: plane stress  1: plain strain, 2:axisymmetric
	public double rpm=0;


	public Model(){}


	public Model(int nRegions, int nElements,int nNodes, String elType){

		this.numberOfRegions=nRegions;
		this.numberOfElements=nElements;

		this.numberOfNodes=nNodes;

		this.setElType(elType);

		region=new Region[this.numberOfRegions+1];
		for(int i=1;i<=this.numberOfRegions;i++)
			region[i]=new Region(dim);

		element=new Element[this.numberOfElements+1];
		for(int i=1;i<=this.numberOfElements;i++)
			element[i]=new Element(elType);

		node=new Node[this.numberOfNodes+1];
		for(int i=1;i<=this.numberOfNodes;i++)
			node[i]=new Node(i,dim);
	}


	public Model(String bun){
		new Model();
		loadMesh(bun);
	}


	public void setFemCalc(){
		femCalc=new Calculator(this);

	}
	public void setForceCalc(){
		forceCalc=new Force(this);

	}

	public void alloc(int nRegions, int nElements,int nNodes, String elType){

		this.numberOfRegions=nRegions;
		this.numberOfElements=nElements;

		this.numberOfNodes=nNodes;

		this.setElType(elType);


		region=new Region[this.numberOfRegions+1];
		for(int i=1;i<=this.numberOfRegions;i++)
			region[i]=new Region(dim);

		element=new Element[this.numberOfElements+1];
		for(int i=1;i<=this.numberOfElements;i++)
			element[i]=new Element(elType);

		node=new Node[this.numberOfNodes+1];
		for(int i=1;i<=this.numberOfNodes;i++)
			node[i]=new Node(i,dim);
	}
	public void loadMesh(String bunFilePath){

		loader.loadMesh(this, bunFilePath);

	}

	public void loadData(String dataFilePath)
	{
		//if(dim==3)
			loader.loadData(this,dataFilePath);
	//	else
		//	loader.loadData2D(this,dataFilePath);
		this.femCalc=new Calculator(this);
		this.forceCalc=new Force(this);	
		setMagMech();

	}

	public Model deepCopy(){

		Model copy=new Model(this.numberOfRegions,this.numberOfElements,this.numberOfNodes,this.elType);

		for(int i=1;i<=this.numberOfRegions;i++){
			copy.region[i].setFirstEl(this.region[i].getFirstEl());
			copy.region[i].setLastEl(this.region[i].getLastEl());
		}


		for(int i=1;i<=this.numberOfElements;i++){
			copy.element[i].setVertNumb(this.element[i].getVertNumb());
			copy.element[i].setEdgeNumb(this.element[i].getEdgeNumb());
			copy.element[i].setRegion(this.element[i].getRegion());
		}


		for(int i=1;i<=this.numberOfNodes;i++){
			copy.node[i].setCoord(this.node[i].getCoord());
		}

		if(this.edge!=null){
			copy.edge=new Edge[this.numberOfEdges+1];
			for(int i=1;i<=this.numberOfEdges;i++){
				//copy.edge[i]=new Edge(this.edge[i].endNodeNumber[0],this.edge[i].endNodeNumber[1]);
				copy.edge[i]=new Edge(this.edge[i].node[0],this.edge[i].node[1]);
				copy.edge[i].A=this.edge[i].A;
				copy.edge[i].Ap=this.edge[i].Ap;
			}
		}

		copy.elCode=this.elCode;
		copy.scaleFactor=this.scaleFactor;

		return copy;
	}

	public Model fill(int nRegs,int nEls,int nNodes, String elType){

		Model copy=new Model(nRegs,nEls,nNodes,elType);
		for(int i=1;i<=this.numberOfRegions;i++){
			copy.region[i].setFirstEl(this.region[i].getFirstEl());
			copy.region[i].setLastEl(this.region[i].getLastEl());
			copy.region[i].setName(this.region[i].getName());
		}

		for(int i=this.numberOfRegions+1;i<=nRegs;i++){
			copy.region[i].setFirstEl(1);
			copy.region[i].setLastEl(0);
		}

		for(int i=1;i<=this.numberOfElements;i++){
			copy.element[i].setVertNumb(this.element[i].getVertNumb());
			copy.element[i].setRegion(this.element[i].getRegion());
		}


		for(int i=1;i<=this.numberOfNodes;i++){
			copy.node[i].setCoord(this.node[i].getCoord());
		}


		copy.scaleFactor=this.scaleFactor;
		return copy;
	}


	public double edgeLength(int i){
	//	double length=node[edge[i].endNodeNumber[1]].getCoord().sub(node[edge[i].endNodeNumber[0]].getCoord()).norm();
		double length=edge[i].node[0].getCoord().sub(edge[i].node[1].getCoord()).norm();

		return length;
	}

	public void setMinEdge(){
		double minEdge=1e40;
		for(int i=1;i<=this.numberOfEdges;i++){
			if(edge[i].edgeLength<minEdge) minEdge=edge[i].edgeLength;
		}
		this.minEdgeLength=minEdge;
	}
	public void setMaxEdge(){
		double maxEdge=0;
		for(int i=1;i<=this.numberOfEdges;i++)
			if(edge[i].edgeLength>maxEdge) maxEdge=edge[i].edgeLength;

		this.maxEdgeLength=maxEdge;
	}

	public  void setMaxDim(){

		maxDim=getRmax();
	}

	public void setHasJ(){

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].hasJ){
			hasJ=true;
			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

					element[i].setHasJ(true);
			}
			}
			else for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setHasJ(false);
		}
			}
	}
	public void setHasM(){
		for(int i=1;i<=numberOfRegions;i++)
			if(region[i].hasM){
				hasM=true;
				break;
			}
	}

	public void setNonLin(boolean nonlin){

		nonLin=nonlin;
		
		if(nonlin) nonlinearSolver=new StaticNonlinearMagSolver();

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].isNonLinear)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
					element[i].setNonlin(true);
		}	

	}
	public void setNonLinToElememts(){

		if(nonLin)
			for(int ir=1;ir<=numberOfRegions;ir++){
				if(region[ir].BHnumber>0){
					region[ir].isNonLinear=true;
					for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
						element[i].setNonlin(true);
				}
			}	

	}

	public void setHasMS(){
		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].isNonLinear){
				region[ir].MS=true;
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
					this.element[i].setHasMS(true);
			}
			else
				region[ir].MS=false;		
		}		
	}

	public void setHasThermal(){		

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].thermal)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

					element[i].setHasThermal(true);

				}

		}	
	}



	public void setMagMech(){
		magMat=new MagMatAssembler(this);
		mechMat=new MechMatrix(this);
	}

	public void femSolver(){
		this.femSolver=new FEMsolver(this);

	}


	public void setMagMat(){
		magMat.setMagMat(this);

	}


	public void setStiffMat(boolean massNeeded){
		mechMat.setStiffMat(this,massNeeded);
	}

	public void setStiffMat(){
		mechMat.setStiffMat(this,false);
	}



	public void setNodalMass(){
		mechMat.setNodalMass(this);

	}

	public Vect getDeformation(){


		return mechMat.getDeformation(this, solver,defMode);
	}


	public Vect getVibration(int iter){
		if(this.timeIntegMode==1)
			return mechMat.getVibration(this, solver,defMode, iter);
		else if(this.timeIntegMode==2) 
			return mechMat.getVibrationCentral(this, solver,defMode,iter);
		else if(this.timeIntegMode==3) 
			return mechMat.getVibrationNewmark(this, solver,defMode,iter);
		else if(this.timeIntegMode==4) 
			return mechMat.getVibrationSS22(this, solver,defMode,iter);
		else if(this.timeIntegMode==5) 
			return mechMat.XEuler(this, solver,defMode,iter);
		else if(this.timeIntegMode==6) 
			return mechMat.IEuler(this, solver,defMode,iter);
		else if(this.timeIntegMode==7) 
			return mechMat.getVibrationVerlet(this, solver,defMode, iter);
		else if(this.timeIntegMode==8) 
			return mechMat.Runge_Kutta(this, solver,defMode,iter);
		else if(this.timeIntegMode==9) 
			return mechMat.IRunge_Kutta(this, solver,defMode,iter);
		else if(this.timeIntegMode==10) 
			return mechMat.leapFrog(this, solver,defMode,iter);
		else

			return null;
	}
	public Vect getVibrationContact(int iter){
		
		return contact.getVibration(this, solver,defMode,iter);
	}

	public Vect getDeformationContact(int iter){
		
		return contact.getDeformation(this, solver,defMode,iter);
	}
	
	public void writeMesh(String bunFilePath){
		writer.writeMesh(this,bunFilePath);
	}

	public void writeMesh(String bunFilePath, boolean deform){
		writer.writeMesh(this,bunFilePath,deform);
	}

	public void writeMeshOneNode(String bunFilePath){
		writer.writeMeshOneNode(this,bunFilePath);
	}

	public void writeModeShape(String bunFilePath , double a){
		writer.writeModeShape(this,bunFilePath,a);

	}

	public void writeMeshTriToQuad2(String bunFilePath , boolean deformed){
		writer.writeMeshTriToQuad2(this,bunFilePath ,deformed);
	}

	public void writeMeshq2h(String bunFilePath , boolean deformed){
		writer.writeMeshq2h(this,bunFilePath ,deformed);
	}

	public void hexaToPyramid(String bunFilePath){
		writer.writeMeshHexaToPyramid(this,bunFilePath);
	}

	public void pyramidToTetra(String bunFilePath){
		writer.writeMeshPyramidToTetra(this,bunFilePath);
	}

	public void writeMeshTriToTri(String bunFilePath,double r1,double r2){
		//	writer.writeMesh323(this,bunFilePath,r1,r2);
		writer.writeMeshTriToTriW(this,bunFilePath);

	}

	public void writeMeshTriToQuad(String bunFilePath){
		writer.writeMeshTriToQuad(this,bunFilePath);
	}

	public void writeMeshQuadToTri(String bunFilePath){
		writer.writeMeshQuadToTri(this,bunFilePath);

	}

	public void writeMeshHexaToPrism(String bunFilePath){
		writer.writeMeshHexaToPrism(this,bunFilePath );
	}
	public void writeMeshHexaToPrism(String bunFilePath,int dir){
		writer.writeMeshHexaToPrism(this,bunFilePath,dir );
	}

	public void prismToHexa(String bunFilePath){
		writer.writeMeshPrismToHexa(this,bunFilePath );
	}

	public void writeData(String dataFilePathOut) {
		try{
			this.writer.copyFile(this.dataFilePath,dataFilePathOut);
			System.out.println(" Data file was written to "+dataFilePathOut);

		}
		catch(IOException e){System.out.println("writing data file failed.");}

	}

	public boolean loadFlux(String fluxFilePath){

		return loader.loadFlux(this,fluxFilePath);

	}	

	public boolean loadFlux(String fluxFilePath, double angDeg){

		return loader.loadFlux(this,fluxFilePath,angDeg);

	}	

	public boolean loadPotential(String file){

		return loader.loadPotential(this,file);

	}
	public boolean loadStress(String stressFilePath){
		return loader.loadStress(this,stressFilePath);
	}



	public boolean loadNodalField(String filePath,int mode){

		return loader.loadNodalField(this,filePath,mode);


	}

	public boolean loadNodalField(String filePath,int mode,double mechAng){

		return loader.loadNodalField(this,filePath,mode,mechAng);


	}


	public boolean loadEdgeField(String edgeFilePath,int mode){
		return loader.loadNodalField(this,edgeFilePath,mode);
	}

	public void reportData(){
		writer.reportData(this);
	}

	public Node[] elementNodes(int i){
		Node[] elementNode=new Node[nElVert];
		int[] vertNumb=element[i].getVertNumb();
		for(int j=0;j<nElVert;j++){

			elementNode[j]=node[vertNumb[j]];
		}

		return elementNode;
	}
	public Edge[] elementEdges(int i){
		Edge[] elementEdge=new Edge[nElEdge];
		int[] edgeNumb=element[i].getEdgeNumb();
		for(int j=0;j<nElEdge;j++)
			elementEdge[j]=edge[edgeNumb[j]];
		return elementEdge;
	}

	public int[] getContainingElement(Vect P){
		if(elCode==0) return getContaining3angElement(P);
		else if(elCode==1) return getContainingQuadElement(P);
		int[] elResult;
		int[] result=new int[7];

		for(int ir=1;ir<=this.numberOfRegions;ir++){

			//	if(ir!=1 && ir!=2) continue;
			//if(ir!=1 ) continue;


			for(int i=this.region[ir].getFirstEl();i<=this.region[ir].getLastEl();i++){

				elResult=elContains(i,P);
				result[0]=0;
				for(int j=0;j<6;j++)
					if(elResult[j]<0){
						result[0]=-1;
						break;
					}

				if(result[0]==-1) continue;		


				result[0]=i;
				for(int j=0;j<6;j++)
					result[j+1]=elResult[j];
				break;

			}
		}



		return result;
	}

	public int[] getContaining3angElement(Vect P){
		int[] ne=new int[1];


		double[] S;
		double S0,St;
		for(int i=1;i<=numberOfElements;i++){
			S0=el3angArea(i);
			S=subS3ang(i,P);
			St=S[0]+S[1]+S[2];

			if(abs(1-St/S0)<1e-6){ ne[0]= i; break;}

		}

		return ne;
	}



	public int[] elContains(int i,Vect P){
		int n;
		int[] result=new int [nBoundary];
		for(int ns=0;ns<nBoundary;ns++){
			n=pointOnFace(i,ns,P);
			result[ns]=n;

		}

		return result;
	}

	public int pointOnFace(int i, int ns,Vect P){
		Vect v0,v1,v2,vP;

		v0=node[element[i].getVertNumb(facetVert[ns][0])].getCoord();
		v1=node[element[i].getVertNumb(facetVert[ns][1])].getCoord();
		v2=node[element[i].getVertNumb(facetVert[ns][2])].getCoord();
		vP=P.sub(v0);

		if(vP.norm()==0)
			return 0;
		double crossdot=v1.sub(v0).cross(v2.sub(v0)).dot(P.sub(v0));
		if(crossdot==0)
			return 0;
		else if(crossdot>1e-10)
			return -1;
		else
			return 1;

	}

	public int[] getContainingQuadElement(Vect P){

		int[] ne=new int[1];

		Vect[] v=new Vect[4];

		for(int i=1;i<=numberOfElements;i++){
			int[] vertNumb=element[i].getVertNumb();
			for(int j=0;j<4;j++)
				v[j]=node[vertNumb[j]].getCoord();

			double S0=v[1].sub(v[0]).cross(v[3].sub(v[0])).norm()+
					v[1].sub(v[2]).cross(v[3].sub(v[2])).norm();
			double S=0;
			for(int j=0;j<4;j++)
			{
				S=S+v[j].sub(P).cross(v[(j+1)%4].sub(P)).norm();
			}
			if(abs(1-S/S0)<1e-6) {ne[0]=i; break;}

		}

		return ne;

	}

	public double getRegionArea(int ir)

	{

		if(this.dim>2) throw new NullPointerException(" Region is not 2D. ");

		double S=0;


		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
			S+=getElementArea( i);

		return S;

	}

	public double getRegionAreaXY(int ir)

	{

		if(this.dim!=3) throw new NullPointerException(" Region is not 3D. ");

		double S=0;


		S=3.08E-4;


		/*for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
				S+=getElementArea( i);
		 */
		return S;

	}

	public double getRegionVolume(int ir)

	{

		if(this.dim!=3) throw new NullPointerException(" Region is not 3D. ");

		double V=0;


		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
			V+=getElementVolume( i);

		return V;

	}

	public double getElementArea(int i){
		if(elCode==0) return el3angArea(i);
		else 	if(elCode==1) return elQuadArea(i);
		else throw new NullPointerException(" Element is not 2D. ");

	}

	public double el3angArea(int i){
		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;
		return S;
	}

	public double el3angMaxCosAng(int i){
		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double dd1= abs(v1.normalized().dot(v2.normalized()));
		v1=v1.times(-1);
		v2=vertexNode[2].getCoord().sub(vertexNode[1].getCoord());
		double dd2= abs(v1.normalized().dot(v2.normalized()));
		v2=v2.times(-1);
		v1=vertexNode[0].getCoord().sub(vertexNode[2].getCoord());
		double dd3= abs(v1.normalized().dot(v2.normalized()));
		double dm=max(dd3,max(dd1,dd2));

		return dm;
	}

	public double elQuadArea(int i){
		Node[] vertexNode=elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[3].getCoord().sub(vertexNode[0].getCoord());
		Vect v3=vertexNode[1].getCoord().sub(vertexNode[2].getCoord());
		Vect v4=vertexNode[3].getCoord().sub(vertexNode[2].getCoord());
		double S=(v1.cross(v2).norm()+v4.cross(v3).norm())/2;
		return S;
	}

	public void setSolvedAL(Vect x){

		for(int i=1;i<=numberOfUnknownEdges;i++){


			edge[unknownEdgeNumber[i]].setSolvedAL(x.el[i-1]);	

		}


		if(this.hasPBC || this.hasTwoNodeNumb){

			for(int i=1;i<=numberOfEdges;i++){
				if(edge[i].map>0 && !edge[i].edgeKnown)
				{
					if(edge[i].aPBC){

						edge[i].setSolvedAL(this.cpb*edge[edge[i].map].A);	

					}


					else{

						edge[i].setSolvedAL(edge[edge[i].map].A);	


					}
				}

			}	
		}




	}


	public void setU(Vect u){

		int ix=0;
		for(int i=1;i<=numberOfUnknownU;i++){
			int nodeNumb=unknownUnumber[i];

			if(node[nodeNumb].isDeformable()){
				for(int k=0;k<dim;k++){
					if(!node[ nodeNumb].is_U_known(k)){
						node[ nodeNumb].setU(k, u.el[ix]);
						ix++;
					}
				}

			}
		}

		if(this.hasPBC || this.hasTwoNodeNumb){

			Mat R=new Mat();
			if(this.dim==2)
				R=util.rotMat2D(this.cpm);		
			else R=util.rotEuler(new Vect(0,0,1), this.cpm);


			for(int i=1;i<=numberOfNodes;i++){
				if( node[i].isDeformable() && node[i].getMap()>0 && !node[i].is_U_known())
				{
					if(this.coordCode==0 && node[i].hasPBC()){

						Vect v=R.mul(node[node[i].getMap()].getU());

						for(int k=0;k<dim;k++){
							node[ i].setU(k, v.el[k]);
						}
					}
					else{
						node[i].setU(node[node[i].getMap()].getU());	

					}
				}

			}	
		}


	}



	public Vect getU(int comp){
		Vect u=new Vect(this.numberOfNodes);
		for(int i=1;i<this.numberOfNodes;i++){
			u.el[i-1]=node[i].u.el[comp];

		}

		return u;
	}


	public void deformMesh(){

		for(int i=1;i<=this.numberOfNodes;i++){
			if(node[i].isDeformable() && node[i].u.norm()!=0)
				node[i].setCoord(node[i].getCoord().add(node[i].getU()));

		}

	}


	private void setElementB(int i){

		if(elCode==0&& !this.axiSym)set3angElementB(i);
		else if(elCode==0) set3angElementAxiB(i);
		else if(elCode==1&& !this.axiSym)setQuadElementB(i);
		else if(elCode==1) setQuadElementAxiB(i);
		else if(elCode==3) setPrismElementB(i);
		else{
			boolean[] edgeDir=element[i].getEdgeReverse();

			Node[] vertexNode=elementNodes(i);
			Vect zero=new Vect(3);
			Mat jac=femCalc.jacobian(vertexNode,zero);
			Vect B;
			Vect[] rotNe=femCalc.rotNe(jac,zero,edgeDir);
			B=getElementB(i,rotNe);
			element[i].setB(B);


		}
	}	



	private void setQuadElementB(int i){

		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] rotNe;

		jac=femCalc.jacobian(vertexNode,zero);
		rotNe=femCalc.rotNeQuad(jac,zero);

		Vect B=getElementB(i,rotNe);

		element[i].setB(B);

	}


	private void setQuadElementAxiB(int i){

		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] rotNe;

		jac=femCalc.jacobian(vertexNode,zero);
		rotNe=femCalc.rotNeQuad(jac,zero);


		double rr=0;
		for(int j=0;j<4;j++){
			double r1=vertexNode[j].getCoord(0)+1e-5;
			rr+=.25/r1;
		}

		//	rr=1;


		Vect B=getElementB(i,rotNe).times(rr);

		//	if(this.getElementCenter(i).el[1]<.052) B=new Vect(0,1);

		element[i].setB(B);

	}

	private void set3angElementAxiB(int i){

		Edge[] elEdge=elementEdges(i);
		double[] A=new double[3];
		for(int j=0;j<3;j++)
			A[j]=elEdge[j].A;


		Node[] vertexNode=elementNodes(i);


		double rr=0;
		for(int j=0;j<3;j++){
			double r1=vertexNode[j].getCoord(0)+1e-6;
			rr+=.33333333/r1;
		}

		Vect[] rotNe=femCalc.rotNe3ang(this.elementNodes(i));
		Vect B=new Vect(2);		
		for(int j=0;j<3;j++)
			B=B.add(rotNe[j].times(A[j]*rr));


		element[i].setB(B);

	}

	private void setPrismElementB(int i){
		
		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(3);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		boolean[] edgeDir=element[i].getEdgeReverse();

		Vect[] rotNe=femCalc.rotNePrism(jac,zero,edgeDir);
		Vect B=getElementB(i,rotNe);

		//B=this.getElementCenter(i).times(new Vect(1,1,0)).normalized().times(-1);

		element[i].setB(B);

	}




	private void set3angElementB(int i){

		Edge[] elEdge=elementEdges(i);
		double[] A=new double[3];
		for(int j=0;j<3;j++)
			A[j]=elEdge[j].A;

		Vect[] rotNe=femCalc.rotNe3ang(this.elementNodes(i));
		Vect B=new Vect(2);		
		for(int j=0;j<3;j++)
			B=B.add(rotNe[j].times(A[j]));

		//	B=new Vect(this.getElement3angA(i),0);

		element[i].setB(B);
	}


	public void setElType(String type){
		elType=type;
		if(type.equals("triangle") ){
			elCode=0;
			nElVert=3;
			nElEdge=3;
			this.dim=2;
		}
		else if(type.equals("quadrangle") ){
			elCode=1;
			nElVert=4;
			nElEdge=4;
			this.dim=2;
		}
		else if(type.equals("tetrahedron") ){
			elCode=2;
			nElVert=4;
			nElEdge=6;
			dim=3;
		}
		else if(type.equals("prism") ){
			elCode=3;
			nElVert=6;
			nElEdge=9;
			dim=3;
		}
		else if(type.equals("hexahedron") ){
			elCode=4;
			nElVert=8;
			nElEdge=12;
			dim=3;
		}

		else if(type.equals("pyramid") ){
			elCode=5;
			nElVert=5;
			nElEdge=8;
			dim=3;
		}
		nBoundary=2*dim;
	}


	private Vect getVeclocity(int i, Vect[] gradN){

		int[] vert=element[i].getVertNumb();
		Vect Kt=element[i].getSigma();
		Vect B=new Vect(dim);
		for(int j=0;j<nElVert;j++)		{
			int nn=vert[j];
			B=B.add(gradN[j].times(node[nn].T));
		}

		return Kt.times(B);

	}


	private Vect getElementB(int i, Vect[] rotNe){

		Edge[] edge=elementEdges(i);
		Vect B=new Vect(dim);
		for(int j=0;j<nElEdge;j++)		{
		B=B.add(rotNe[j].times(edge[j].A));
		}
	

		return B;

	}

	public double getElementQuadA(int i){
		Edge[] elEdge=elementEdges(i);
		double[] Ae=new double[4];
		for(int j=0;j<4;j++)
			Ae[j]=elEdge[j].A;
		double A=0;
		Vect zero=new Vect(2);

		double[] Ne=femCalc.NeQuad(zero);

		for(int j=0;j<nElEdge;j++)	{		
			A= A+Ne[j]*Ae[j];
		}
		return  A;	

	}



	public void setReluctForce(){
		forceCalc.setReluctForce(this);
	}

	public void resetReluctForce(){
		for(int i=1;i<=this.numberOfNodes;i++)
			if(this.node[i].F!=null) this.node[i].F=new Vect(dim);
	}

	public void setMSForce(){
		forceCalc.setMSForce(this);
	}

	public void setThermalForce(){
		forceCalc.setThermalForce(this);
	}

	public void setStressForce(){
		forceCalc.setStressForce(this);
	}

	public double getElementVolume(int i){

		if( this.elCode==2) return getElementVolumeTetra( i);
		if( this.elCode==3) return getElementVolumePrism( i);

		if(dim==2) return 0;
		double vol=0;
		Node[] vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(dim);
		ws=8;
		jac=femCalc.jacobian(vertexNode,localCo);
		//detJac=abs(jac.determinant());
		detJac=jac.determinant();

		vol=detJac*ws;

		return vol;

	}

	public double getElementVolumePrism(int i){

		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=abs(v1.cross(v2).norm())/2;


		double h=vertexNode[3].getCoord().sub(vertexNode[0].getCoord()).norm();
		double vol=S*h;



		return vol;

	}


	public double getElementVolumeTetra(int i){

		double vol=0;
		Node[] vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(dim);
		jac=femCalc.jacobian(vertexNode,localCo);
		detJac=abs(jac.determinant())/6;

		vol=detJac;

		return vol;

	}



	public void magBC(){
		this.bcSet.magBC(this);
	}

	public void setMagBC(){
		this.bcSet.setMagBC(this);
	}

	public void setMechBC(){
		this.bcSet.mechBC(this);
	}

	public void setSeepBC(){
		this.bcSet.seepBC(this);
	}

	public void setEdge(){

		EdgeSet edgeSet=new EdgeSet();
		edgeSet.setEdge(this);

	}


	public void setSliceBounds(){
		this.bcSet.setSliceBounds(this);
	}

	public void setBounds(){
		if(this.coordCode==1)
			this.bcSet.setSliceBounds(this);
		else if(coordCode==0){

			double[] spb=new double[2*this.dim];


			for(int i=0;i<this.dim;i++){
				spb[2*i]=1e10;
				spb[2*i+1]=-1e10;
			}


			for(int i=1;i<=this.numberOfNodes;i++){

				if(this.node[i].getCoord(0)<spb[0]) spb[0]=this.node[i].getCoord(0);
				else if(this.node[i].getCoord(0)>spb[1]) spb[1]=this.node[i].getCoord(0);

				if(this.node[i].getCoord(1)<spb[2]) spb[2]=this.node[i].getCoord(1);
				else if(this.node[i].getCoord(1)>spb[3]) spb[3]=this.node[i].getCoord(1);
				if(this.dim==3){
					if(this.node[i].getCoord(2)<spb[4]) spb[4]=this.node[i].getCoord(2);
					else if(this.node[i].getCoord(2)>spb[5]) spb[5]=this.node[i].getCoord(2);
				}


			}

			this.spaceBoundary=spb;

		}
	}



	public void mapCommonNodes(){
		this.bcSet.mapCommonNodes(this);
	}

	public void setNodeOnBound(){
		this.bcSet.setNodeOnBound(this);
	}

	public void mapPBC(){
		this.bcSet.mapPBC(this);
	}


	public void setInUseNodes(){
		//=====================	identifying the used nodes
		for(int i=1;i<=this.numberOfElements;i++){

			int[] vert=this.element[i].getVertNumb();
			for(int j=0;j<this.nElVert;j++)	{	
				if(vert[j]>0)
					this.node[vert[j]].inUse=true;
			}
		}
	}


	public double getElementMass(int i){

		double elementMass=0;
		double elementRo=element[i].getRo();
		if(elCode==0) return el3angArea(i)*elementRo;
		else if(elCode==1) return elQuadArea(i)*elementRo;
		Node[] vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(dim);
		ws=8;
		jac=femCalc.jacobian(vertexNode,localCo);
		detJac=abs(jac.determinant());

		elementMass=elementRo*detJac*ws;

		return elementMass;

	}

	private void setElementJ0(int i){

		if(this.elCode!=4) return;

		boolean[] edgeDir=element[i].getEdgeReverse();

		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(3);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		Vect[] rotNe=femCalc.rotNe(jac,zero,edgeDir);

		Edge[] edge=elementEdges(i);
		Vect J=new Vect(dim);
		for(int j=0;j<nElEdge;j++)		{
			J=J.add(rotNe[j].times(edge[j].T));
		}

		element[i].setJ(J);

	}



	public void setJ0(){

		if(this.Tmethod){
			TMatrix tMat=new TMatrix(this);
			this.bcSet.setTBC(this);

			tMat.setTMat(this);

			Vect Ci=this.Cs.scale(this.bT);
			SpMat L=this.Cs.ichol();
			Vect x=this.solver.ICCG(this.Cs, L, this.bT,1e-6,2000).times(1);

			x.timesVoid(Ci);

			for(int i=1;i<=numberOfEdges;i++)
				if(this.T_unknownIndex[i]!=0)
					edge[i].T=x.el[this.T_unknownIndex[i]-1];

			for(int ir=1;ir<=numberOfRegions;ir++){
				if(!region[ir].hasJ) continue;

				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
				{
					this.setElementJ0(i);

				}
			}

			return;

		}


		double Jn2=0,Jmax2=0,Jmin2=0;
		double eps=1e-3;


		Vect regJ=new Vect(3);
		for(int ir=1;ir<=numberOfRegions;ir++){
			if(!region[ir].hasJ) continue;
			regJ=region[ir].getJ();
			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				if(this.coordCode==1 && this.dim==3)
				{

					Mat R=new Mat();


					Mat R2D=util.rotMat2D(util.getAng(this.getElementCenter(i).v2()));
					R=new Mat(dim,dim);
					for(int m=0;m<2;m++)
						for(int n=0;n<2;n++)
							R.el[m][n]=R2D.el[m][n];

					R.el[2][2]=1;

					element[i].setJ(R.mul(regJ));


				}

				else{
					element[i].setJ(regJ);
				}

				Jn2=regJ.dot(regJ);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}
		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);


	}

	public void setM(){

		Vect MReg=new Vect(dim);

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(!region[ir].hasM) continue;


			MReg=region[ir].getM().times(region[ir].getNu());

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setM(MReg);
			}
		}

	}

	public void setFreq(double f){
		this.freq=f;

	}
	public void setDt(double dt){

		this.dt=dt;
	}

	public void setnTsteps(int N){

		this.nTsteps=N;
	}

	public void setDefMode(int defMode){
		this.defMode=defMode;

	}

	public Vect  setDeformation(int step){
		return setDeformation(step,false);
	}


	public Vect  setDeformation(int step,boolean stress){

		boolean dyn=false;

		if(this.Ks==null)
			this.setStiffMat(dyn);



		Vect u=new Vect();
		
		

		if(this.numberOfUnknownU!=0){

			if(this.nonLin){
				
				MatNonAnalysis matNon=new MatNonAnalysis();
				u= matNon.getDeformation(this, solver,defMode,step);
				this.setU(u);	
				if(stress)
					this.setStress();
				return u;
			}
			if(contact!=null) return getDeformationContact(step);

			u=this.mechMat.getDeformation(this, solver,defMode);

			this.setU(u);	

			double un,umax=0; int im=1;
			for(int i=1;i<=this.numberOfNodes;i++){
				if(this.node[i].u!=null ){
					un=this.node[i].u.norm();

					if(un>umax) {umax=un; im=i;}
				}
			}

			System.out.println("Deformation Max: "+umax+" at node "+im);
		}


		if(stress)
			this.setStress();

		return u;


	}
	
	public Vect  setDeformationNonlin(int step,boolean stress){

		boolean dyn=false;

		if(this.Ks==null)
			this.setStiffMat(dyn);


		Vect u=new Vect();

		if(this.numberOfUnknownU!=0){

			if(contact!=null) return getDeformationContact(step);

			u=this.mechMat.getDeformation(this, solver,defMode);

			this.setU(u);	

			double un,umax=0; int im=1;
			for(int i=1;i<=this.numberOfNodes;i++){
				if(this.node[i].u!=null ){
					un=this.node[i].u.norm();

					if(un>umax) {umax=un; im=i;}
				}
			}

			System.out.println("Deformation Max: "+umax+" at node "+im);
		}


		if(stress)
			this.setStress();

		return u;


	}

	public void writeNodalScalar(String file){
		writer.writeNodalScalar(this,file);
	}
	
	public void writePhi(String file){
		writer.writePhi(this,file);
	}

	public Vect  setVibration(int iter){

		Vect u=new Vect();

		boolean dyn=true;

		if(iter==0){
			/*			 {
				  Runtime runtime = Runtime.getRuntime();
				  long totalMemory = runtime.totalMemory(); // current heap allocated to the VM process
				  long freeMemory = runtime.freeMemory(); // out of the current heap, how much is free
				  long maxMemory = runtime.maxMemory(); // Max heap VM can use e.g. Xmx setting
				  long usedMemory = totalMemory - freeMemory; // how much of the current heap the VM is using
				  long availableMemory = maxMemory - usedMemory; // available memory i.e. Maximum heap size minus the current amount used
				  util.pr( availableMemory/1024/1024);
				}*/

			long mx = Runtime.getRuntime().maxMemory();
			util.pr(mx/1024/1024);

			this.setStiffMat(dyn);

			/*			SpMat P=Ks.timesNew(1e-9);
		//	P.showcl();
			int ix=0;
			Vect v=new Vect(100000);
			for (int i=0;i<P.nRow;i++)
				for(int j=0;j<P.row[i].nzLength;j++)
					v.el[ix++]=P.row[i].el[j];
			v.quickSort();
			//v.show();
			Vect v2=new Vect(90);
			ix=0;
			v2.el[ix++]=v.el[0];
			for (int i=1;i<v.length;i++)
				if(abs(v.el[i]-v.el[i-1])>1e-3)
					v2.el[ix++]=v.el[i];
			v2.show();*/



			double a1=this.rayAlpha;
			double a2=this.rayBeta;

			this.Cs=this.Ms.timesNew(a1).addNew(this.Ks.timesNew(a2));

			if(this.timeIntegMode==3)
			{	

				ud=new Vect(Ks.nRow);
				udd=new Vect(Ks.nRow);
			}

		}


		if(this.numberOfUnknownU!=0){
			
			if(contact==null)
			u=this.getVibration(iter);	
			else
			u=this.getVibrationContact(iter);		
			
			this.setU(u);	

			double un,umax=0; int im=1;
			for(int i=1;i<=this.numberOfNodes;i++){
				if(this.node[i].u!=null ){
					un=this.node[i].u.norm();

					if(un>umax) {umax=un; im=i;}
				}
			}

			System.out.println("Deformation Max: "+umax+" at node "+im);

		}

		return u;
	}



	public void setTorque(double r){
		setTorque(new Vect(dim),0,r,1);

	}

	public void setTorque(double r1, double r2){
		setTorque(new Vect(dim),r1,r2,1);
	}

	public void setTorque(double r1, double r2, int mode){

		setTorque(new Vect(dim),r1,r2,mode);

	}


	public void setTorque(Vect origin, double r1,double r2, int mode){

		int nx=0;

		double trq=0;
		Vect F=new Vect();
		for(int i=1;i<=this.numberOfNodes;i++)
		{

			//if(node[i].getMap()>0) {continue;}	
			
			
			Vect r2d=this.node[i].getCoord().v2();
	

			double rn=r2d.norm();
			if(rn<r1 ||rn>r2) continue;

			if(mode>=0) F=this.node[i].getNodalVect(mode);


			if(F!=null){
				nx++;
				double nodeTrq=0;
				if(this.coordCode==1){
					nodeTrq=rn*F.el[1];

				}
				else if(this.coordCode==0){
					nodeTrq=F.v2().cross(r2d).el[2];
				}
				trq=trq+nodeTrq;
			}

		}


		if(this.dim==2)
			this.TrqZ=trq*this.height;
		else this.TrqZ=trq;


		if(this.hasPBC)
			if(this.alpha2-this.alpha1>0){
				this.TrqZ*=2*PI/(alpha2-alpha1);
			}


	}

	public void setTorque(int ir){

		int mode=1;

		int[] nn=this.getRegNodes(ir);
		double trq=0;
		Vect F=new Vect();
		for(int i=0;i<nn.length;i++)
		{
			//if(node[i].getMap()>0) {continue;}		
			int nx=nn[i];
			double rn=this.node[nx].getCoord().v2().norm();

			//if(rn<r1 ||rn>r2) continue;

			if(mode>=0) F=this.node[nx].getNodalVect(mode);


			if(F!=null){

				trq=trq+rn*F.el[1];
			}

		}

		if(this.dim==2)
			this.TrqZ=trq*this.height;
		else this.TrqZ=trq;


		if(this.hasPBC)
			if(this.alpha2-this.alpha1>0)
				this.TrqZ*=2*PI/(alpha2-alpha1);


	}

	public void setJe(){
		double Jn2=0,Jmax2=0,Jmin2=0;
		for(int i=1;i<=numberOfElements;i++){

			if(!element[i].isConductor()) continue;

			setElementJe(i);

			Vect Je=element[i].getJe();

			if(element[i].isConductor()){
				Jn2=Je.dot(Je);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}

		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}
	
	public void setJStatic(){
		double Jn2=0,Jmax2=0,Jmin2=0;
		for(int i=1;i<=numberOfElements;i++){

			if(!element[i].isConductor()) continue;

			setElementJStatic(i);

			Vect Je=element[i].getJe();
			//Je.hshow();
			if(element[i].isConductor()){
				Jn2=Je.dot(Je);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}

		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}


	private void setElementJe(int i){
		if(dim==2) {setElement2DJe(i); return;}

		Node[] vertex=elementNodes(i);
		Edge[] edge=elementEdges(i);
		double[] dAe=new double[nElEdge];
		Vect dA=new Vect(dim);
		for(int j=0;j<nElEdge;j++)
			dAe[j]=edge[j].getDiffA();

		dA=getElementdA(vertex,dAe,i);

		double rdt=1.0/dt;
		
		Vect Je=dA.deepCopy();


		if(analysisMode==2){
		double[] nodePhi=new double[nElVert];
		Vect gradPhi=new Vect(dim);

		for(int j=0;j<nElVert;j++)
			nodePhi[j]=vertex[j].getPhi();
		
		gradPhi=femCalc.gradPhi(vertex,nodePhi);
		
		
		Je=Je.add(gradPhi);
	}
		
		Vect sigma=element[i].getSigma();
				
		Je=Je.times(sigma.times(-rdt));

		element[i].setJe(Je);

	}

	private void setElementJStatic(int i){

		if(this.dim==2){
			setElementJStaticJe2D(i);
			return;
		}
		Node[] vertexNode=elementNodes(i);
		Edge[] elemEdges=elementEdges(i);
		boolean[] edgeDir=element[i].getEdgeReverse();

		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();

			double[] nodePhi=new double[nElVert];
			Vect gradPhi=new Vect(dim);

			for(int j=0;j<nElVert;j++){
				nodePhi[j]=vertexNode[j].getPhi();
			gradPhi=femCalc.gradPhi(vertexNode,nodePhi);
			}
			
			Vect localCo=new Vect(this.dim);
			
			Mat jac=femCalc.jacobian(vertexNode,localCo);

			
			Vect[] Ne=femCalc.Ne(jac,localCo,edgeDir);
			
			Vect Adot=new Vect(dim);
			
			for(int j=0;j<this.nElEdge;j++){
				Adot=Adot.add(Ne[j].times(A[j]));
			
			}
			
			element[i].setJe(gradPhi.add(Adot).times(element[i].getSigma()).times(-1));
			
}

	private void setElement2DJe(int i){

		double rdt=1.0/dt;

		Edge[] edge=elementEdges(i);
		double[] dAe=new double[nElEdge];
		double dA=0;
		for(int j=0;j<nElEdge;j++){
			dAe[j]=edge[j].getDiffA();
		}

		if(elCode==0)
			dA=getElement3angdA(i,dAe);
		else if(elCode==1)
			dA=getElementQuaddA(dAe);

		if(analysisMode==1)
			element[i].setJe(new Vect(0,0,-dA*element[i].getSigmaZ()*rdt));

		else if(analysisMode==2){
			double gradPhiz=0;
			if(elCode==0) gradPhiz=getElement3angPhi(i);
			else if(elCode==1) gradPhiz=getElementQuadPhi(i);
			element[i].setJe(new Vect(0,0,-(dA+gradPhiz)*element[i].getSigmaZ()*rdt));
		}
	}
	
	private void setElementJStaticJe2D(int i){

		double rdt=1.0/this.height;

		Edge[] edge=elementEdges(i);

		double[] Ae=new double[nElEdge];
		double A=0;

		if(elCode==0)
			A=getElement3angA(i);
		else if(elCode==1)
			A=getElementQuadA(i);
	
		
			double gradPhiz=0;
			if(elCode==0) gradPhiz=getElement3angPhi(i);
			else if(elCode==1) gradPhiz=getElementQuadPhi(i);
			element[i].setJe(new Vect(0,0,-(A+gradPhiz)*element[i].getSigmaZ()*rdt));
			
		}	


	public Vect getElementA(int ie){

		boolean[] edgeDir=element[ie].getEdgeReverse();

		Vect A=new Vect(dim);
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		Edge[] edge=this.elementEdges(ie);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		Vect[] Ne=femCalc.Ne(jac,zero,edgeDir);

		for(int j=0;j<nElEdge;j++)	{		
			A= A.add(Ne[j].times(edge[j].A));
		}
		return  A;	

	}

	public double getElementPhi(int ie){

		double phi=0;
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		double[] N=femCalc.N(zero);

		for(int j=0;j<nElVert;j++)	{		
			phi+=N[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}



	public double getElementT(int ie){

		double T=0;
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		double[] N=femCalc.N(zero);

		for(int j=0;j<nElVert;j++)	{		
			T+=N[j]*vertexNode[j].T;
		}
		return  T;	

	}

	public Vect getElementdA(Node[] vertexNode,double[] Ae,int i){

		boolean[] edgeDir=element[i].getEdgeReverse();

		Vect dA=new Vect(dim);
		Vect zero=new Vect(dim);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		Vect[] Ne=femCalc.Ne(jac,zero,edgeDir);

		for(int j=0;j<nElEdge;j++)	{		
			dA= dA.add(Ne[j].times(Ae[j]));
		}
		return  dA;	

	}
	public double getElementQuadA(double[] dAe){

		double dA=0;
		Vect zero=new Vect(2);
		double[] Ne=femCalc.NeQuad(zero);
		for(int j=0;j<nElEdge;j++)	{		
			dA+=Ne[j]*dAe[j];
		}
		return  dA;	

	}


	public double getElementQuadPhi(int ie){

		double phi=0;
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		double[] N=femCalc.NQuad(zero);

		for(int j=0;j<nElVert;j++)	{		
			phi+=N[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}

	public double getElementQuaddA(double[] dAe){

		double dA=0;
		Vect zero=new Vect(2);
		double[] Ne=femCalc.NeQuad(zero);
		for(int j=0;j<nElEdge;j++)	{		
			dA+=Ne[j]*dAe[j];
		}
		return  dA;	

	}

	public double getElement2DA(int ie){
		double A=0;
		if(this.elCode==0)
			A= getElement3angA(ie);
		else if(this.elCode==1)
			A= getElementQuadA(ie);


		return  A;	

	}

	public double getElement3angA(int ie){

		int[] en=this.element[ie].getEdgeNumb();

		double a=1.0/3;
		double[] localNe={a,a,a};
		double A=0;
		for(int j=0;j<3;j++)	{		
			double Aj=0;
			if(this.edge[en[j]].map==0)
				Aj=this.edge[en[j]].A;
			else
				Aj=this.edge[this.edge[en[j]].map].A;

			A+=localNe[j]*Aj;
		}
		return  A;	

	}


	public double getElement3angPhi(int ie){

		double phi=0;
		Node[] vertexNode=this.elementNodes(ie);
		double a=1.0/3;
		double[] localNe={a,a,a};
		for(int j=0;j<nElVert;j++)	{		
			phi+=localNe[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}

	public double getElement3angdA(int ie,double[] dAe){

		double a=1.0/3;
		double[] localNe={a,a,a};
		double dA=0;
		for(int j=0;j<3;j++)	{		
			dA+=localNe[j]*dAe[j];
		}
		return  dA;	

	}

	public void setElementsParam(){



		int nx=0,ny=0;
		for(int ir=1;ir<=this.numberOfRegions;ir++){
			if(this.region[ir].isNonLinear){				
				nx++;
				this.region[ir].BHnumber=nx;
			}else{
				this.region[ir].BHnumber=0;
				}

			if(this.region[ir].MS){				
				ny++;
				this.region[ir].lamBNumber=ny;
			}

		}
		this.nBH=nx;
		this.nLam=ny;

		if(coupled){
			this.BHS=new BHSCurve[this.nBH+1];
			this.lamBS=new LamBSCurve[this.nLam+1];		
		}
		else{
			this.BH=new BHCurve[this.nBH+1];
			this.lamB=new LamBCurve[this.nLam+1];
		}

		try {

			for(int ir=1;ir<=this.numberOfRegions;ir++){

				if(this.region[ir].lamBNumber>0)		
					if(coupled){
						this.lamBS[this.region[ir].lamBNumber]=new LamBSCurve() ;
						/*	Curve cv=new Curve(this.lamBS[this.region[ir].lamBNumber],600,600);
						cv.show(true);*/
					}
					else{
						this.lamB[this.region[ir].lamBNumber]=new LamBCurve(this.region[ir].getMaterial()) ;
					}


				if(this.region[ir].BHnumber>0)	
					if(coupled){
						if(cpvms && calCurve)
							this.BHS[this.region[ir].BHnumber]=new BHSCurve(this, ir, this.BHstress,this.region[ir].lamBNumber);
						else
							this.BHS[this.region[ir].BHnumber]=new BHSCurve(this.region[ir].getMaterial(),cpvms);

						/*	String file = System.getProperty("user.dir") + "\\BHSH350.xls";
						this.BHS[this.region[ir].BHnumber].writexls(file);*/


						/*				String file2 = System.getProperty("user.dir") + "\\LamBSbelahcen.xls";
						this.lamBS[this.region[ir].lamBNumber].writexls(file2);*/

						//util.show(this.BHS[this.region[ir].BHnumber].BH[this.BHS[this.region[ir].BHnumber].nZeroStress].BH);
						if(ir==1){/*
						this.BH[this.region[ir].BHnumber].showCurve();


							Curve cv=new Curve(this.BHS[this.region[ir].BHnumber],600,600);
							cv.show(true);
						 */}
					}
					else{
						this.BH[this.region[ir].BHnumber]=new BHCurve(this.region[ir].getMaterial()) ;
						double nux=this.BH[this.region[ir].BHnumber].getNu(1.0);							
						region[ir].setNu(nux);

					}


			}


		} catch (Exception e) {
			System.err.println("Can't open magnetization or magnetostriction data file.");
		}


		setM();

		nEdgeHasJ=0;


		for(int ir=1;ir<=numberOfRegions;ir++){

			boolean regCond=region[ir].isConductor;

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setRo(region[ir].getRo());
				element[i].setNu(region[ir].getNu());
				element[i].setSigma(region[ir].getSigma());

				if(regCond)
					element[i].setSigma(region[ir].getSigma());

				element[i].setRegion(ir);

				element[i].setYng(region[ir].getYng());

				element[i].setPois(region[ir].getPois());

				element[i].setShear(region[ir].getShear());


			}
		}

		for(int ir=1;ir<=numberOfRegions;ir++)
			if(this.region[ir].rotor)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
					element[i].rotor=true;		

		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}


	public void setElementsParamMech(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setRo(region[ir].getRo());
				element[i].setRegion(ir);
				element[i].setYng(region[ir].getYng());
				element[i].setPois(region[ir].getPois());
				element[i].setShear(region[ir].getShear());

				element[i].setDeltaT(region[ir].getDeltaT());

			}
		}

		for(int ir=1;ir<=numberOfRegions;ir++)
			if(this.region[ir].rotor)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
					element[i].rotor=true;		

		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}

	public void setElementsParamSeep(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setSigma(region[ir].getSigma());

			}
		}


		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}


	public void setElementsParamPC(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setSigma(region[ir].getEr());
				element[i].setNu(region[ir].getMur());

			}
		}


		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}



	public void scaleKnownEdgeAL(double c){
		for(int i=1;i<=numberOfKnownEdges;i++)
			if(knownEdgeValue[i]!=0)
				edge[knownEdgeNumber[i]].setSolvedAL(knownEdgeValue[i]*c);	
	}


	public void saveAp(){
		for(int i=1;i<=this.numberOfEdges;i++){
			edge[i].saveAp();	


		}

	}



	public Vect solveMagLin(int step,Vect x_init){

		return femSolver.solveMagLin(this,step,x_init);
	}

	public Vect solveNonLinear(Vect x, boolean b,int step){
		
		return femSolver.solveMagNonlin(this,x, b,step);
	}


	public void solveCoupled(Vect x){

		femSolver.solveCoupled(this,x);
	}


	public void setNodePhi(Vect x){
		int nodeNumber;
		util.pr("-------------------------  "+x.norm());

		for(int i=1;i<=numberOfVarNodes;i++){
			nodeNumber=varNodeNumber[i];	
			if(nodeNumber>0){
				node[varNodeNumber[i]].setPhi(x.el[i+numberOfUnknownEdges-1]);	
			}
		}

	}

	public void setSolution(Vect x){

		setSolvedAL(x);

		if(analysisMode==2)
			setNodePhi(x);
		
		if(this.Q!=null){

			int kp=this.Q.nCol;
			Vect vp=new Vect(kp);
			for(int k=0;k<vp.length;k++)
				vp.el[k]=x.el[x.length-vp.length+k];

		Vect interfaceA=this.Q.mul(vp);
		
		for(int i=1;i<=this.numberOfEdges;i++){
			
			if(this.edgeOnFSIndices[i]>=0){	
				this.edge[i].setA(interfaceA.el[this.edgeOnFSIndices[i]]);
			}
		}

		}
		
/*		for(int i=1;i<=this.numberOfEdges;i++){
			Vect v1=this.edge[i].node[0].getCoord();
			Vect v2=this.edge[i].node[1].getCoord();
			Vect edv=v2.sub(v1);
			double xx=v1.add(v2).el[0]/2;
			Vect A=new Vect(0,0,xx);
		//	double a=A.dot(edv);
			//this.edge[i].setA(a);
			util.pr(i+" , "+this.edge[i].A);

		}*/
			
		this.setB();	

	}

	public Vect getUnknownA(){
		Vect x=new Vect(this.numberOfUnknownEdges);
		for(int i=0;i<x.length;i++){
			x.el[i]=edge[unknownEdgeNumber[i+1]].A;	
		}

		return x;
	}

	public Vect getUnknownAp(){
		Vect x=new Vect(this.numberOfUnknownEdges);
		for(int i=0;i<x.length;i++)
			x.el[i]=edge[unknownEdgeNumber[i+1]].Ap;	

		return x;
	}

	public int[] getRegNodes(int ir){

		boolean[] nc=new boolean[1+this.numberOfNodes];
		int[] nn=new int[this.numberOfNodes];
		int ix=0;
		for(int i=this.region[ir].getFirstEl();i<=this.region[ir].getLastEl();i++){		
			int[] vertNumb=this.element[i].getVertNumb();
			for(int j=0;j<nElVert;j++){
				int nx=vertNumb[j];
				if(!nc[nx]){

					nc[nx]=true;
					nn[ix]=nx;
					ix++;
				}
			}
		}

		int[] regNodes=new int[ix];
		for(int i=0;i<ix;i++)
			regNodes[i]=nn[i];

		return regNodes;


	}
	

	public void setB(){

		double Bn2,Bmax2=0,Bmin2=0;

		for(int i=1;i<=numberOfElements;i++){
			
			setElementB(i);

			Bn2=element[i].getB().dot(element[i].getB());
			if(Bn2>Bmax2)
				Bmax2=Bn2;
			if(Bn2<Bmin2)
				Bmin2=Bn2;}

		Bmax=sqrt(Bmax2);
		Bmin=sqrt(Bmin2);



	}


	public void setHead(Vect h){
		for(int i=1;i<=this.numberOfNodes;i++){
			if(this.T_unknownIndex[i]!=0)
			{
				this.node[i].T=h.el[this.T_unknownIndex[i]-1];

			}
		}
	}
	public void setVelocity(){

		double Bn2,Bmax2=0,Bmin2=0;

		for(int i=1;i<=numberOfElements;i++){

			if (!element[i].isConductor()) continue;

			setElementVelocity(i);

			Bn2=element[i].getB().dot(element[i].getB());
			if(Bn2>Bmax2)
				Bmax2=Bn2;
			if(Bn2<Bmin2)
				Bmin2=Bn2;}

		Bmax=sqrt(Bmax2);
		Bmin=sqrt(Bmin2);


	}


	public void setElementVelocity(int i){
		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] gradN;
		Vect B=new Vect(2);
		Vect Kt=element[i].getSigma();

		jac=femCalc.jacobian(vertexNode,zero);
		gradN=femCalc.gradN(jac,zero);

		int[] vn=element[i].getVertNumb();

		for(int j=0;j<this.nElVert;j++){
			double T=node[vn[j]].T;
			B=B.add(gradN[j].times(Kt).times(-T));
		}

		element[i].setB(B);

	}

	public Vect getBAt(Vect P){
		return femCalc.getBAt(this, P);
	}

	public double[] getApAnAt(Vect P){

		return femCalc.getApAnAt(this, P);
	}

	public Vect getStressAt(Vect P){
		if(this.elCode==3) return new Vect(dim);
		Vect na=new Vect(3*(this.dim-1));
		na.el[0]=1e10;
		na.el[1]=-1e10;
		int[] m=this.getContainingElement(P);
		if(m[0]<=0) return na;

		if(this.elCode==0) return this.element[m[0]].getStress();

		Vect lc=this.femCalc.localCo(this,m,P);

		return getStress(m[0]);

	}

	public Vect[] getAllB(){	

		Vect[] B=new Vect[numberOfElements];

		for(int ie=0;ie<numberOfElements;ie++)
			B[ie]=element[ie+1].getB();
		return B;
	}

	public double getDiffMax(Vect[] u, Vect[] v){
		double diff;
		double diffMax=0;
		for(int i=0;i<u.length;i++){


			diff=u[i].sub(v[i]).norm();
			if(diff>diffMax)
				diffMax=diff;
		}

		return diffMax;
	}

	public double getBmax(){	

		double Bmax=0;

		for(int i=1;i<=numberOfElements;i++){
			double Bn=element[i].getB().norm();
			if(Bn>Bmax) Bmax=Bn;

		}
		return Bmax;
	}
	
	public double getFluxErrSquared(Vect[] u, Vect[] v){	

		double err_sqared=getErrorSquared(u,v)/getSumB_Squared(u,v);
		
		return err_sqared;
	}
	
	public double getErrorSquared(Vect[] u, Vect[] v){	

		double sum=0;

		for(int i=0;i<u.length;i++){


			sum+=u[i].sub(v[i]).norm2();
		
		}
		return sum;
	}

	public double getSumB_Squared(Vect[] u, Vect[] v){	

		double sum=0;
		for(int i=0;i<u.length;i++){

			double Bn2=u[i].norm2()+v[i].norm2();
			sum+=Bn2;

		}
		return sum;
	}

	public double getSumB_Squared(){	

		double sum=0;

		for(int i=1;i<=numberOfElements;i++){
			double Bn2=element[i].getB().norm2();
			sum+=Bn2;

		}
		return sum;
	}

	public double getRmax(){	

		double rmax=0;

		for(int i=1;i<=numberOfNodes;i++){

			double rn=node[i].getCoord().v2().norm();
			if(rn>rmax) rmax=rn;

		}
		return rmax;
	}

	public void setAuMax(){	

		AuMax=0;

		for(int i=1;i<=numberOfEdges;i++){
			double aum=abs(edge[i].A);
			if(aum>AuMax) AuMax=aum;

		}

	}

	public void setuMax(){	

		uMax=0;

		for(int i=1;i<=numberOfNodes;i++){
			if(node[i].u==null) continue;
			double a=node[i].u.norm();
			if(a>uMax) uMax=a;

		}

	}

	public Vect[] getAllU(){	

		Vect[] u=new Vect[this.numberOfUnknownU];

		for(int i=0;i<this.numberOfUnknownU;i++){

			u[i]=node[unknownUnumber[i+1]].u.deepCopy();
		}
		return u;
	}


	public void setJ(double t)
	{
		setJ(t,1.0);
	}

	public void setJ(double t,double factor){


		for(int ir=1;ir<=numberOfRegions;ir++){
			if(!region[ir].circuit){
				if(region[ir].hasJ  ) {
					setJfJ(ir,t);	
				}
				else {
					//setJfV(ir,t,factor);	<==== temporary fir J given
				}

			}
		}


		if(this.threePhaseRegs[0]!=0 && this.threePhaseRegs[1]!=0 && this.threePhaseRegs[2]!=0)
			setJ3phaseStrandedCircuit(t);
		else if(this.circuit)
			setJStrandedCircuit(t);

	}

	public void setJfJ(int ir,double t){
		Vect J=new Vect();

		J=this.region[ir].getJ().times(cos(this.region[ir].omega*t+this.region[ir].phase0+this.region[ir].beta));
		//============

		double ju=0;



		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

			this.element[i].setJ(J);

		}
	}



	public void setdJfromCurrentWaves(double t,double t1, double t2,double alpha,double factor){
		Vect J=new Vect();

		for(int ir=1;ir<=numberOfRegions;ir++){
	
			

				 double curr=0;

			if(ir==9 ) curr=this.ia.getI(t)-(this.ia.getI(t1)*(1-alpha)+this.ia.getI(t2)*alpha);
			else if(ir==10 ) curr=-(this.ia.getI(t)-(this.ia.getI(t1)*(1-alpha)+this.ia.getI(t2)*alpha));
			else if(ir==11 ) curr=this.ib.getI(t)-(this.ib.getI(t1)*(1-alpha)+this.ib.getI(t2)*alpha);
			else if(ir==12 ) curr=-(this.ib.getI(t)-(this.ib.getI(t1)*(1-alpha)+this.ib.getI(t2)*alpha));
			else if(ir==13 ) curr=this.ic.getI(t)-(this.ic.getI(t1)*(1-alpha)+this.ic.getI(t2)*alpha);
			else if(ir==14 ) curr=-(this.ic.getI(t)-(this.ic.getI(t1)*(1-alpha)+this.ic.getI(t2)*alpha));



			double jz=curr*this.region[ir].NtS*factor;
			//if(this.dim==2)
			J=new Vect(0,0,jz);
	/*		else
				J=new Vect(0,jz,0);*/
	//	}

		
		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){
	

			this.element[i].setJ(J);
		}
		
	}
		
		}
	
	public void setJfV(int ir,double t)
	{
		setJfV(ir,t,1);
	}
	
	public void setJfV(int ir,double t,double factor){
		Vect J=new Vect();
		//if(this.dim==2 ){
			double curr=this.region[ir].terminalVoltage0/this.region[ir].getWireRes();
			
			double kf=0;
			
			kf=cos(this.region[ir].omega*t+this.region[ir].phase0+this.region[ir].beta);
	
			
			curr*=kf;
			
			
			
if(this.numberOfRegions>8){
				curr=0;

			if(ir==9 ) curr=this.ia.getI(t);
			else if(ir==10 ) curr=-this.ia.getI(t);
			else if(ir==11 ) curr=this.ib.getI(t);
			else if(ir==12 ) curr=-this.ib.getI(t);
			else if(ir==13 ) curr=this.ic.getI(t);
			else if(ir==14 ) curr=-this.ic.getI(t);

		}

curr*=factor;//*(1+.2*(.5-Math.random()));

			double jz=curr*this.region[ir].NtS;
			//if(this.dim==2)
			J=new Vect(0,0,jz);
	/*		else
				J=new Vect(0,jz,0);*/
	//	}

		
		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){
	

			this.element[i].setJ(J);
		}
		
		}


	public void setJ3phaseStrandedCircuit(double t){


		if(dim==3){
			throw new IllegalArgumentException("stranded 3D coil not defined.");
		}

		double fc=7200;
		//	fc=5450;
		for(int j=0;j<3;j++){
			int ir=this.threePhaseRegs[j];
			double kf=0;
			double vr=spwmLevel*cos(region[ir].phase0+this.region[ir].beta+region[ir].omega*t);

			double vt=util.triangWave(1,fc*t);

			if(vr>vt)
				kf=1;

			else

				kf=0;

			//	kf*=0.5;

			//	kf=vr;

			region[ir].terminalVoltagep=region[ir].terminalVoltage;
			region[ir].terminalVoltage=kf*this.region[ir].terminalVoltage0;

			/*		if(ir==this.threePhaseRegs[0])
			region[ir].terminalVoltage=this.va.getI(t);
			else if(ir==this.threePhaseRegs[1])
				region[ir].terminalVoltage=this.vb.getI(t);
			else if(ir==this.threePhaseRegs[2])
				region[ir].terminalVoltage=this.vc.getI(t);*/

		}			

	}

	public void setJStrandedCircuit(double t){


		/*if(dim==3){
			throw new IllegalArgumentException("stranded 3D coil not defined.");
			}*/

		for(int i=0;i<unCurRegNumb.length;i++){
			int ir=this.unCurRegNumb[i];

			double kf=cos(region[ir].phase0+this.region[ir].beta+region[ir].omega*t);

			/*	if(t<.03)
		kf=200*t;
		else kf=6;*/

			//kf=1-exp(-1000*t);

			region[ir].terminalVoltagep=region[ir].terminalVoltage;
			region[ir].terminalVoltage=kf*this.region[ir].terminalVoltage0;

		}

		/*		int ir=this.unCurRegNumb[0];

			//*****************
			region[ir].terminalVoltagep=region[ir].terminalVoltage;
			region[ir].terminalVoltage=this.va.getI(t);*/



	}


	public void setBeta(double beta){

		for(int ir=1;ir<=numberOfRegions;ir++){

			if(region[ir].hasJ) {
				region[ir].setBeta(beta);
			}

		}

	}

	public Vect getFOf(Vect globalCo){

		int m[]=getContainingElement(globalCo);

		if(m[0]<=0) throw new NullPointerException("given point outside the space ");

		if(element[m[0]].getNu().norm()==nu0)
			return new Vect(dim);
		Vect F=new Vect(dim);
		for(int j=0;j<nElVert;j++)
			F=F.add(node[element[m[0]].getVertNumb(j)].F);

		return  F.times(1.0/nElVert);	

	}

	public Vect getJeAt(Vect globalCo){

		int[] m=getContainingElement(globalCo);
		if(m[0]<=0) throw new NullPointerException("given point outside the space ");
		else if(!this.element[m[0]].isConductor()) return new Vect(3);
		return this.element[m[0]].getJe();

	}

	/*	public Vect getJeAt(Vect globalCo){

		int[] m=getContainingElement(globalCo);

		if(m[0]<=0) throw new NullPointerException("given point outside the space ");

		if(!element[m[0]].isConductor())
			return new Vect(3);

		Node[] vertex=elementNodes(m[0]);
		Edge[] edge=elementEdges(m[0]);
		double[] dAe=new double[nElEdge];

		for(int j=0;j<nElEdge;j++)
			dAe[j]=edge[j].getDiffAu();


		if(elCode==0){
			dA=getElement3angA(m[0]);

		Vect dA=getElementA(vertex,dAe);
		if(analysisMode==1)
		return dA.times(element[m[0]].getSigma().times(-rdt));
		else if(analysisMode==2){

			double[] nodePhi=new double[nElVert];
			Vect gradPhi=new Vect(dim);
			for(int j=0;j<nElVert;j++)
				nodePhi[j]=vertex[j].getPhi();
			gradPhi=femCalc.gradPhi(vertex,nodePhi);
			return (dA.times(rdt).add(gradPhi)).times(element[m[0]].getSigma().times(-1));

		}
		else
			return new Vect(dim);


	}*/


	public Vect getJeOf(int m){
		if(!element[m].isConductor())
			return new Vect(3);

		Vect Jn=new Vect(3);
		if(m>0)
			Jn=element[m].getJe();

		else if(m<-1)
			Jn=element[-m-1].getJe();

		return Jn;
	}

	public void writeB(String file){
		writer.writeB(this,file);

	}


	public void writeA(String file){
		writer.writeA(this, file);

	}
	public void writeNodalA2D(String file){
		

		for(int n=1;n<=this.numberOfEdges;n++)
		{	
			this.edge[n].node[0].T=this.edge[n].A;
		}		
		
		writer.writeNodalScalar(this, file);

	}


	public void writeNodalField(String nodalForceFile,int mode){

		writer.writeNodalField(this, nodalForceFile, mode);
	}

	public void appendNodalField(String nodalForceFile,int mode,int[]nn, int it){

		writer.appendNodalField(this, nodalForceFile, mode,nn,it);
	}

	public void appendScalarField(String file,int mode,int[]nn, int it){

		writer.appendScalarField(this, file, mode,nn,it);
	}

	public void writeJe(String eddyFile){
		writer.writeJe(this, eddyFile);

	}

	public void writeJ0(String j0filr){
		writer.writeJ0(this, j0filr);

	}


	public void writeStress(String stressFile){

		writer.writeStress(this, stressFile);

	}


	public void rotate(double rad){
		MeshManipulator mf=new MeshManipulator();
		mf.rotate(this,rad,false);	
	}


	public void setRotorPosition(double deg){
		Combiner cmb=new Combiner();
		cmb.setRotorPosition(this,deg);	
	}

	public void setRotorIndex(int step){
		Combiner cmb=new Combiner();
		cmb.setRotorIndex(this,step);	
	}

	public void combine(String rotf, String statf,int ang){
		combiner.combineTwoNode(this,rotf,statf,ang);	
	}

	public void combineFull(String rotf, String statf,int ang){
		combiner.combineFull(this,rotf,statf,ang);	
	}

	private double[] subS3ang(int ie, Vect P){

		Node[] vertexNode=elementNodes(ie);
		Vect v1=vertexNode[0].getCoord().sub(P);
		Vect v2=vertexNode[1].getCoord().sub(P);
		Vect v3=vertexNode[2].getCoord().sub(P);

		double[] N=new double[3];

		N[0]=v2.cross(v3).norm()/2;
		N[1]=v3.cross(v1).norm()/2;
		N[2]=v1.cross(v2).norm()/2;
		return N;
	}


	public Vect getElementCenter(int ie){


		Vect center=new Vect(dim);
		int[] vertNumb=element[ie].getVertNumb();
		for(int j=0;j<nElVert;j++)
			center=center.add(node[vertNumb[j]].getCoord());


		return center.times(1.0/nElVert);
	}

	public void resetMotor(){


		for(int i=1;i<=numberOfNodes;i++){
			if(node[i].getMap()>0){
				node[i].common=false;
				node[i].sPBC=false;
				node[i].aPBC=false;
				node[i].setMap(0);
			}
		}

		for(int i=1;i<=numberOfEdges;i++){
			if(edge[i].map>0){
				edge[i].common=false;
				edge[i].sPBC=false;
				edge[i].aPBC=false;
				edge[i].map=0;
			}
		}

		int L=0;
		if(commonNodes!=null){
			L=commonNodes[0][0].length;

		for(int iz=0;iz<commonNodes.length;iz++)
		for(int i=0;i<L;i++){
			int nn=commonNodes[iz][0][i];
			node[nn].setMap(commonNodes[iz][1][i]);
			node[nn].common=true;
		}
		}

		if(this.hasPBC)
			this.mapPBC();
	}


	public void setStress(){
		setStress(true);
	}

	public void setStress(boolean add){
		stressMax=0;stressMin=1e40;
		for(int i=1;i<=this.numberOfElements;i++)
			if(element[i].isDeformable()){
				setElementStress(i,add);
				double sn=element[i].getStress().norm();
				if(sn<stressMin) stressMin=sn;
				if(sn>stressMax) stressMax=sn;
			}

	}




	public void setElementStress(int ie,boolean add){


		Vect stress=this.getStress(ie, new Vect(dim));

		if(!add || element[ie].getStress()==null )
			element[ie].setStress(stress);
		else
			element[ie].setStress(stress.add(element[ie].getStress()));


	}



	public Vect getStress(int ie,Vect lc){

		boolean MS=((this.element[ie].hasMS() || this.element[ie].isThermal()) && (this.defMode>1));

		Vect strain= femCalc.getStrain(this,ie,lc);

		Mat D =new Mat();
		if(this.dim==3)
			D=femCalc.hook3D(this,ie);
		else
			D=femCalc.hook(this,ie);

		Vect stress=D.mul(strain);

		stress=stress.times(1e-6);


		if(MS)
		{
			stress=stress.add(element[ie].getStress());

		}

		return stress;


	}

	public Vect getStress(int ie){

		return element[ie].getStress();


	}

	public void setNodalStress(){
		double c=1.0/nElVert;
		double eps=1.0;
		int[] count=new int[this.numberOfNodes+1];
		double[] str=new double[this.numberOfNodes+1];
		stressMax=-1e40;
		stressMin=1e40;
		Mat S;
		Vect sv=null;
		Eigen eg=new Eigen();
		double se=0;
		for(int i=1;i<=this.numberOfElements;i++){

			//se=element[i].getStress().norm();

			sv=element[i].getStress();

			if(sv!=null){


				S=util.tensorize(sv);
				double s1,s2,s12,s3;
				s1=sv.el[0];
				s2=sv.el[1];
				s12=sv.el[2];
				s3=element[i].getPois().el[0]*(s1+s2);

				se=pow(s1-s2,2)+pow(s2-s3,2)+pow(s1-s3,2)+6*s12;
				se/=2;
				se=sqrt(se);

			}

			//double se=element[i].getStress().sum();
			//	double se=element[i].getStress().el[2];
			if(abs(se)<eps) continue;
			if(se>stressMax)stressMax=se;
			if(se<stressMin)stressMin=se;
			int[] vertNumb=element[i].getVertNumb();
			for(int j=0;j<nElVert;j++){
				count[vertNumb[j]]++;
				str[vertNumb[j]]+=c*se;
			}

		}

		nodalStressMax=-1e40;
		nodalStressMin=1e40;

		for(int i=1;i<=this.numberOfNodes;i++)
		{
			if(count[i]>0)
				str[i]/=count[i];


			node[i].stress=str[i]*nElVert;
			if(node[i].stress>nodalStressMax)
				nodalStressMax=node[i].stress;
			if(node[i].stress<nodalStressMin)
				nodalStressMin=node[i].stress;
		}

		//nodalStressMin=nodalStressMax*.1;

		util.pr("Stress Min: "+this.stressMin+"  Stress Max: "+this.stressMax);
		util.pr("Nodal Stress Min: "+this.nodalStressMin+"  Nodal Stress Max: "+this.nodalStressMax);
	}

	public void resultAt(Vect P){
		System.out.println();
		Vect Bn=this.getBAt(P);
		if(Bn.el[0]==1e10 && Bn.el[1]==-1e10) util.pr(" >>>>  Given point is located outside the analyzed space!");
		{
			System.out.print("At : ");	P.hshow();
			System.out.print("B(muT) : ");	Bn.times(1e6).hshow();

			System.out.println();

			if(this.analysisMode>0){

				Vect Jen=this.getJeAt(P);
				System.out.print("J : ");	Jen.hshow();
				System.out.println();
				System.out.println();
			}
		}
	}

	public double getSliceAngle(){
		return this.alpha2-this.alpha1;
	}

	public double getStressB(int ie){
		return getStressB(ie,this.element[ie].getB());
	}

	public double getStressB(int ie,Vect B){

		double Bn=B.norm();
		if(Bn==0)
			return 0;

		Vect dirB=B.times(1.0/Bn);
		Mat S=this.element[ie].getStressTensor();
		return S.mul(dirB).dot(dirB);	

	}

	public Vect getStressToB(int ie,Vect B){

		double Bn=B.norm();
		if(Bn==0)
			return new Vect(dim);

		Vect dirB=B.times(1.0/Bn);
		Mat S=this.element[ie].getStressTensor();
		return S.mul(dirB);	

	}


	public void modalQR(){

		this.setStiffMat(true);

		Eigen eg=new Eigen();
		eg.generalSym(Ks.matForm(true),Ms.matForm(true));

		//Mat Q=eg.V;
		Vect lam=eg.lam;


		Vect lam2=lam.sqrt().times(.5/PI);
		this.lams=lam2.deepCopy();
		this.eigVects=eg.V;


	}

	public void modalAnalysis(int p, double er)
	{


		modalAnalysis(p,er,0.0,false);
	}


	public void modalAnalysis(int p, double er,boolean write)
	{


		modalAnalysis(p,er,0.0,write);
	}

	public void modalAnalysis(int p, double er, double f0,boolean write){





		if(this.numberOfUnknownUcomp<200){
			modalQR();
			return;
		}


		this.setStiffMat(true);


		double shft=pow(2*PI*f0,2);
		Eigen eg=new Eigen();
		SpMat K1=Ks.deepCopy();

		//----------
		if(p>K1.nRow) p=K1.nRow;
		//----------

		double cm=1e-8;
		cm=1;
		SpMat M1=Ms.deepCopy();

		M1.times(shft);

		K1.addSmaller(M1);

		K1.times(cm);

		M1=null;

		Mat Q=new Mat(K1.nRow,p);
		Vect lam=eg.subspace(K1, Ms,p,Q, er,this.iterMax, solver);

		lam=lam.times(1.0/cm);

		lam=lam.add(-shft);

		for(int j=0; j<lam.length;j++)
			if(lam.el[j]<0) lam.el[j]=0;

		Vect lam2=lam.sqrt().times(.5/PI);

		this.lams=lam2.deepCopy();
		this.lams.show();
		//lam2.hshow();
		/*Mat LV=new Mat(Q.nRow+1,p);
		LV.setRow(lam, 0);
		for(int i=1; i<LV.nRow;i++)
			LV.setRow(Q.rowVect(i-1), i);*/

		this.eigVects=Q.deepCopy();

		boolean writeShape=false;

		if(write)
			for(int i=0;i<lam.length;i++){

				Vect q=Q.getColVect(i);
				this.setU(q);
				String nodalfile=this.resultFolder+"\\modalShape"+i+".txt";

				double kk=.5*this.maxEdgeLength/q.abs().max();
				kk=.1*this.getRmax()/q.abs().max();

				if(writeShape)
					this.writeModeShape(nodalfile, kk);
				else
					this.writeNodalField(nodalfile, -1);

			}


		String qfile=this.resultFolder+"\\eigVects.txt";
		Q.normalizeColumns();
		writer.writeMat(Q,qfile);


		//return LV;

	}

	public void modalAnalysisAll(){


		this.setStiffMat(true);
		Eigen eg=new Eigen();

		eg.general(this.Ks.matForm(), this.Ms.matForm());
		eg.lam.show();

	}



	public Vect getbUt(int mode){

		Vect bUt=new Vect(this.numberOfUnknownUcomp);
		int[][] index=new int[1+this.numberOfNodes][this.dim];
		Mat R=new Mat();
		if(hasPBC){
			if(this.coordCode==1){
				if(this.dim==2)
					R=util.rotMat2D(this.cpm);		
				else{
					Mat R2D=util.rotMat2D(this.cpm);
					R=new Mat(this.dim,this.dim);
					for(int m=0;m<2;m++)
						for(int n=0;n<2;n++)
							R.el[m][n]=R2D.el[m][n];

					R.el[2][2]=1;
				}
			}
		}
		else 
			R.eye(this.dim);

		int ix=0,nn;
		Vect F,Fp5;

		if(mode>0){

			for(int i=1;i<=this.numberOfNodes;i++){

				if(!this.node[i].isDeformable() || 
						this.node[i].is_U_known() || 			
						this.node[i].getMap()>0) continue;

				for(int p=0;p<this.dim;p++)
					if(!this.node[i].is_U_known(p)){

						index[i][p]=ix++;
					}
			}


			for(int i=1;i<=this.numberOfNodes;i++){
				if(!this.node[i].isDeformable() || this.node[i].is_U_known()) continue;

				nn=i;

				F=this.node[i].getNodalVect(mode);

				if(F==null) continue;

				Vect Fstr=this.node[i].Fstr;


				if(this.node[i].hasPBC())
				{
					nn=this.node[i].getMap();
					F=R.mul(F);
				}

				if(!this.node[i].hasPBC()){
					for(int p=0;p<this.dim;p++)
						if(!this.node[i].is_U_known(p)){
							bUt.el[index[nn][p]]+=F.el[p];
							if(this.timeIntegMode==4 && Fstr!=null ){
								bUt.el[index[nn][p]]+=Fstr.el[p];
								bUt.el[index[nn][p]]*=0.5;

							}


						}

				}
			}

		}

		return bUt;
	}




	public void transfer2DTo3D(String file,  boolean flux3DNeeded,boolean force3DNeeded){

		//m2d.loadNodalField(file,1);


		// the following part enforces the periodic bc from 2d to the 3d model 
		for(int j=1;j<=m2d.numberOfNodes;j++){
			Vect v=m2d.node[j].getCoord();
			if(m2d.node[j].F!=null && v.el[0]<1e-5){

				m2d.node[j].F=new Vect(0,0);
			}
		}

		for(int i=1;i<=this.numberOfNodes;i++){
			this.node[i].F=null;
		}


		int ir=8;
		int p=0;
		int nr=1;
		int nL=6;
		int nP=1;

		int nLh=nL;
		/*	if(this.tag==17 ) {nL=6; nLh=nL;nP=4;}
	else
		if(this.tag==19 )*/ {nL=6; nLh=nL/2;}

		//=================
		if(this.numberOfElements<10000) 
		{nLh=3;
		nP=1;

		}
		else if(this.numberOfElements<20000) 
		{nLh=6;
		nP=1;
		}
		else if(this.numberOfElements<50000) 
		{nLh=6;
		nP=4;
		}

		else if(this.numberOfElements>60000){
			nLh=12;
			nP=4;
		}


		//if(m2d.getSliceAngle()<1.6 && this.getSliceAngle()>6) nP=4;


		if(flux3DNeeded){

			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++)

				{
					Mat R=util.rotMat2D(t*PI/2);
					for(int j=m2d.region[ir].getFirstEl();j<=m2d.region[ir].getLastEl();j++){

						Vect B=m2d.element[j].getB();
						if(t!=0) B=R.mul(B);
						B=B.v3();
						this.element[this.region[nr].getFirstEl()+p].setB(B);
						p++;
					}
				}
		}

		if(force3DNeeded){

			double spaceFaqcor=.96;

			p=0;

			int he=m2d.nElVert;


			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++){

					Mat R=util.rotMat2D(t*PI/2);			

					{
						boolean[] nc=new boolean[m2d.numberOfNodes+1];
						for(int j=m2d.region[ir].getFirstEl();j<=m2d.region[ir].getLastEl();j++){
							int[] vn=m2d.element[j].getVertNumb();

							int[] vn2=this.element[this.region[nr].getFirstEl()+p].getVertNumb();
							double kf=abs(this.node[vn2[vn2.length-1]].getCoord(2)-this.node[vn2[0]].getCoord(2))/2;

							//============
							//============
							kf*=spaceFaqcor;
							//============
							//============


							p++;

							for(int kk=0;kk<m2d.nElVert;kk++){
								if(nc[vn[kk]]) continue;

								double hx=this.node[vn2[kk]].getCoord(2);


								int ih=(int)Math.round(hx/.025*6);
								if(ih<4)
									ih=0;
								else
									ih=1;
								ih=0;

								Vect F=this.forceLamin[ih][this.mapnr[vn[kk]]];

								//	Vect F=m2d.node[vn[kk]].F;


								if(F!=null) 
								{

									if(t!=0 && this.coordCode==0)
									{
										F=R.mul(F);  
									}

									F=F.v3().times(kf);

									if(this.node[vn2[kk]].F!=null){

										this.node[vn2[kk]].setF(F.add(this.node[vn2[kk]].F));
									}
									else
										this.node[vn2[kk]].setF(F);

									if(this.node[vn2[kk+he]].F!=null)
										this.node[vn2[kk+he]].setF(F.add(this.node[vn2[kk+he]].F));
									else
										this.node[vn2[kk+he]].setF(F);



								}


								nc[vn[kk]]=true;


							}
						}
					}

				}


		}



	}

	public void transfer2DTo3DLayered(String file,  boolean flux3DNeeded,boolean force3DNeeded){

		// the following part enforces the periodic bc from 2d to the 3d model 
		for(int j=1;j<=m2d.numberOfNodes;j++){
			Vect v=m2d.node[j].getCoord();
			if(m2d.node[j].F!=null && v.el[0]<1e-5){

				m2d.node[j].F=new Vect(0,0);
			}
		}

		for(int i=1;i<=this.numberOfNodes;i++){
			this.node[i].F=null;
		}


		int ir=8;
		int p=0;
		int nr=1;
		int nL=6;
		int nP=1;

		int nLh=nL;
		/*	if(this.tag==17 ) {nL=6; nLh=nL;nP=4;}
	else
		if(this.tag==19 )*/ {nL=6; nLh=nL/2;}

		//=================
		if(this.numberOfElements<10000) 
		{nLh=3;
		nP=1;

		}
		else if(this.numberOfElements<20000) 
		{nLh=6;
		nP=1;
		}
		else if(this.numberOfElements<50000) 
		{nLh=6;
		nP=4;
		}

		else if(this.numberOfElements>60000){
			nLh=12;
			nP=4;
		}


		//if(m2d.getSliceAngle()<1.6 && this.getSliceAngle()>6) nP=4;


		if(flux3DNeeded){

			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++)

				{
					Mat R=util.rotMat2D(t*PI/2);
					for(int j=m2d.region[ir].getFirstEl();j<=m2d.region[ir].getLastEl();j++){

						Vect B=m2d.element[j].getB();
						if(t!=0) B=R.mul(B);
						B=B.v3();
						this.element[this.region[nr].getFirstEl()+p].setB(B);
						p++;
					}
				}
		}

		if(force3DNeeded){

			double spaceFaqcor=.96;

			p=0;

			int he=m2d.nElVert;


			for(int k=0;k<nLh;k++)	
				for(int t=0;t<nP;t++){

					Mat R=util.rotMat2D(t*PI/2);			

					{
						boolean[] nc=new boolean[m2d.numberOfNodes+1];
						for(int j=m2d.region[ir].getFirstEl();j<=m2d.region[ir].getLastEl();j++){
							int[] vn=m2d.element[j].getVertNumb();

							int[] vn2=this.element[this.region[nr].getFirstEl()+p].getVertNumb();
							double kf=abs(this.node[vn2[vn2.length-1]].getCoord(2)-this.node[vn2[0]].getCoord(2))/2;

							//============
							//============
							kf*=spaceFaqcor;
							//============
							//============


							p++;

							for(int kk=0;kk<m2d.nElVert;kk++){
								if(nc[vn[kk]]) continue;

								double hx=this.node[vn2[kk]].getCoord(2);

								int ih=(int)Math.floor((hx-1e-6)/.025*6);

								Vect F=this.forceLamin[ih][this.mapnr[vn[kk]]];


								if(F!=null) 
								{

									if(t!=0 && this.coordCode==0)
									{
										F=R.mul(F);  
									}

									F=F.v3().times(kf);

									if(this.node[vn2[kk]].F!=null){

										this.node[vn2[kk]].setF(F.add(this.node[vn2[kk]].F));
									}
									else
										this.node[vn2[kk]].setF(F);

									if(this.node[vn2[kk+he]].F!=null)
										this.node[vn2[kk+he]].setF(F.add(this.node[vn2[kk+he]].F));
									else
										this.node[vn2[kk+he]].setF(F);



								}


								nc[vn[kk]]=true;


							}
						}
					}

				}


		}



	}

}

