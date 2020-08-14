package io;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Scanner;

import javax.swing.JOptionPane;

import materialData.BHCurve;
import materialData.BHSCurve;
import materialData.CurrentWaveForm;
import math.Complex;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;
import fem.*;
import fem.Network.ElemType;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 15, 2012.
 */
public class Loader {

	private String regex="[:; ,=\\t]+";
	public String regex2="[\\[\\]\\s: )(,=\\t]+";


	

	public void loadMesh(Model model, String bunFilePath){

		model.meshFilePath=bunFilePath;
		
		try{
			FileReader fr=new FileReader(bunFilePath);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			String elType=br.readLine();
			model.setElType(elType);

			br.readLine();	
			
			line=br.readLine();
			sp=line.split(regex);
			if(!sp[0].equals(""))		
				model.numberOfNodes=Integer.parseInt(sp[0]);
			else
				model.numberOfNodes=Integer.parseInt(sp[1]);


			br.readLine();		
			line=br.readLine();

			sp=line.split(regex);
			if(!sp[0].equals(""))		
				model.numberOfElements=Integer.parseInt(sp[0]);
			else
				model.numberOfElements=Integer.parseInt(sp[1]);
		
			model.element=new Element[model.numberOfElements+1];
			for(int i=1;i<=model.numberOfElements;i++){
				model.element[i]=new Element(elType);
			}
			
		
			model.node=new Node[model.numberOfNodes+1];

			for(int i=1;i<=model.numberOfNodes;i++)
				model.node[i]=new Node(i, model.dim);

			br.readLine();		
			line=br.readLine();
			sp=line.split(regex);

			int nRegs=0;
			if(!sp[0].equals(""))		
				nRegs=Integer.parseInt(sp[0]);
			else
				nRegs=Integer.parseInt(sp[1]);
			
			if(model.numberOfRegions<nRegs) model.numberOfRegions=nRegs;

			model.region=new Region[model.numberOfRegions+1];
			for(int i=1;i<=model.numberOfRegions;i++)
				model.region[i]=new Region(model.dim);

			br.readLine();

			line=br.readLine();


			model.scaleFactor=Double.parseDouble(line);

			double factor=1.0/model.scaleFactor;

			for(int i=1;i<=model.numberOfElements;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;
				for(int j=0;j<sp.length;j++){
					if(!sp[j].equals(""))		
						model.element[i].setVertNumb(k++, Integer.parseInt(sp[j]));		
				}
			}

			
			Vect z=new Vect(model.dim);
			

			
			for(int i=1;i<=model.numberOfNodes;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;
				
				for(int j=0;j<sp.length;j++)
					if(!sp[j].equals(""))
						z.el[k++]=Double.parseDouble(sp[j])*factor;
				
	
				model.node[i].setCoord(z);

	
				}
			
			
			model.setBounds();
			
		
			
				for(int i=1;i<=nRegs;i++){
					
					line=br.readLine();
					if(line==null) line="1,0,x"+i;
					sp=line.split(regex);
					String[] str=new String[10];
					int k=0;
					for(int j=0;j<sp.length;j++){
						if(!sp[j].equals(""))
							str[k++]=sp[j];
						
					}
					model.region[i].setFirstEl(Integer.parseInt(str[0]));	
					
					model.region[i].setLastEl(Integer.parseInt(str[1]));
					model.region[i].setName(str[2]);
					model.region[i].setMaterial(str[2]);

			}
				
				
			for(int i=nRegs+1;i<=model.numberOfRegions;i++){
	
				model.region[i].setFirstEl(1);	
					
					model.region[i].setLastEl(0);
					model.region[i].setName("xxx");
					model.region[i].setMaterial("xmat");

			}
				
				
				line=br.readLine();
				if(line!=null)
				model.motor=getBooleanData(line);

				if(model.motor){
					line=br.readLine();
					if(line==null) return;
					sp=line.split(regex);
					String[] str=new String[3];
					int k=0;
					for(int j=0;j<sp.length;j++){
						if(!sp[j].equals(""))
							str[k++]=sp[j];
					}

					int rotBegin,rotEnd;

					rotBegin=Integer.parseInt(str[0]);
					rotEnd=Integer.parseInt(str[1]);

					for(int ir=rotBegin;ir<=rotEnd;ir++){
						model.region[ir].rotor=true;
					model.nRotReg++;
					}
				}
				
				//line=br.readLine();
				//if(line==null|| line.equals(""))
				//model.hasTwoNodeNumb=true;
				//else
				//model.hasTwoNodeNumb=getBooleanData(line);
				
				line=br.readLine();
				if(line==null|| line.equals(""))
				model.rotateConnect=false;
				else
				model.rotateConnect=getBooleanData(line);
		
				//==============
			
				if(model.motor) model.hasTwoNodeNumb=true;
				//=========
				
			System.out.println();
			System.out.println("Loading mesh file completed.");

			br.close();
			fr.close();

			for(int ir=1;ir<=nRegs;ir++)
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
					model.element[i].setRegion(ir);
			
			model.setMaxDim();
			
			model.setFemCalc();

						
				
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading model file.");
		}


	}	

	public void loadData(Model model,String dataFilePath){
		int dim0=model.dim;
	
		try{
			BufferedReader br = new BufferedReader(new FileReader(dataFilePath));
			String line;
			line=getNextDataLine(br);
			util.pr("// data type (0: Magnetic)");
			util.pr(line);
			int dataType =getIntData(line);
			model.dataType=dataType;
			
		
			util.pr("// dimension (2: 2D, 3:3D ");
			line=getNextDataLine(br);

			String[] sp=line.split(this.regex);	
		
			int ib=0;
			if(sp[0].equals("")) ib=1;
			int dim=Integer.parseInt(sp[ib++]);
		
			model.struc2D=0;
			if(dim==2 && ib<sp.length){
				model.struc2D=Integer.parseInt(sp[ib]);;	
			}

			line=getNextDataLine(br);
			util.pr("// coordCode (0: Cartesian, 1: Cylindrical ");
			int coordCode =getIntData(line);
			model.coordCode=coordCode;
			
			if(dim!=dim0){

				System.err.println("Mesh and Data do not match in dimension: "+dim0+" and "+dim);
			}
		


			if(dataType==0)
				setDataMag( model,br);
			else if(dataType==1)
				setDataMech( model,br);
			else if(dataType==3)
				setDataCLN( model,br);
		
				br.close();
		}

		catch(IOException e){
			e.printStackTrace();
		}

		System.out.println();
		System.out.println("Loading data file completed.");

	}
	
	public void setDataMag(Model model,BufferedReader br){
		String line;
		String s;
		int dim=model.dim;

		try {
	
			line=getNextDataLine(br);
			
			util.pr("// ANALYSIS MODE (0: Magnetostatic, 1:  A-method,  2: A-fi-method ");
			
			int am =getIntData(line);
			model.analysisMode=am;
			
			line=getNextDataLine(br);
			util.pr("// NONLINEAR (0: Linear , 1: Nonliear ");

			boolean nonlin=getBooleanData(line);;
			
			model.setNonLin(nonlin);

			line=getNextDataLine(br);
			model.AC=getBooleanData(line);
			
			line=getNextDataLine(br);;
			int nRegions =getIntData(line);		
			if(nRegions!=model.numberOfRegions){
				System.err.println("Mesh and Data do not match in the number of regions: "+model.numberOfRegions+" and "+nRegions);
			}
	

		
			
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				line=getNextDataLine(br);
			readAndSetRegMagPropery(model,ir,line);
			}
			
		
			model.BCtype=new int[model.nBoundary];
			for(int j=0;j<model.nBoundary;j++)
				model.BCtype[j]=-1;
			model.PBCpair=new int[model.nBoundary];

			Vect B;
			int[] bcData=new int[2];
			
			for(int j=0;j<model.nBoundary;j++){
				line=getNextDataLine(br);

				if(model.BCtype[j]>-1) continue;
				bcData=getBCdata(line);
				model.BCtype[j]=bcData[0];
				
				model.PBCpair[j]=bcData[1];
				if(model.BCtype[j]>1){
				
					model.BCtype[model.PBCpair[j]]=bcData[0];
					model.PBCpair[model.PBCpair[j]]=j;

				 model.hasPBC=true;

				}

								
			}
			
			line=getNextDataLine(br);

			int numbRegsWithJ=Integer.parseInt(line);
			for(int j=0;j<numbRegsWithJ;j++){
				line=getNextDataLine(br);
				String[] sp=line.split(this.regex);	

				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nr=Integer.parseInt(sp[ib++]);
				Vect J=new Vect(3);
				
				for(int k=0;k<3;k++)
					J.el[k]=Double.parseDouble(sp[ib++]);
				
				model.region[nr].setJ(J);

				
			}

			line=getNextDataLine(br);
			model.hasBunif=getBooleanData(line);
			if(model.hasBunif){
				line=getNextDataLine(br);
				double[] array=getCSV(line);
				
				model.unifB=new Vect(array);
			}
			
			line=getNextDataLine(br);
		
			model.eddyTimeIntegMode=getIntData(line);	
			
			line=getNextDataLine(br);
			model.dt=getScalarData(line);	

			line=getNextDataLine(br);
			model.rotSpeed=getScalarData(line);	
			
			line=getNextDataLine(br);
			model.meshAngStep=getScalarData(line);	

			line=getNextDataLine(br);
			if(line!=null){
				int nSteps=1;
				String sp[]=line.split(regex);
				int L=sp.length;
				
				int n1=0,n2=0, d=1;
		
					n1=Integer.parseInt(sp[L-3]);
					n2=Integer.parseInt(sp[L-2]);
					d=Integer.parseInt(sp[L-1]);
					if(d!=0)
					 nSteps=(n2-n1)/d+1;
			
			
			
					model.setnTsteps(nSteps);
					model.nBegin=n1;
					model.nEnd=n2;
					model.nInc=d;



				}
			
			
			line=getNextDataLine(br);
			
			int numBHdata=getIntData(line);	
			
			int jx=0;
			for(int j=0;j<numBHdata;j++){
				line=getNextDataLine(br);
				
				String[] sp=line.split(regex);
				int ib=0;
				if(sp[ib].equals("")) ib++;
				int bhID=Integer.parseInt(sp[ib]);
				for(int ir=1;ir<=model.numberOfRegions;ir++)
					if(model.region[ir].BHnumber==bhID)
						model.region[ir].setMaterial(sp[ib+1]);
				
				
			}
				
			
		


	for(int j=0;j<30+model.numberOfRegions;j++){
		
			line=br.readLine();
			if(line==null) continue;
			if(line.startsWith("MS")){
						
				String sp[]=line.split(regex);
				int L=sp.length;
								
					int nr=Integer.parseInt(sp[L-3]);
					model.region[nr].MS=true;

					Vect E=new Vect(Double.parseDouble(sp[L-2]),0,0);
					model.region[nr].setYng(E);
					Vect v=new Vect(Double.parseDouble(sp[L-1]),0,0);
					model.region[nr].setPois(v);	
					
					if(model.region[nr].isotElast){
					Vect sh=new Vect(1);
						v.el[0]=.5*E.el[0]/(1+v.el[0]);
						
						model.region[nr].setShear(sh);

					}
					

				}
			else if(line.startsWith("loadFlux")) 
				model.loadFlux=this.getBooleanData(line);
			else if(line.startsWith("fluxFolder")){

			line=br.readLine();
				model.fluxFolderIn=line;
			}
			else if(line.startsWith("saveFlux")) model.saveFlux=this.getBooleanData(line);
			else if(line.startsWith("saveForce")) model.saveForce=this.getBooleanData(line);
			else if(line.startsWith("mag")) model.magAnalysis=this.getBooleanData(line);
			else if(line.startsWith("trans")) model.transfer2DTo3D=this.getBooleanData(line);	
			else if(line.startsWith("forceCal")) model.forceCalcMode=this.getIntData(line);	
			else if(line.startsWith("rotate")) model.rotateRotor=this.getBooleanData(line);	
			else if(line.startsWith("loadPrev")) model.loadPrevMag=this.getBooleanData(line);
			else if(line.startsWith("axi")) model.axiSym=this.getBooleanData(line);	
			else if(line.startsWith("height")) model.height=this.getScalarData(line);
			else if(line.startsWith("POD")) model.POD=this.getIntData(line);
			else if(line.startsWith("snapShot")) model.snapShot=this.getIntData(line);
		}

	//util.pr(model.fluxFolderIn);

	if(model.axiSym) model.height=2*Math.PI;
		
			model.setFreq(100);
			model.setHasJ();

			model.setHasM();

			model.setHasMS();


			model.setNonLinToElememts();
			
			model.setEdge();


			model.setElementsParam();

			
			model.setBounds();

		
		if(model.coordCode==1) {
			
			model.cpb=1;

			for(int j=0;j<model.nBoundary;j++){

				if(model.BCtype[j]==3) model.cpb=-1;
			}
			if(model.hasTwoNodeNumb && model.rotateConnect && model.nRotReg>0){

			 model.mapCommonNodes();	
			}
		}
		
		//=====================


		model.setNodeOnBound();
		
		if(model.hasPBC) model.mapPBC();

		
			if(model.deform)
				model.setMechBC();
			
			
			
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].circuit ) {
					model.circuit=true;
					break;
				}
				
			}

		
			int ix=0;
			int iy=0;

			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].circuit ) {
					model.region[ir].currentIndex=ix;
			ix++;
			}

			}
			
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].circuit && model.region[ir].curMap1==0)
				{
				model.region[ir].unknownCurrentIndex=iy;
			iy++;
			}

			}
			
			
			model.numberOfCurrents=ix;
			
			model.numberOfUnknownCurrents=iy;
			
			
			if(ix>0){/*
				model.analysisMode=1;
				model.eddyTimeIntegMode=-2;
			*/}
		//--------------
			
			
			model.threePhaseRegs=new int[3];
			if(model.circuit){
	

				if(model.numberOfRegions>8){
					model.threePhaseRegs[0]=9;
					model.threePhaseRegs[1]=11;
					model.threePhaseRegs[2]=13;
				}
					
					 if(model.numberOfRegions==8){
						model.threePhaseRegs[0]=1;
						model.threePhaseRegs[1]=3;
						model.threePhaseRegs[2]=5;
					}
					 if(model.numberOfRegions==4){
							model.threePhaseRegs[0]=1;
							model.threePhaseRegs[1]=0;
							model.threePhaseRegs[2]=0;
						}

				}
			
		//)))))))))))))))))))))))
	
			if(model.threePhaseRegs[0]!=0 &&model.threePhaseRegs[1]!=0&& model.threePhaseRegs[2]!=0)
				model.nNeutral=1;
			
			//model.nNeutral=0;

		
			int nCur=model.numberOfUnknownCurrents;
			
			model.unCurRegNumb=new int[nCur];
		
			for(int ir=1;ir<=model.numberOfRegions;ir++)
				if(model.region[ir].circuit && model.region[ir].curMap1==0)
					model.unCurRegNumb[model.region[ir].unknownCurrentIndex]=ir;
			


			try {
				model.ia=new CurrentWaveForm("emf//Ian.txt");
				model.ib=new CurrentWaveForm("emf//Ibn.txt");
				model.ic=new CurrentWaveForm("emf//Icn.txt");

				/*	model.va=new CurrentWaveForm("emf//Va.txt");
				model.vb=new CurrentWaveForm("emf//Vb.txt");
				model.vc=new CurrentWaveForm("emf//Vc.txt");*/
			} catch (Exception e) {
				e.printStackTrace();
			}
		
		
		
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		//model.setForceCalc();

}

	public void setDataCLN(Model model,BufferedReader br){

		String line;
		String s;
		try {

			model.analysisMode=0;
			util.pr(" ---0000000000-----> "+model.dim);
			line=getNextDataLine(br);
			util.pr(" ------------> "+line);
			int nRegions =getIntData(line);		

			if(nRegions!=model.numberOfRegions){
				System.err.println("Mesh and Data do not match in the number of regions: "+model.numberOfRegions+" and "+nRegions);
			}
	
			line=getNextDataLine(br);;
			int nStage =getIntData(line);		
			
			model.nCLNstages=nStage;
			
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				line=getNextDataLine(br);
			readAndSetRegMagPropery(model,ir,line);
			}

			model.setElementsParam();

			
			model.setBounds();
			
			
			
			model.BCtype=new int[model.nBoundary];
			for(int j=0;j<model.nBoundary;j++)
				model.BCtype[j]=-1;
			model.PBCpair=new int[model.nBoundary];

			Vect B;
			int[] bcData=new int[2];
			
			for(int j=0;j<model.nBoundary;j++){
				line=getNextDataLine(br);

				if(model.BCtype[j]>-1) continue;
				bcData=getBCdata(line);
				model.BCtype[j]=bcData[0];
				
				model.PBCpair[j]=bcData[1];
				if(model.BCtype[j]>1){
				
					model.BCtype[model.PBCpair[j]]=bcData[0];
					model.PBCpair[model.PBCpair[j]]=j;

				 model.hasPBC=true;

				}

								
			}
			
			line=getNextDataLine(br);
		
			int numbRegsWithJ=Integer.parseInt(line);
		
			model.phiCoils=new PhiCoil[numbRegsWithJ];
			
			
			for(int j=0;j<numbRegsWithJ;j++){
				line=getNextDataLine(br);
				String[] sp=line.split(this.regex);	

				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nr=Integer.parseInt(sp[ib++]);
				model.phiCoils[j]=new PhiCoil(nr);
				model.phiCoils[j].index=j;
				double regSigma=0;
				if(model.region[nr].getSigma()!=null )
					regSigma=model.region[nr].getSigma().el[0];
						
						
				if(regSigma>0) model.phiCoils[j].setSigma(regSigma);
				
				else{
				if(ib<=sp.length){
				double sigma=Double.parseDouble(sp[ib++]);

				model.phiCoils[j].setSigma(sigma);
				}
				}
		
				
				
				double[][] boxdata=new double[2][6];
				
				for(int k=0;k<2;k++){
				
				line=getNextDataLine(br);
				
				ib=0;
				if(sp[0].equals("")) ib=1;
				sp=line.split(this.regex);	

				String ss=sp[ib++];
				if(ss.equals("x") || ss.equals("r"))
				{
					boxdata[k][0]=Double.parseDouble(sp[ib++]);
					boxdata[k][1]=Double.parseDouble(sp[ib++]);
					
				}
				ss=sp[ib++];
				if(ss.equals("y") || ss.equals("t"))
				{
					boxdata[k][2]=Double.parseDouble(sp[ib++]);
					boxdata[k][3]=Double.parseDouble(sp[ib++]);
					
				}
				ss=sp[ib++];
				if(ss.equals("z"))
				{
					boxdata[k][4]=Double.parseDouble(sp[ib++]);
					boxdata[k][5]=Double.parseDouble(sp[ib++]);
				}
				model.phiCoils[j].faceBox[k]=boxdata[k];
				}
				
			}

			line=getNextDataLine(br);
		
			model.hasBunif=getBooleanData(line);
			if(model.hasBunif){
				line=getNextDataLine(br);
				double[] array=getCSV(line);
				
				model.unifB=new Vect(array);
			}
		
		line=getNextDataLine(br);
		line=util.dropLeadingSpaces(line);
		if(line.equals("NETWORK")){
			Network network=new Network();
			network.read(this, br);
			
			model.network=network;
			
			for(int j=0;j<network.numElements;j++){
				if(network.elems[j].type==ElemType.FEM){
					for(int k=0;k<model.phiCoils.length;k++){
						if(model.phiCoils[k].regNo==network.elems[j].fem_id)
							network.elems[j].fem_index=k;
					}
				}
			}
		}
		line=getNextDataLine(br);

		if(line==null|| line.equals("")) 
			model.open_vps=true;
		else
		model.open_vps=getBooleanData(line);
		


	if(model.axiSym) model.height=2*Math.PI;


			model.setHasMS();

			
			model.setEdge();

		
		if(model.coordCode==1) {
			
			model.cpb=1;

			for(int j=0;j<model.nBoundary;j++){

				if(model.BCtype[j]==3) model.cpb=-1;
			}
			if(model.hasTwoNodeNumb && model.rotateConnect && model.nRotReg>0)
			 model.mapCommonNodes();	
		}
		
		//=====================


		model.setNodeOnBound();
		
		if(model.hasPBC) model.mapPBC();

			
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].circuit ) {
					model.circuit=true;
					break;
				}
				
			}
		
		
		
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


}

	public void setDataMech(Model model,BufferedReader br){

		String line;

		try {


			line=getNextDataLine(br);
		
			model.motor=getBooleanData(line);
	
			line=getNextDataLine(br);
			int nRegions =getIntData(line);		
			if(nRegions!=model.numberOfRegions){
				System.err.println("Mesh and Data do not match in  number of regions: "+model.numberOfRegions+" and "+nRegions);
			}

			
		
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				line=getNextDataLine(br);
			readAndSetRegDataMech(model,ir,line);
			}
			
					
			model.BCtype=new int[model.nBoundary];
			for(int j=0;j<model.nBoundary;j++)
				model.BCtype[j]=-1;
			model.PBCpair=new int[model.nBoundary];


			int[] bcData=new int[2];
			
			for(int j=0;j<model.nBoundary;j++){

				line=getNextDataLine(br);

				if(model.BCtype[j]>-1) continue;
				bcData=getBCdataMech(line);
				model.BCtype[j]=bcData[0];
				
				if(model.BCtype[j]>1){
					model.PBCpair[j]=bcData[1];
					model.BCtype[model.PBCpair[j]]=bcData[0];
					model.PBCpair[model.PBCpair[j]]=j;


				}
		
				
				if(model.BCtype[j]>1) model.hasPBC=true;

			}
			
			line=getNextDataLine(br);

			int numBC=getIntData(line);		
			
			model.mechBC=new String[numBC][50];

			for(int j=0;j<numBC;j++){
				line=getNextDataLine(br);
			String sp[]=line.split(regex);
			for(int k=0;k<sp.length;k++){
				model.mechBC[j][k]=sp[k];
			}
			}
			
			String[] forceFile=null;
			
			for(int j=0;j<5;j++){
				
			
				line=br.readLine();

				if( line.startsWith("force")) {
					
					line=br.readLine();
					model.forceFolder=line;

					break;
				}
			}
		
			for(int j=0;j<5;j++){

				if( line.startsWith("force")) {
			
					String sp[]=line.split(regex);
			
					int L=sp.length;
					model.nBegin=Integer.parseInt(sp[L-3]);
					model.nEnd=Integer.parseInt(sp[L-2]);
					model.nInc=Integer.parseInt(sp[L-1]);
					
					int N=0;
					for(int i=model.nBegin;i<=model.nEnd;i+=	model.nInc)
						N++;
						
					forceFile =new String[N];
					int ix=0;
					for(int i=model.nBegin;i<=model.nEnd;i+=	model.nInc){
						
						int nf=i%1800;
						if(model.forceFolder.startsWith("cent")){
							model.centrigForce=true;
							 sp=model.forceFolder.split(regex);	
							 model.rpm=Double.parseDouble(sp[1]);
						forceFile[ix]=null;
						}
						else
						forceFile[ix]=model.forceFolder+"\\force"+nf+".txt";
					
						
						ix++;
					}
					

					break;
				}
				line=br.readLine();

				}
			
		
			if(forceFile!=null){
			
			model.forceFile=forceFile;
			model.nTsteps=forceFile.length;
			}

			
			
			if(model.nTsteps>0){
			line=br.readLine();

			line=br.readLine();
			
				model.timeIntegMode=getIntData(line);
			
			line=br.readLine();

			line=br.readLine();

				model.setDt(model.nInc*getScalarData(line));
			}
			else
				model.timeIntegMode=0;
			
		//	line=br.readLine();
		//	line=br.readLine();
			line=getNextDataLine(br," /* CONTACT");
		//	util.pr(line);
			 if(line.startsWith("contact")) {
				model.contact=new ContactAnalysis();
			model.contact.readContact(this, br,model);

			}

	for(int j=0;j<30;j++){
		
		line=br.readLine();
	
		if(line==null) continue;
	
		if(line.startsWith("loadFlux")) model.loadFlux=this.getBooleanData(line);
		else if(line.startsWith("loadForce")) model.loadForce=this.getBooleanData(line);
		else if(line.startsWith("loadDisp")) model.loadDisp=this.getBooleanData(line);
		else if(line.startsWith("saveFlux")) model.saveFlux=this.getBooleanData(line);
		else if(line.startsWith("saveForce")) model.saveForce=this.getBooleanData(line);
		else if(line.startsWith("saveDisp")) model.saveDisp=this.getBooleanData(line);
		else if(line.startsWith("saveStress")) model.saveStress=this.getBooleanData(line);
		else if(line.startsWith("mag")) model.magAnalysis=this.getBooleanData(line);
		else if(line.startsWith("mec")) model.mechAnalysis=this.getBooleanData(line);
		else if(line.startsWith("trans")) model.transfer2DTo3D=this.getBooleanData(line);
		else if(line.startsWith("loadPrev")) model.loadPrevMech=this.getBooleanData(line);

		else if(line.startsWith("modal")) model.modal=this.getBooleanData(line);
		else if(line.startsWith("rayA")) model.rayAlpha=this.getScalarData(line);
		else if(line.startsWith("rayB")) model.rayBeta=this.getScalarData(line);



	}
	
	//model.deform=!model.loadDisp;
	model.deform=true;
			
			model.setHasThermal();
			

			model.setElementsParamMech();

		
		if(model.coordCode==1) {
			
			model.cpb=1;

			for(int j=0;j<model.nBoundary;j++){

				if(model.BCtype[j]==3) model.cpb=-1;
			}
			model.setSliceBounds();
			if(model.hasTwoNodeNumb && model.rotateConnect && model.nRotReg>0)
			 model.mapCommonNodes();	
		}
		
		
	
		//=====================


		model.setNodeOnBound();
		
		if(model.hasPBC) model.mapPBC();


			if(model.deform)
				model.setMechBC();
	
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		


}
	
	
		
	public void loadBH(BHCurve BH,String mateName) throws Exception{
		double[][] BH1=new double[200][2];

		String file = System.getProperty("user.dir") + "\\BH\\"+mateName+".txt";

		try{
			Scanner scr=new Scanner(new FileReader(file));
			while(scr.hasNext()){

				while(!scr.next().equals("begin")){}
				int j=0;
				String s=scr.next();
				while(!s.equals("end")){
					BH1[j][0]=Double.parseDouble(s);
					s=scr.next();
					BH1[j++][1]=Double.parseDouble(s);
					s=scr.next();
				}
				BH.length=j;
			}

			BH.BH=new double[BH.length][2];
			for(int k=0;k<BH.length;k++)
				BH.BH[k]=BH1[k];

			scr.close();
		}
		catch(IOException fnf){
			throw new Exception(fnf);
		}
		BH.setGradBH();

	}

	public void loadBHS(BHSCurve BHS,String mateName, boolean vms) throws Exception{
		double[][] BH1=new double[2000][2];
		BHSCurve BHS1=new BHSCurve(100);


		String file = System.getProperty("user.dir") + "\\BH\\"+mateName+".txt";

		try{
			BufferedReader br = new BufferedReader(new FileReader(file));

			String line=null;

			int i=0;

			int nn=0;
			
			for( i=0;i<BHS1.nstr;i++){

				BH1= new double[200][2];
				nn=0;
				line="";
				while( nn<10 &&  !line.startsWith("str")){ 
					line=br.readLine();
					if(line==null) { line=""; nn++;}
				
					}
				
				if(nn>=10) break;

				BHS1.stress[i]=getScalarData(line);
				line="";
				while(!line.startsWith("begin")){ line=br.readLine();}
				int j=0;
				line=br.readLine();
				while(!line.startsWith("end")){ 

					if(line.length()!=0 )
						{
					BH1[j]=getPair(line);
					j++;
					line=br.readLine();
						}
					else
						line=br.readLine();



				}


				double[][] BH2= new double[j][2];

				for(int k=0;k<j;k++){
					BH2[k]=BH1[k];

				}
				
				BHS1.BH[i]=new BHCurve(BH2);
			}
	
		
			if(vms){
				int nz=-1;
				for(int p=0;p<BHS1.nstr;p++)
					if(abs(BHS1.stress[p])==0) {nz=p; break;}
				if(nz==-1) throw new IllegalArgumentException(file+"contains NO BH data for zero stress");
				BHS.nstr=1;
				BHS.stress=new double[BHS.nstr];
				BHS.BH=new BHCurve[BHS.nstr];
				BHS.BH[0]=new BHCurve(BHS1.BH[nz].BH);
				BHS.nZeroStress=0;
				return;
			}
			
			
			
				
			BHS.nstr=i;
			BHS.BH=new BHCurve[BHS.nstr];
			BHS.stress=new double[BHS.nstr];
			
			for(int k=0;k<i;k++)
				BHS.stress[k]=BHS1.stress[k];
			
			Vect ss=new Vect(BHS.stress);
			int[] indx=ss.bubble();
			
			for(int k=0;k<i;k++)
				BHS.BH[k]=new BHCurve(BHS1.BH[indx[k]].BH);
		
	
			for(int p=0;p<BHS.nstr;p++)
				if(abs(BHS.stress[p])<1.0) {BHS.nZeroStress=p; break;}


			
			br.close();
		}

		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading BHS file.");
		}


	}


	

	public boolean loadFlux(Model model,String fluxFilePath){
			return loadFlux(model,fluxFilePath,0);
	}
	
public boolean loadFlux(Model model,String fluxFilePath, double angDeg){

	boolean rotating=true;
		if(angDeg==0) rotating=false;
		
		Mat R=new Mat();
		if(rotating)
			R=util.rotMat2D(angDeg*Math.PI/180);
		
		try{
			Scanner scr=new Scanner(new FileReader(fluxFilePath));

			scr.next();
			int dim=Integer.parseInt(scr.next());

			int nElements=Integer.parseInt(scr.next());
			

			if(nElements!=model.numberOfElements) {

			
				String msg="Flux file does not match the mesh";
				System.err.println(msg);
				return false;
				
			}
	
			
			double Bn2,Bmax2=-1e40,Bmin2=1e40;
			
			for(int ir=1;ir<=model.numberOfRegions;ir++)
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				
				Vect B1=new Vect(dim);
			
				for(int j=0;j<dim;j++)
					B1.el[j]=Double.parseDouble(scr.next());
				
				
				if(rotating && model.element[i].rotor) B1=R.mul(B1);
				
				//B1=model.getElementCenter(i).normalized();

				model.element[i].setB(B1);
				
				Bn2=B1.dot(B1);
				
				if(Bn2>Bmax2)
					Bmax2=Bn2;
				if(Bn2<Bmin2)
					Bmin2=Bn2;

			}


			model.Bmax=sqrt(Bmax2);
			model.Bmin=sqrt(Bmin2);
			scr.close();

			System.out.println("Flux was loaded from "+fluxFilePath+" to the model.");
			model.fluxLoaded=true;
			return true;
		}
		catch(IOException e){System.err.println("Error in loading flux file.");
		return false;
		}

	}	

public void average(String bun1, String bun2,String bun3){
	

	try{

		FileReader fr=new FileReader(bun1);
		BufferedReader br = new BufferedReader(fr);
		String line=br.readLine();
		 line=br.readLine();
			int dim1=Integer.parseInt(line);
		 line=br.readLine();
		int ne=Integer.parseInt(line);
		
		Vect[] flx1=new Vect[ne];
		for(int i=0;i<ne;i++){
			line=br.readLine();
			flx1[i]=new Vect(getTabedData(line));
		}
  

		 fr=new FileReader(bun2);
		 br = new BufferedReader(fr);
		 line=br.readLine();
		 line=br.readLine();
			int dim2=Integer.parseInt(line);
		 line=br.readLine();
		int ne2=Integer.parseInt(line);
		
		Vect[] flx2=new Vect[ne2];
		for(int i=0;i<ne2;i++){
			line=br.readLine();
			flx2[i]=new Vect(getTabedData(line));
		}                    

		for(int i=0;i<ne2;i++){
			flx2[i]=flx2[i].add(flx1[i]).times(.5);
		}   
		
		fr.close();
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(bun3)));	
		pw.println("flux");
		pw.println(dim1);
		pw.println(ne);

		for(int i=0;i<ne2;i++){
			for(int k=0;k<dim1;k++)
					pw.format("%E\t",flx2[i].el[k]);
			pw.println();
		} 
		
		pw.close();

	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading flux file.");
	}
	
	System.out.println("The output was written to "+bun3);

}	

	public Vect[] fluxSub(String file1,String file2){
	
		String f3 = System.getProperty("user.dir") + "//fluxdiff.txt";

		try{

			FileReader fr=new FileReader(new File(file1));
			BufferedReader br = new BufferedReader(fr);
			String line=br.readLine();
			 line=br.readLine();
				int dim1=Integer.parseInt(line);
			 line=br.readLine();
			int ne=Integer.parseInt(line);
			
			Vect[] flx1=new Vect[ne];
			for(int i=0;i<ne;i++){
				line=br.readLine();
				flx1[i]=new Vect(getTabedData(line));
			}
	  

			 fr=new FileReader(new File(file2));
			 br = new BufferedReader(fr);
			 line=br.readLine();
			 line=br.readLine();
				int dim2=Integer.parseInt(line);
			 line=br.readLine();
			int ne2=Integer.parseInt(line);
			
			Vect[] flx2=new Vect[ne2];
			for(int i=0;i<ne2;i++){
				line=br.readLine();
				flx2[i]=new Vect(getTabedData(line));
			}                    

			for(int i=0;i<ne2;i++){
				flx2[i]=flx2[i].sub(flx1[i]);
			}   
			
			fr.close();
			
			return flx2;
			/*PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f3)));	
			pw.println("flux");
			pw.println(dim1);
			pw.println(ne);

			for(int i=0;i<ne2;i++){
				for(int k=0;k<dim1;k++)
						pw.format("%E\t",flx2[i].el[k]);
				pw.println();
			} 
			
			pw.close();*/

		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading flux file.");
		}
		
		//System.out.println("The output was written to "+f3);
 return null;
	}	
	
	
	public boolean loadPotential(Model model,String vPotFile){

		try{
			Scanner scr=new Scanner(new FileReader(vPotFile));

			scr.nextLine();
			int dim=Integer.parseInt(scr.next());

			int nEdges=Integer.parseInt(scr.next());
			if(nEdges!=model.numberOfEdges) {
				String msg="Vector potential file does not match the mesh";
				JOptionPane.showMessageDialog(null, msg," ", JOptionPane. INFORMATION_MESSAGE);
				return false;
			}

			for(int i=1;i<=nEdges;i++){

				double A=Double.parseDouble(scr.next());

				model.edge[i].setA(A);

			}


			scr.close();

			System.out.println("Vector potential was loaded from "+vPotFile+" to the model.");
			model.potentialLoaded=true;
			return true;
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading vector potential file.");
		return false;
		}

	}	
	
	public void loadPrevMag(Model model,String vPotFile){

		try{
			Scanner scr=new Scanner(new FileReader(vPotFile));
	

			for(int i=0;i<model.numberOfUnknownEdges;i++){

				double A=Double.parseDouble(scr.next());

				int ix=model.unknownEdgeNumber[i+1];

				model.edge[ix].Ap=A;
				
				 A=Double.parseDouble(scr.next());
				model.edge[ix].setA(A);


			}

			for(int i=0;i<model.numberOfUnknownCurrents;i++){

				int nr=model.unCurRegNumb[i];
				double cur=Double.parseDouble(scr.next());

				model.region[nr].currentp=cur;
				 cur=Double.parseDouble(scr.next());
					model.region[nr].current=cur;
			}
	
			
			model.vNeutral=Double.parseDouble(scr.next());
			model.tet=Double.parseDouble(scr.next());
			model.tetp=Double.parseDouble(scr.next());

			model.tetpp=Double.parseDouble(scr.next());


			scr.close();

			System.out.println("Vector potential was loaded from "+vPotFile+" to the model.");
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading vector potential file.");
		}

	}
	
	public void loadPrevMech(Model model,String uFile){

		try{
			Scanner scr=new Scanner(new FileReader(uFile));
	

			model.up=new Vect(model.numberOfUnknownUcomp);
			model.upp=new Vect(model.numberOfUnknownUcomp);
			model.ud=new Vect(model.numberOfUnknownUcomp);
			model.udd=new Vect(model.numberOfUnknownUcomp);
				
			for(int i=0;i<model.numberOfUnknownUcomp;i++){

				model.up.el[i]=Double.parseDouble(scr.next());
				model.upp.el[i]=Double.parseDouble(scr.next());
				model.ud.el[i]=Double.parseDouble(scr.next());
				model.udd.el[i]=Double.parseDouble(scr.next());
			
			}

		
			scr.close();

			System.out.println("initial displacement potential was loaded from "+uFile+" to the model.");
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading vector potential file.");
		}

	}

	public boolean loadStress(Model model,String stressFilePath){

		try{
			Scanner scr=new Scanner(new FileReader(stressFilePath));

			scr.next();
			int dim=Integer.parseInt(scr.next());
			int L=3*(dim-1);

			int nElements=Integer.parseInt(scr.next());
			if(nElements!=model.numberOfElements) {
				String msg="Stress file doesnt match the mesh";
				JOptionPane.showMessageDialog(null, msg," ", JOptionPane. INFORMATION_MESSAGE);
				return false;
			}
			scr.next();
			scr.next();
			scr.next();
			scr.next();

			while(scr.hasNext()){
				
				int ne=Integer.parseInt(scr.next());

					Vect ss=new Vect(L);
					for(int j=0;j<L;j++)
						ss.el[j]=Double.parseDouble(scr.next());
					model.element[ne].setDeformable(true);
					model.element[ne].setStress(ss);
			}
	


			scr.close();

			System.out.println("Stress was loaded from "+stressFilePath+" to the model.");
			return true;
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading stress file.");
		return false;
		}



	}

	public boolean loadNodalField(Model model,String nodalFilePath,int mode){
		
		if(nodalFilePath==null) return false;
		return loadNodalField(model,nodalFilePath,mode,1e10);
	}

	public boolean loadNodalField(Model model,String nodalFilePath,int mode,double angDeg){
		
		boolean rotating =true;
		if(angDeg==1e10) rotating=false;
		
		Mat R=new Mat();
		if(rotating)
			R=util.rotMat2D(angDeg*Math.PI/180);
	
		String str;

		try{
	

			Scanner scr=new Scanner(new FileReader(nodalFilePath));
			str=scr.next();
			if(!str.startsWith("nodal")) {
				String msg="File does not contain nodal field.";
				System.err.println(msg);
				//JOptionPane.showMessageDialog(null, msg," ", JOptionPane. INFORMATION_MESSAGE);
				return false;
			}

			str=scr.next();
			int dim=Integer.parseInt(str);
			//util.pr(str);
			str=scr.next();
			int nNodes=Integer.parseInt(str);
			//util.pr(str);
			if(nNodes!=model.numberOfNodes || dim!=model.dim) {
				String msg="Nodal field file does not match the mesh";
				System.err.println(msg);
				//JOptionPane.showMessageDialog(null, msg," ", JOptionPane. INFORMATION_MESSAGE);
				return false;
			}
			int nn;
			double sn2=0,smax2=0,smin2=0;
		
			while(scr.hasNext()){
				str=scr.next();
				//util.pr(str);
				nn=Integer.parseInt(str);
		
					Vect v=new Vect(dim);
					for(int j=0;j<dim;j++){
						str=scr.next();
					//	util.pr(str);
						v.el[j]=Double.parseDouble(str);		
			}
			
				
					if(model.coordCode==1){

						Mat R2=new Mat();

						double ang=-util.getAng(model.node[nn].getCoord().v2());
						
						R2=util.rotMat2D(ang);
						
						if(model.dim==2){
							v=R2.mul(v);		

						}
						else {
													
							Mat R3=new Mat(dim,dim);
							for(int p=0;p<2;p++)
								for(int q=0;q<2;q++)
									R3.el[p][q]=R2.el[p][q];

							R3.el[2][2]=1;		
							v=R3.mul(v);	
						

						
						}

					
					}
						
						if(mode==1){
						
						if(model.timeIntegMode==4)
							if(model.node[nn].F!=null){
								model.node[nn].Fstr=model.node[nn].F.deepCopy();
							}
					///v.hshow();	
							model.node[nn].setF(v);
					
						}
						else if(mode==2)
							model.node[nn].setFms(v);
						else if(mode==3)
							model.node[nn].setF(v);
						else if(mode==4){
							model.node[nn].Fstr=v.deepCopy();
						}
						else if(mode==-1){
							v.timesVoid(1e-9);
						model.node[nn].setU(v);
						}
					
						sn2=v.dot(v);
			
						
					if(sn2>smax2)
						smax2=sn2;
					if(sn2<smin2)
						smin2=sn2;
			}

		

			if(mode==-1){
				model.uMax=sqrt(smax2);

			}
			else	if(mode==1){
				model.FreluctMax=sqrt(smax2);
			}

			else	if(mode==2){
				model.FmsMax=sqrt(smax2);

			}
			if(mode>0)
				model.forceLoaded=true;


			scr.close();


			System.out.println("Force was loaded from "+nodalFilePath+" to the model.");
			return true;
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading force.");


		return false;
		}
	}
	

	

	public boolean loadNodalFieldTemporal(Model model,String nodalFilePath,int mode){
	
	

		try{
			Scanner scr=new Scanner(new FileReader(nodalFilePath));

			scr.next();
			int dim=Integer.parseInt(scr.next());
			model.dim=dim;
			int nNodes=Integer.parseInt(scr.next());

			if(nNodes!=model.numberOfNodes) {
				String msg="Nodal field file does not match the mesh";
				JOptionPane.showMessageDialog(null, msg," ", JOptionPane. INFORMATION_MESSAGE);
				return false;
			}
			int nn;
			double sn2=0,smax2=0,smin2=0;

			while(scr.hasNext()){
	
				nn=Integer.parseInt(scr.next());

			
					Vect v=new Vect(dim);
					for(int j=0;j<dim;j++)
						v.el[j]=Double.parseDouble(scr.next());			
			
					
					if(model.coordCode==1){

						Mat R2=new Mat();

						double ang=-util.getAng(model.node[nn].getCoord().v2());
						
						R2=util.rotMat2D(ang);
						
						if(model.dim==2){
							v=R2.mul(v);		

						}
						else {
													
							Mat R3=new Mat(dim,dim);
							for(int p=0;p<2;p++)
								for(int q=0;q<2;q++)
									R3.el[p][q]=R2.el[p][q];

							R3.el[2][2]=1;		
							v=R3.mul(v);	
						

						
						}

					
					}
						
						if(mode==1){
							model.node[nn].setF(v);

						}
						else if(mode==2)
							model.node[nn].setFms(v);
						else if(mode==3)
							model.node[nn].setF(v);
						else if(mode==-1){
							v.timesVoid(1e-9);
						model.node[nn].setU(v);
						}
						sn2=v.dot(v);
			
						
					if(sn2>smax2)
						smax2=sn2;
					if(sn2<smin2)
						smin2=sn2;
			}

		

			if(mode==-1){
				model.uMax=sqrt(smax2);

			}
			else	if(mode==1){
				model.FreluctMax=sqrt(smax2);
			}

			else	if(mode==2){
				model.FmsMax=sqrt(smax2);

			}
			if(mode>0)
				model.forceLoaded=true;


			scr.close();


			System.out.println("Force was loaded from "+nodalFilePath+" to the model.");
			return true;
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading force.");


		return false;
		}
	}
	


	private int[] getBCdata(String line){
	
		
		int[] bctp=new int[2];
		String[] sp=line.split(regex);	
		int pair=-1,bct=0;
	
		int ib=0;
		if(sp[0].equals("")) ib++;
		if(line.startsWith("N")){
				bct=0;
				
		}
		else if(line.startsWith("D")){
			bct=1;
			
	}
		else if(line.startsWith("PS")){
			bct=2;
			pair=Integer.parseInt(sp[ib+1])-1;
			
	}

		else if(line.startsWith("PA")){

				bct=3;
				pair=Integer.parseInt(sp[ib+1])-1;
			}
		else{

			bct=0;
		}
		
		bctp[0]=bct;
		bctp[1]=pair;
		
		return bctp;
	}


	private int[] getBCdata2D(String line){
		
		int[] bctp=new int[2];
		String[] sp=line.split("");	
		int pair=-1,bct=0;
		int k=0;
		while(!sp[k].equals(":")){ k++;}	
		k++;
		if(k==sp.length) return bctp;
		if(sp[k].equals(" ")) k++;
		String s=sp[k];

		if(s.equals("P")){

			if(sp[k+1].equals("S")){
				bct=2;
			}

			else if(sp[k+1].equals("A")){

				bct=3;
			}

			sp=line.split(regex);
			String s2=sp[sp.length-1];
			if(s2.equals(","))
				s2=sp[sp.length-2];
			pair=Byte.parseByte(s2)-1;

		}

		else if(s.equals("N"))
			bct=0;

		else if(s.equals("D")){
			sp=line.split(regex);
			bct=1;
			k=0;
			
		}
			
		bctp[0]=bct;
		bctp[1]=pair;

		return bctp;
	}
	
	private int[] getBCdataMech(String line){
	
		int[] bctp=new int[2];
		String[] sp=line.split("");	
		int pair=-1,bct=0;
		int k=0;
		while(!sp[k].equals(":")){ k++;}	
		k++;
		if(k==sp.length) return bctp;
		if(sp[k].equals(" ")) k++;
		String s=sp[k];

		if(s.equals("P")){

			if(sp[k+1].equals("S")){
				bct=2;
			}

			else if(sp[k+1].equals("A")){

				bct=3;
			}

			sp=line.split(regex);
			String s2=sp[sp.length-1];
			if(s2.equals(","))
				s2=sp[sp.length-2];
			pair=Byte.parseByte(s2)-1;

		}

		else if(s.equals("N"))
			bct=0;

		else if(s.equals("D")){
			sp=line.split(regex);
			bct=1;
			k=0;
			
		}
			
		bctp[0]=bct;
		bctp[1]=pair;
		
		return bctp;
	}

	
	public double[] getTabedData(String line){
		String[] sp=line.split(regex);	
		int L=sp.length;
		double[] v=new double[L];
		for( int p=0;p<L;p++)
			v[p]=Double.parseDouble(sp[p]);

		return v;
	}
	
	private double[] getPair(String line){
		String[] sp=line.split(regex);	
		double[] v=new double[2];
		int k=0;
		if(sp[k].equals("")) k++;

		for( int p=0;p<2;p++)
			v[p]=Double.parseDouble(sp[k+p]);

		return v;
	}
	
	public double getScalarData(String line){
		String[] sp=line.split(regex);	
		return Double.parseDouble(sp[sp.length-1]);
	}

	public int getIntData(String line){
		String[] sp=line.split(regex);	
		return Integer.parseInt(sp[sp.length-1]);
	}

	private boolean getBooleanData(String line){
		boolean b=false;
		String[] sp=line.split(regex);	
		
		if(sp[sp.length-1].startsWith("t") || sp[sp.length-1].equals("1"))	
			b=true;
		else 	if(sp[sp.length-1].startsWith("f") || sp[sp.length-1].equals("0"))	
			b=false;
		else {
			System.err.println(line);
			System.err.println(" Bad input!");
			
			wait(10000*10000);

			return false;
		}
		return b;

	}

	private String getStringData(String line){
		String[] sp=line.split(regex);	
		String[] sp2=sp[sp.length-1].split(" ");
		return sp2[sp2.length-1];
	}

	private void readAndSetRegMagPropery(Model model,int i,String line){
		int is=0;
		String[] sp=line.split(this.regex);	

		int ib=0;
		if(sp[0].equals("")) ib=1;
		int ir=Integer.parseInt(sp[ib++]);
	
		int BH_id=Integer.parseInt(sp[ib++]);
	
		model.region[ir].BHnumber=BH_id;
		
		Vect v=new Vect(model.dim);
		
		double mu=Double.parseDouble(sp[ib++]);
		
		for(int k=0;k<v.length;k++){
			v.el[k]=mu;
		}
	
	//	if(model.region[ir].BHnumber==0)
		model.region[ir].setMur(v);
		
		double sigma=Double.parseDouble(sp[ib++]);
		Vect v3=new Vect(3);
		for(int k=0;k<v3.length;k++){
			v3.el[k]=sigma;
		}

		model.region[ir].setSigma(v3);

		if(ib<sp.length){
			Vect M=new Vect(model.dim);
			for(int k=0;k<v.length;k++){
				M.el[k]=Double.parseDouble(sp[ib++]);;
			}
			model.region[ir].setM(M);
		}
		
		//=========== 
		if(model.region[ir].stranded) {


			model.stranded=true;
			
			if(model.dim==2)
			model.region[ir].windingSurf=model.getRegionArea(ir);
			else
			{
				model.region[ir].nloop=93;
			//	model.region[ir].windingSurf=model.getRegionAreaXY(ir);
				model.region[ir].windingSurf=model.getRegionVolume(ir)/.05;
			}
				
			/*
			model.region[ir].terminalVoltage0=Double.parseDouble(sp[is++]);
			model.region[ir].setFreq(Double.parseDouble(sp[is++]));

			model.region[ir].phase0=Double.parseDouble(sp[is++])*Math.PI/180;
			model.region[ir].setWireRes(Double.parseDouble(sp[is++]));

			if(is<sp.length)
				model.region[ir].nloop=Double.parseDouble(sp[is++]);
				else
					model.region[ir].nloop=100;

			
			if(is<sp.length)
				model.region[ir].circuit=sp[is++].startsWith("t");
			
			if(model.region[ir].circuit){

				if(is<sp.length)
				model.region[ir].curMap1=Integer.parseInt(sp[is++]);
			
			if(is<sp.length)
				model.region[ir].currCoef1=Double.parseDouble(sp[is++]);
			else
				model.region[ir].currCoef1=1;
			}*/
			
			model.region[ir].NtS=model.region[ir].nloop/model.region[ir].windingSurf;
		//	model.region[ir].NtS=462357.4142194;
 
			

		}

		
	}
	
	
private void readAndSetRegDataMech(Model model,int ir,String line){
		
		int is=0;
		String[] sp=line.split(this.regex);	


		is=1;
		
		model.region[ir].isotElast=this.getBooleanData(sp[is++]);
		
		model.region[ir].setRo(Double.parseDouble(sp[is++]));

		
		int mdim=model.dim;
		if(model.region[ir].isotElast)
			mdim=1;
		
		Vect E=new Vect(mdim);
		for(int k=0;k<mdim;k++)
			E.el[k]=Double.parseDouble(sp[is++]);
		
		model.region[ir].setYng(E);
		
		Vect v=new Vect(mdim);
		for(int k=0;k<mdim;k++)
			v.el[k]=Double.parseDouble(sp[is++]);
		
		model.region[ir].setPois(v);
		
		if(model.region[ir].isotElast){
			Vect sh=new Vect(1);
				sh.el[0]=.5*E.el[0]/(1+v.el[0]);
				
				model.region[ir].setShear(sh);

			}
		else{
			Vect sh=new Vect(model.dim);
			for(int k=0;k<mdim;k++)
				sh.el[k]=Double.parseDouble(sp[is++]);
				model.region[ir].setShear(sh);

			}


		if(is<sp.length)
			model.region[ir].setThermalCoef(Double.parseDouble(sp[is++]));
		if(is<sp.length)
			model.region[ir].setDeltaT(Double.parseDouble(sp[is++]));

		if(is<sp.length){
			model.nonLin=true;

			model.region[ir].setYield(Double.parseDouble(sp[is++]));

		}
		
		if(model.region[ir].getYield()>0)
		model.region[ir].setTangYoung(Double.parseDouble(sp[is++]));



	}


private void readAndSetRegDataSeep(Model model,int ir,String line){
	
	int is=0;
	String[] sp=line.split(this.regex);	

	model.region[ir].setName(sp[0]);

	is=1;
	Vect v=new Vect(model.dim);
	for(int k=0;k<v.length;k++){
		v.el[k]=Double.parseDouble(sp[is++]);
	}
	
	model.region[ir].setSigma(v);
}

public double[] loadArray(){
	String file=util.getFile();
	if(file==null || file.equals("") )  throw new NullPointerException("file not found.");
	return loadArray(file);
}
public int[][] loadIntArray(String arrayPath,int nRow,int nCol,int skip){
	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;

		
		int[][] data=new int[nRow][nCol];
		
		line=br.readLine();

		int i=0;
		while(line!=null && i<skip){
			i++;
			line=br.readLine();
			
		}
			
	//	util.pr(i);
		//util.pr(line);
		String[] sp;
		
		i=0;
		while(line!=null){
			if(i>nRow-1) break;
			sp=line.split(regex2);
			for(int k=0;k<sp.length;k++)
				data[i][k]=Integer.parseInt(sp[k]);
			
			i++;
			line=br.readLine();
			
		}
		
	br.close();
	fr.close();
			return data;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}

public double[] loadArray(String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;

		int N=100000;
		
		double[] x1=new double[N];
		
		int i=0;
		line=br.readLine();
		while(line!=null){
			if(i>N) break;
			x1[i++]=Double.parseDouble(line);
			line=br.readLine();
			
		}

		double[] x=Arrays.copyOf(x1, i);
		
	br.close();
	fr.close();
			return x;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}	


public Mat loadMatSymm(int n,String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
	
		Mat A=new Mat(n,n);
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<=i;j++){
				A.el[i][j]=x[j];	
				if(i!=j)
					A.el[j][i]=x[j];
			}
	
		}

		br.close();
		fr.close();
		
		return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}


public SpMat loadSparseMat(int dim,int K,String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;

		String[] sp;
	
		double val=0;
		SpMat As=new SpMat(dim);
		As.lower=true;
		int row=-1,col=0;
		int rowp;
		int k=0;

		while(true){
			line=br.readLine();
			sp=line.split(regex2);
			int ib=0;
			if(sp[0].equals("")) ib=1;
			rowp=row;
			row=Integer.parseInt(sp[ib]);
			if(row>rowp){
				
				if(rowp>=0){
				As.row[rowp].trim(k);
				}
				As.row[row]=new SpVect(dim,K);
				k=0;
			}
			col=Integer.parseInt(sp[ib+1]);
			val=Double.parseDouble(sp[ib+2]);
			As.row[row].index[k]=col;
			As.row[row].el[k]=val;
			k++;
			
			if(row==dim-1 && col==row){
				As.row[row].trim(k);
				 break;
			}
			}

		br.close();
		fr.close();
		
		return As;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}

public double[][] loadArrays(int n, int m,String arrayPath){
return loadArrays( n,  m, arrayPath,  0);
}

public double[][] loadArrays(int n, int m,String arrayPath, int skip){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		for(int i=0;i<skip;i++){
			line=br.readLine();
		}
		
		double[][] A=new double[n][m];
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<m;j++)
				A[i][j]=x[j];
	
			
			
		}

		br.close();
		fr.close();
			return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}

public Complex[] loadFFT(String file){

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
	
		line=br.readLine();
		int N=Integer.parseInt(line);;
		
		Complex[] x=new Complex[N];
		

		for(int i=0;i<N;i++){
			line=br.readLine();
			double[] v=this.getCSV(line);
	
			x[i]=new Complex(v[0],v[1]);
	
			
		}

		
			return x;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}	

public double[] getCSV(String line){
	
	String[] sp=line.split(regex);	

	int p0=0;
	if(sp[0].equals(""))
	{
		p0=1;
	}
	int L=sp.length-p0;

	double[] v=new double[L];

	for( int p=0;p<L;p++){

		v[p]=Double.parseDouble(sp[p+p0]);
	}

	return v;
}

public int[] getCSInt(String line){
	String[] sp=line.split(regex);	
	int L=sp.length;
	int[] v=new int[L];
	for( int p=0;p<L;p++)
				v[p]=Integer.parseInt(sp[p]);

	return v;
}

public String getNextDataLine(BufferedReader br) throws IOException{
	String line="";
	while(true){
		line=br.readLine();
		if(line==null) break;
		if(!line.startsWith("/")) break;
	}
return line;
}

public String getNextDataLine(BufferedReader br,String title) throws IOException{
	String line="";

	util.pr(title);

	
	while(true){
		line=br.readLine();
		if(line==null) break;
		
		if(line.equals("")) continue;

		if(!line.startsWith("/")) break;
	}
	util.pr(line);
	
return line;
}

public void setDataMag2D(Model model,BufferedReader br){

	String line;
	String s;
	int dim=model.dim;

	try {
		
		line=br.readLine();
		line=br.readLine();
		int coordCode =getIntData(line);
		model.coordCode=coordCode;
		line=br.readLine();
		line=br.readLine();
		int am =getIntData(line);
		model.analysisMode=am;

		line=br.readLine();
		line=br.readLine();
		model.motor=getBooleanData(line);
		line=br.readLine();
		line=br.readLine();
		
		int nRegions =getIntData(line);		
		if(nRegions!=model.numberOfRegions){
			System.err.println("Mesh and Data do not match in the number of regions: "+model.numberOfRegions+" and "+nRegions);
		}


		line=br.readLine();
		line=br.readLine();
	
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){
		line=br.readLine();
		line=br.readLine();
		readAndSetRegDataMag2D(model,ir,line);
		}

		model.BCtype=new int[model.nBoundary];
		for(int j=0;j<model.nBoundary;j++)
			model.BCtype[j]=-1;
		model.PBCpair=new int[model.nBoundary];
		line=br.readLine();
		line=br.readLine();

		int[] bcData=new int[2];

		for(int j=0;j<model.nBoundary;j++){
			line=br.readLine();


			if(model.BCtype[j]>-1) continue;
			bcData=getBCdata2D(line);
			model.BCtype[j]=bcData[0];

			model.PBCpair[j]=bcData[1];
		
			if(model.BCtype[j]>1){
			
				model.BCtype[model.PBCpair[j]]=bcData[0];
				model.PBCpair[model.PBCpair[j]]=j;

			 model.hasPBC=true;

			}

		
		}
		
		line=getNextDataLine(br);
		int numbRegsWithJ=Integer.parseInt(line);
		util.pr(numbRegsWithJ);
		for(int j=0;j<numbRegsWithJ;j++){
			line=getNextDataLine(br);
			String[] sp=line.split(this.regex);	

			int ib=0;
			if(sp[0].equals("")) ib=1;
			int nr=Integer.parseInt(sp[ib++]);
			Vect J=new Vect(dim);
			
			for(int k=0;k<dim;k++)
				J.el[k]=Double.parseDouble(sp[ib++]);
			util.pr(nr);
			model.region[nr].setJ(J);
			
		}

		line=getNextDataLine(br);
		model.hasBunif=getBooleanData(line);
		if(model.hasBunif){
			line=getNextDataLine(br);
			double[] array=getCSV(line);
			
			model.unifB=new Vect(array);
		}
		
		

		double f=0,dt=0;
		//int N=1;
	
		
		line=br.readLine();
		line=br.readLine();
		if(line!=null){
			 f =getScalarData(line);		
		}
		
		line=br.readLine();
		if(line!=null){
			 dt =getScalarData(line);		
		}
		line=br.readLine();
		if(line!=null){
			model.rotSpeed =getScalarData(line);		
		}
		line=br.readLine();
		if(line!=null){
			model.meshAngStep =getScalarData(line);		
		}
		line=br.readLine();

		if(line!=null){
			
			int nSteps=1;
			String sp[]=line.split(regex);
			int L=sp.length;
			
			int n1=0,n2=0, d=1;
	
				n1=Integer.parseInt(sp[L-3]);
				n2=Integer.parseInt(sp[L-2]);
				d=Integer.parseInt(sp[L-1]);
				if(d!=0)
				 nSteps=(n2-n1)/d+1;
		
		
		
				model.setnTsteps(nSteps);
				model.nBegin=n1;
				model.nEnd=n2;
				model.nInc=d;

			}
		
	
	
		line=br.readLine();
		model.eddyTimeIntegMode=getIntData(line);	
		

	


for(int j=0;j<50+model.numberOfRegions;j++){
	
		line=br.readLine();
	
		if(line==null) continue;
		if(line.startsWith("MS")){
					
			String sp[]=line.split(regex);
			int L=sp.length;
							
				int nr=Integer.parseInt(sp[L-3]);
				model.region[nr].MS=true;

				Vect E=new Vect(Double.parseDouble(sp[L-2]),0,0);
				model.region[nr].setYng(E);
				Vect v=new Vect(Double.parseDouble(sp[L-1]),0,0);
				model.region[nr].setPois(v);	
				
				if(model.region[nr].isotElast){
				Vect sh=new Vect(1);
					v.el[0]=.5*E.el[0]/(1+v.el[0]);
					
					model.region[nr].setShear(sh);

				}
				

			}
		else if(line.startsWith("loadFlux")) 
			model.loadFlux=this.getBooleanData(line);
		else if(line.startsWith("fluxFolder")){
			line=br.readLine();
			model.fluxFolderIn=line;
		}
		else if(line.startsWith("saveFlux")) model.saveFlux=this.getBooleanData(line);
		else if(line.startsWith("saveForce")) model.saveForce=this.getBooleanData(line);
		else if(line.startsWith("mag")) model.magAnalysis=this.getBooleanData(line);
		else if(line.startsWith("trans")) model.transfer2DTo3D=this.getBooleanData(line);	
		else if(line.startsWith("forceCal")) model.forceCalcMode=this.getIntData(line);	
		else if(line.startsWith("rotate")) model.rotateRotor=this.getBooleanData(line);	
		else if(line.startsWith("loadPrev")) model.loadPrevMag=this.getBooleanData(line);
		else if(line.startsWith("axi")) model.axiSym=this.getBooleanData(line);	
		else if(line.startsWith("height")) model.height=this.getScalarData(line);
		else if(line.startsWith("POD")) model.POD=this.getIntData(line);
		else if(line.startsWith("snapShot")) model.snapShot=this.getIntData(line);
		

	}








if(model.axiSym) model.height=2*Math.PI;
			
		model.setFreq(f);
		model.setDt(dt*model.nInc);
		
		model.setHasJ();

		model.setHasM();

		model.setHasMS();


		model.setNonLinToElememts();
		
		model.setEdge();


		model.setElementsParam();

		
		model.setBounds();

	
	if(model.coordCode==1) {
		
		model.cpb=1;

		for(int j=0;j<model.nBoundary;j++){

			if(model.BCtype[j]==3) model.cpb=-1;
		}
		if(model.hasTwoNodeNumb && model.rotateConnect && model.nRotReg>0)
		 model.mapCommonNodes();	
	}
	
	//=====================


	model.setNodeOnBound();
	
	if(model.hasPBC) model.mapPBC();

	
		if(model.deform)
			model.setMechBC();
		
		
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(model.region[ir].circuit ) {
				model.circuit=true;
				break;
			}
			
		}

	
		int ix=0;
		int iy=0;

		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(model.region[ir].circuit ) {
				model.region[ir].currentIndex=ix;
		ix++;
		}

		}
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(model.region[ir].circuit && model.region[ir].curMap1==0)
			{
			model.region[ir].unknownCurrentIndex=iy;
		iy++;
		}

		}
		
		
		model.numberOfCurrents=ix;
		
		model.numberOfUnknownCurrents=iy;
		
		
		if(ix>0){/*
			model.analysisMode=1;
			model.eddyTimeIntegMode=-2;
		*/}
	//--------------
		
		
		model.threePhaseRegs=new int[3];
		if(model.circuit){


			if(model.numberOfRegions>8){
				model.threePhaseRegs[0]=9;
				model.threePhaseRegs[1]=11;
				model.threePhaseRegs[2]=13;
			}
				
				 if(model.numberOfRegions==8){
					model.threePhaseRegs[0]=1;
					model.threePhaseRegs[1]=3;
					model.threePhaseRegs[2]=5;
				}
				 if(model.numberOfRegions==4){
						model.threePhaseRegs[0]=1;
						model.threePhaseRegs[1]=0;
						model.threePhaseRegs[2]=0;
					}

			}
		
	//)))))))))))))))))))))))

		if(model.threePhaseRegs[0]!=0 &&model.threePhaseRegs[1]!=0&& model.threePhaseRegs[2]!=0)
			model.nNeutral=1;
		
		//model.nNeutral=0;

	
		int nCur=model.numberOfUnknownCurrents;
		
		model.unCurRegNumb=new int[nCur];
	
		for(int ir=1;ir<=model.numberOfRegions;ir++)
			if(model.region[ir].circuit && model.region[ir].curMap1==0)
				model.unCurRegNumb[model.region[ir].unknownCurrentIndex]=ir;
		


		try {
			model.ia=new CurrentWaveForm("emf//Ian.txt");
			model.ib=new CurrentWaveForm("emf//Ibn.txt");
			model.ic=new CurrentWaveForm("emf//Icn.txt");
		/*	model.ia=new CurrentWaveForm("emf//IaMath.txt");
			model.ib=new CurrentWaveForm("emf//IbMath.txt");
			model.ic=new CurrentWaveForm("emf//IcMath.txt");

			/*	model.va=new CurrentWaveForm("emf//Va.txt");
			model.vb=new CurrentWaveForm("emf//Vb.txt");
			model.vc=new CurrentWaveForm("emf//Vc.txt");*/
		} catch (Exception e) {
			e.printStackTrace();
		}
	
	
	
		
		
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	
	if(model.nonLin)
		model.setNonLin(true);
	
	
	
	//model.setForceCalc();

}

private void readAndSetRegDataMag2D(Model model,int ir,String line){
	int is=0;
	String[] sp=line.split(this.regex);	
	


	model.region[ir].setName(sp[0]);
	model.region[ir].setMaterial(sp[1]);
	model.region[ir].setNonLinear(sp[2].startsWith("t"));
	
	if(model.region[ir].isNonLinear) model.nonLin=true;

	is=3;
	Vect v=new Vect(model.dim);
	for(int k=0;k<v.length;k++){
		v.el[k]=Double.parseDouble(sp[is++]);
	}
	
	if(!model.region[ir].isNonLinear)
	model.region[ir].setMur(v);

	if(is==sp.length) return;

	v=new Vect(model.dim);
		for(int k=0;k<v.length;k++)
		v.el[k]=Double.parseDouble(sp[is++]);
	model.region[ir].setM(v);
	

	
	if(is==sp.length) return;
	

	model.region[ir].stranded=sp[is++].startsWith("t");
	


if(model.region[ir].stranded) {


	model.stranded=true;
	
	if(model.dim==2)
	model.region[ir].windingSurf=model.getRegionArea(ir);
	else{
	}
	
	model.region[ir].terminalVoltage0=Double.parseDouble(sp[is++]);
	model.region[ir].setFreq(Double.parseDouble(sp[is++]));

	model.region[ir].phase0=Double.parseDouble(sp[is++])*Math.PI/180;
	model.region[ir].setWireRes(Double.parseDouble(sp[is++]));

	if(is<sp.length)
		model.region[ir].nloop=Double.parseDouble(sp[is++]);
		else
			model.region[ir].nloop=100;

	
	if(is<sp.length)
		model.region[ir].circuit=sp[is++].startsWith("t");
	
	if(model.region[ir].circuit){

		if(is<sp.length)
		model.region[ir].curMap1=Integer.parseInt(sp[is++]);
	
	if(is<sp.length)
		model.region[ir].currCoef1=Double.parseDouble(sp[is++]);
	else
		model.region[ir].currCoef1=1;
	}
	model.region[ir].NtS=model.region[ir].nloop/model.region[ir].windingSurf;

	

}
else{

	if(is==sp.length) return;

	if(model.dim==3){
		for(int k=0;k<v.length;k++)
			v.el[k]=Double.parseDouble(sp[is++]);

		}
		else{

			v=new Vect(0,0,Double.parseDouble(sp[is++]));
			}
	
		model.region[ir].setJ(v);


/*
	
		if(is==sp.length) return;
		
		model.region[ir].setFreq(Double.parseDouble(sp[is++]));
		model.region[ir].phase0=Double.parseDouble(sp[is++])*Math.PI/180;*/

		
		if(is==sp.length) return;
		model.region[ir].setSigma(new Vect(0,0,Double.parseDouble(sp[is++])));
}


}


public void wait(int ms){
	try {
		Thread.sleep(ms);
	} catch (InterruptedException e) {
		e.printStackTrace();
	}
}


}
