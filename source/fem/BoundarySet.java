package fem;


import math.IntVect;
import math.Mat;
import math.Vect;
import math.util;
import static java.lang.Math.*;

import java.util.Arrays;

public class BoundarySet {

	double epsr=1e-6,epsAng=1e-4;

	public BoundarySet(){	}

	public void magBC(Model model){

		if(model.coordCode==1) {
			setSliceBounds(model);

			mapCommonNodes(model);	

		}

		setNodeOnBound(model);

		if(model.hasPBC) mapPBC(model);


		
		setMagBC(model);



	}

	private void mapEdges(Model model){

		if(model.dim==3) {mapEdges3D(model); return;}
		int[] nodeEdge=new int[1+model.numberOfNodes];
		for(int i=1;i<=model.numberOfEdges;i++){
			nodeEdge[model.edge[i].node[0].id]=i;
		}
		int end0;

	
		for(int i=1;i<=model.numberOfEdges;i++){

			end0=model.edge[i].node[0].id;

			int nmp=model.node[end0].getMap();
			if(nmp>0){
				int nmp2=model.node[nmp].getMap();
				if(nmp2==0)
					model.edge[i].map=nodeEdge[nmp];
				else
					model.edge[i].map=nodeEdge[nmp2];	

				if(model.node[end0].sPBC ) 
					model.edge[i].sPBC=true;
				else if(model.node[end0].aPBC ) {
					model.edge[i].aPBC=true;
				}
				else {
					model.edge[i].common=true;
				}

			}	


		}


	}

	private void mapEdges3D(Model model){

		IntVect[] nodeEdge=new IntVect[1+model.numberOfNodes];
		int[] ic=new int[1+model.numberOfNodes];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodeEdge[i]=new IntVect(12);

		int end0,end1;

		for(int i=1;i<=model.numberOfEdges;i++){
			end0=model.edge[i].node[0].id;
			end1=model.edge[i].node[1].id;
			nodeEdge[end0].el[ic[end0]++]=i;	
			nodeEdge[end1].el[ic[end1]++]=i;
		}




		for(int i=1;i<=model.numberOfEdges;i++){
			end0=model.edge[i].node[0].id;



			int nmp0=model.node[end0].getMap();
			if(nmp0>0){
				end1=model.edge[i].node[1].id;
				int nmp1=model.node[end1].getMap();

				int emap=0;
				for(int j=0;j<ic[nmp0];j++)
					for(int k=0;k<ic[nmp1];k++)
						if(nodeEdge[nmp0].el[j]==nodeEdge[nmp1].el[k]){emap=nodeEdge[nmp0].el[j]; break;}

				if(emap==0) continue;
				//if(model.node[end1].getCoord().v2().sub(model.node[end0].getCoord().v2()).norm()>1e-4) continue;

				model.edge[i].map=emap;	

				if(model.node[end0].sPBC ) 
					model.edge[i].sPBC=true;
				else if(model.node[end0].aPBC ) {
					model.edge[i].aPBC=true;
				}
				else {
					model.edge[i].common=true;
				}

			}	

		}
		



	}

	public void setTBC(Model model){


		if(model.hasPBC || model.hasTwoNodeNumb)
			mapEdges(model);
	


		if(model.dim==3) {


			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].hasJ)
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++) {
						int[] en=model.element[i].getEdgeNumb();
						for(int j=0;j<model.nElEdge;j++){
							if(!model.edge[en[j]].hasJ){
								model.edge[en[j]].hasJ=true;

							}

						}


					}
			}


			double Au0=-1e6;
			
			byte b0=0,b1=1,b2=2;

			
			int[] edgeDir=new int[model.numberOfEdges+1];
			for(int i=1;i<=model.numberOfElements;i++){
				int[] edgeNumb=model.element[i].getEdgeNumb();
				
			
				for(int j=0;j<model.nElEdge;j++)
					if(j<4){
						edgeDir[edgeNumb[j]]=b0;
						}
					else if(j<8)
						edgeDir[edgeNumb[j]]=b1;
					else{
						edgeDir[edgeNumb[j]]=b2;
					}
			}
		

			for(int i=1;i<=model.numberOfEdges;i++){

				if(!model.edge[i].hasJ)continue;

				if(model.edge[i].node[0].getCoord(2)<.01251 && 
						model.edge[i].node[1].getCoord(2)<.01251){
					model.edge[i].edgeKnownT=true;
					if(	edgeDir[i]==1){

						model.edge[i].T=Au0*model.edge[i].edgeLength;


					}
					else
						model.edge[i].T=0;


				}
				else if(model.edge[i].node[0].getCoord(2)>.107 && 
						model.edge[i].node[1].getCoord(2)>.105){
					model.edge[i].edgeKnownT=true;

					model.edge[i].T=0;

				}
				else if(model.edge[i].node[0].getCoord().v2().norm()<.01201 && 
						model.edge[i].node[1].getCoord().v2().norm()<.01201){
					model.edge[i].edgeKnownT=true;
					model.edge[i].T=0;

				}
				else if(model.edge[i].node[0].getCoord().v2().norm()>.03299 && 
						model.edge[i].node[1].getCoord().v2().norm()>.03299){
					model.edge[i].edgeKnownT=true;
					model.edge[i].T=0;

				}
			}
		}

		setTIndice(model);


	}


	public void mechBC(Model model){

		int dim=model.dim;

		if(!model.deform) return;

		for(int ir=1;ir<=model.numberOfRegions;ir++){


			if(model.region[ir].getYng().el[0]<1.0) continue;

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
			

				model.element[i].setDeformable(true);

				if(model.element[i].isDeformable())
				for(int j=0;j<model.nElVert;j++){
					int nNode=model.element[i].getVertNumb(j);
					model.node[nNode].setDeformable(true);

				}
			}

		}

		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].isDeformable() ) {
				model.node[i].u=new Vect(dim);

			}

		int ix=0;


		double[][] bounds=new double[model.mechBC.length][model.nBoundary];

		Vect[] val=new Vect[model.mechBC.length];
		boolean[] regBound=new boolean[model.mechBC.length];
		boolean[] nodeConst=new boolean[model.mechBC.length];

		boolean[][] hasVal=new boolean[model.mechBC.length][dim];

		for(int j=0;j<model.mechBC.length;j++)
			for(int k=0;k<model.nBoundary/2;k++){
				bounds[j][2*k]=-1e10;
				bounds[j][2*k+1]=1e10;
			}
		int[][] boundCoSys=new int[model.mechBC.length][dim];

		for(int j=0;j<model.mechBC.length;j++){

			//util.hshow(model.mechBC[j]);;

			val[j]=new Vect(dim);

			if(model.mechBC[j][0]==null) break;

			for(int k=0;k<model.mechBC[0].length;k++){

				if(model.mechBC[j][k]==null) continue;

				if(model.mechBC[j][k].startsWith("reg"))

				{
					k++;

					int nr=Integer.parseInt(model.mechBC[j][k]);
					int[] nn=model.getRegNodes(nr);

					k++;
					if(model.mechBC[j][k].equals("0") )
					{
						k++;
						val[j].el[0]=Double.parseDouble(model.mechBC[j][k]);

						hasVal[j][0]=true;
			

					}

					k++;
					if(model.mechBC[j][k].equals("1") )
					{
						k++;
						val[j].el[1]=Double.parseDouble(model.mechBC[j][k+1]);
		
						hasVal[j][1]=true;
					}

					k++;
					if(model.mechBC[j][k].equals("2") )
					{
						k++;
						val[j].el[2]=Double.parseDouble(model.mechBC[j][k+1]);
		
						hasVal[j][2]=true;
					}

					for(int m=0;m<nn.length;m++){
						for(int p=0;p<dim;p++)
							if(hasVal[j][p]){
								model.node[nn[m]].setKnownU(p,val[j].el[p]);
							}
						ix++;


					}

					regBound[j]=true;


				}
				if(model.mechBC[j][k].startsWith("node"))

				{
					k++;

					int nid=Integer.parseInt(model.mechBC[j][k]);

					k++;
					if(model.mechBC[j][k].equals("0") )
					{
						k++;
						val[j].el[0]=Double.parseDouble(model.mechBC[j][k]);
			
						model.node[nid].setKnownU(0,val[j].el[0]);
				
						util.pr(nid+" 0 "+val[j].el[0]);
					}
					k++;
					if(model.mechBC[j][k].equals("1") )
					{
			
						k++;
						val[j].el[1]=Double.parseDouble(model.mechBC[j][k]);
						util.pr(nid+" 1 "+val[j].el[1]);
						model.node[nid].setKnownU(1,val[j].el[1]);
					
					}

					k++;
					if(model.dim==3 && model.mechBC[j][k].equals("2") )
					{
						k++;
						val[j].el[2]=Double.parseDouble(model.mechBC[j][k]);
												
						model.node[nid].setKnownU(2,val[j].el[2]);
					}

					ix++;
					nodeConst[j]=true;
				}


				if(!regBound[j] && !nodeConst[j]){
		
					if(model.mechBC[j][k].equals("x") || model.mechBC[j][k].equals("r"))
					{
						bounds[j][0]=Double.parseDouble(model.mechBC[j][k+1]);
						bounds[j][1]=Double.parseDouble(model.mechBC[j][k+2]);
						if(model.mechBC[j][k].equals("x")) 
							boundCoSys[j][0]=0;
						else
							boundCoSys[j][0]=1;

						k+=2;

					}
					else if(model.mechBC[j][k].equals("y") || model.mechBC[j][k].equals("t"))
					{
						bounds[j][2]=Double.parseDouble(model.mechBC[j][k+1]);
						bounds[j][3]=Double.parseDouble(model.mechBC[j][k+2]);

						if(model.mechBC[j][k].equals("y")) 
							boundCoSys[j][1]=0;
						else
							boundCoSys[j][1]=1;

						k+=2;
					}
					else if(dim==3 && model.mechBC[j][k].equals("z"))
					{
						bounds[j][4]=Double.parseDouble(model.mechBC[j][k+1]);
						bounds[j][5]=Double.parseDouble(model.mechBC[j][k+2]);
		
						if(model.mechBC[j][k].equals("z")) 
							boundCoSys[j][2]=0;
						else
							boundCoSys[j][2]=1;

						k+=2;
					}
					else
					{
						if(model.mechBC[j][k].equals("0") )
						{
							val[j].el[0]=Double.parseDouble(model.mechBC[j][k+1]);
							k++;
							hasVal[j][0]=true;
						}
						else if(model.mechBC[j][k].equals("1") )
						{
				
							val[j].el[1]=Double.parseDouble(model.mechBC[j][k+1]);
							k++;
							hasVal[j][1]=true;
						}
						else if(model.mechBC[j][k].equals("2"))
						{
							val[j].el[2]=Double.parseDouble(model.mechBC[j][k+1]);
							k++;
							hasVal[j][2]=true;
						}


					}
				}
			}
		}




		for(int i=1;i<=model.numberOfNodes;i++){

			if(!model.node[i].isDeformable()) continue;
			Vect v=model.node[i].getCoord();
			Vect v2=v.v2();

			Vect vv=new Vect(dim);


			for(int j=0;j<model.mechBC.length;j++){

				if(regBound[j]) continue;

				if(boundCoSys[j][0]==0) 
					vv.el[0]=v.el[0];
				else if(boundCoSys[j][0]==1)					
					vv.el[0]=v2.norm();

				if(boundCoSys[j][1]==0)
					vv.el[1]=v.el[1];
				else if(boundCoSys[j][1]==1) {
					vv.el[1]=util.getAng(v2);
				}

				if(dim==3){
					vv.el[2]=v.el[2];
				}


				boolean contained=true;


				for(int k=0;k<dim;k++)
					if(vv.el[k]<bounds[j][2*k]||vv.el[k]>bounds[j][2*k+1])
						contained=false;

				if(contained){

					for(int k=0;k<dim;k++)
						if(hasVal[j][k]){
							if(!model.node[i].is_U_known(k))

							{
								model.node[i].setKnownU(k,val[j].el[k]);
								ix++;

							}
						}
				}

			}
		}	

	

		//==========================

		System.out.println(" Number of known u components : "+(ix));

		if(model.hasPBC || model.hasTwoNodeNumb)
			for(int i=1;i<=1*model.numberOfNodes;i++){
				int nn=model.node[i].getMap();
				if(nn>0 && model.node[nn].has_U_known())
					for(int k=0;k<dim;k++)
						model.node[i].setU_is_known(k, model.node[nn].is_U_known(k));
			}

		setMechIndice(model);
		System.out.println("Numb of unknown displacements: "+model.numberOfUnknownU);
		System.out.println(" Number of known u components : "+(model.numberOfNodes*dim-ix));

	

	}



	public void seepBC(Model model){

		int dim=model.dim;


		for(int ir=1;ir<=model.numberOfRegions;ir++){


			if(!model.region[ir].isConductor) continue;


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;
				for(int j=0;j<model.nElVert;j++){
					int nNode=model.element[i].getVertNumb(j);
					model.node[nNode].setPhiVar(true);

				}
			}

		}

		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].isPhiVar()) {
				model.node[i].T=0;

			}

		int ix=0;


		double[][] bounds=new double[model.seepBC.length][model.nBoundary];

		double[] val=new double[model.seepBC.length];
		boolean[] hasVal=new boolean[model.seepBC.length];

		for(int j=0;j<model.seepBC.length;j++)
			for(int k=0;k<model.nBoundary/2;k++){
				bounds[j][2*k]=-1e10;
				bounds[j][2*k+1]=1e10;
			}

		int[][] boundCoSys=new int[model.seepBC.length][dim];

		for(int j=0;j<model.seepBC.length;j++){

			if(model.seepBC[j][0]==null) break;

			for(int k=0;k<model.seepBC[0].length;k++){

				if(model.seepBC[j][k]==null) continue;

				if(model.seepBC[j][k].equals("x") || model.seepBC[j][k].equals("r"))
				{
					bounds[j][0]=Double.parseDouble(model.seepBC[j][k+1]);
					bounds[j][1]=Double.parseDouble(model.seepBC[j][k+2]);
					if(model.seepBC[j][k].equals("x")) 
						boundCoSys[j][0]=0;
					else
						boundCoSys[j][0]=1;

					k+=2;

				}
				else if(model.seepBC[j][k].equals("y") || model.seepBC[j][k].equals("t"))
				{
					bounds[j][2]=Double.parseDouble(model.seepBC[j][k+1]);
					bounds[j][3]=Double.parseDouble(model.seepBC[j][k+2]);

					if(model.seepBC[j][k].equals("y")) 
						boundCoSys[j][1]=0;
					else
						boundCoSys[j][1]=1;

					k+=2;
				}
				else if(dim==3 && model.seepBC[j][k].equals("z"))
				{
					bounds[j][4]=Double.parseDouble(model.seepBC[j][k+1]);
					bounds[j][5]=Double.parseDouble(model.seepBC[j][k+2]);

					if(model.seepBC[j][k].equals("z")) 
						boundCoSys[j][2]=0;
					else
						boundCoSys[j][2]=1;

					k+=2;
				}
				else
				{	

					val[j]=Double.parseDouble(model.seepBC[j][k]);
					hasVal[j]=true;
				}
			}
		}

		for(int i=1;i<=model.numberOfNodes;i++){

			if(!model.node[i].isPhiVar()) continue;
			Vect v=model.node[i].getCoord();
			Vect v2=v.v2();
			Vect vv=new Vect(dim);

			for(int j=0;j<model.seepBC.length;j++){

				if(model.seepBC[j][0]==null) continue;

				if(boundCoSys[j][0]==0) 
					vv.el[0]=v.el[0];
				else if(boundCoSys[j][0]==1)					
					vv.el[0]=v2.norm();

				if(boundCoSys[j][1]==0)
					vv.el[1]=v.el[1];
				else if(boundCoSys[j][1]==1) 
					vv.el[1]=util.getAng(v2);

				if(dim==3){
					vv.el[2]=v.el[2];
				}


				boolean contained=true;

				for(int k=0;k<dim;k++)
					if(vv.el[k]<bounds[j][2*k]||vv.el[k]>bounds[j][2*k+1])
						contained=false;


				if(contained){

					if(hasVal[j]){

						if(!model.node[i].is_T_known())

						{
							model.node[i].setKnownT(val[j]);
							ix++;

						}
					}
				}

			}
		}		

		//==========================

		System.out.println(" Number of known T : "+(ix));

		if(model.hasPBC || model.hasTwoNodeNumb)
			for(int i=1;i<=1*model.numberOfNodes;i++){
				int nn=model.node[i].getMap();
				if(nn>0 && model.node[nn].is_T_known())
					model.node[i].setU_is_known(model.node[nn].is_U_known());
			}


		setSeepIndice(model);
		System.out.println("Numb of unknown T: "+model.numberOfUnknownT);

	}


	public void setMagIndice(Model model){
		
		
		if(model.analysisMode==2){
			
		

			model.nodeVarIndex=new int[model.numberOfNodes+1];

			for(int ir=1;ir<=model.numberOfRegions;ir++)
				if(model.region[ir].isConductor)
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
						int[] vertNumb=model.element[i].getVertNumb();
						for(int j=0;j<model.nElVert;j++){
							model.node[vertNumb[j]].setPhiVar(true);
						}
					}
			
			for(int i=1;i<=model.numberOfEdges;i++)
				if(model.edge[i].edgeKnown){
					model.edge[i].node[0].setPhiVar(false);
					model.edge[i].node[1].setPhiVar(false);
				}

			
			
			int nnx=0;	
				for(int i=1;i<=model.numberOfNodes;i++)
					if(model.node[i].isPhiKnown())
						nnx++;
				
			model.numberOfKnownPhis=nnx;
		

			 nnx=0;	
			for(int i=1;i<=model.numberOfNodes;i++)
				if(model.node[i].isPhiVar())
					model.nodeVarIndex[i]=++nnx;

			model.numberOfVarNodes=nnx;
	
			
			model.varNodeNumber=new int[model.numberOfVarNodes+1];

			nnx=0;

			for(int i=1;i<=model.numberOfNodes;i++){
				if(model.node[i].isPhiVar() && !model.node[i].isPhiKnown()){
					model.varNodeNumber[++nnx]=i;
				}
			}

		}

	
		boolean[] rotorEdges=null;
		if(model.motor){
		rotorEdges=new boolean[1+model.numberOfEdges];
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(model.region[ir].rotor){
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					int[] ne=model.element[i].getEdgeNumb();
					for(int k=0;k<ne.length;k++){
						if(!rotorEdges[ne[k]]){
							rotorEdges[ne[k]]=true;
						}
					}
			}
		}
		}
		}

				
		int nex=0;	
		model.edgeUnknownIndex=new int[model.numberOfEdges+1];
		
		
		if(rotorEdges!=null){
			for(int i=1;i<=model.numberOfEdges;i++){
				if(rotorEdges[i] &&!model.edge[i].edgeKnown && model.edge[i].map==0){
					model.edgeUnknownIndex[i]=++nex;
					
						}
			}
			
			for(int i=1;i<=model.numberOfEdges;i++){
				if(!rotorEdges[i] &&!model.edge[i].edgeKnown && model.edge[i].map==0){
					model.edgeUnknownIndex[i]=++nex;
					
						}
			}

		}

		else{

			for(int i=1;i<=model.numberOfEdges;i++){
				if(!model.edge[i].edgeKnown && model.edge[i].map==0){
					model.edgeUnknownIndex[i]=++nex;
						}
			}
			

		}

		for(int i=1;i<=model.numberOfEdges;i++){

			if(model.edge[i].map>0){
				model.edgeUnknownIndex[i]=model.edgeUnknownIndex[model.edge[i].map];

			}
		}

		int kx=0;
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edge[i].map>0) kx++;

		}
		

		//}

		model.numberOfMappedEdges=kx;

		model.numberOfUnknownEdges=nex;


		model.numberOfKnownEdges=model.numberOfEdges-model.numberOfUnknownEdges-kx;


		model.numberOfUnknowns=model.numberOfUnknownEdges+
		model.numberOfVarNodes+model.numberOfUnknownCurrents+model.nNeutral;

		model.unknownEdgeNumber=new int[model.numberOfUnknownEdges+1];
		model.knownEdgeNumber=new int[model.numberOfKnownEdges+1];
		model.knownEdgeValue=new double[model.numberOfKnownEdges+1];

		nex=0;
		int neknown=0;
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edge[i].map>0) continue;
			if(!model.edge[i].edgeKnown )
				model.unknownEdgeNumber[++nex]=i;
			else if(model.edge[i].edgeKnown){
				neknown++;
				model.knownEdgeNumber[neknown]=i;
				model.knownEdgeValue[neknown]=model.edge[i].A;

			}
		}


		

	}


	public void setMagBC(Model model){
		
		detectFarboundaryEdges(model);

		
if(model.hasTwoNodeNumb &&  !model.rotateConnect)
{


	//model.rm=1.6250000;
		//rm=1.8333333;
		//rm=1.6250000;
	double rm=model.rm;//0.0547;
	util.pr("sliding surface radius = "+model.rm);
	//double rm=0.0547;
	int jx=0;
	int[] edgeOnFSIndices=new int[1+model.numberOfEdges];

	for(int i=1;i<=model.numberOfEdges;i++){
		edgeOnFSIndices[i]=-1;

		
		if(abs(model.edge[i].node[0].getR()-rm)<1e-5){
			edgeOnFSIndices[i]=jx;
		
			jx++;
		}
	}


	model.edgeOnFSIndices=edgeOnFSIndices;


	for(int i=1;i<=model.numberOfEdges;i++){

	if(model.edgeOnFSIndices[i]>=0) model.edge[i].setKnownA(0);
	}
	}
	//===========

		if(model.hasPBC || model.hasTwoNodeNumb){

			mapEdges(model);
		}

		if(model.hasJ || model.hasM || model.stranded){
			
			int nNeumann=0;
			int nDirichlet=0;
			for(int i=1;i<=model.numberOfEdges;i++){
				for(int j=0;j<model.nBoundary;j++){

					if(model.BCtype[j]==1 && model.edge[i].node[0].onBound[j] && 
							model.edge[i].node[1].onBound[j]){
						model.edge[i].setKnownA(0);
						nDirichlet++;
					//	break;					
					}		
					else if(model.BCtype[j]==0 && model.edge[i].node[0].onBound[j] && 
							model.edge[i].node[1].onBound[j]){
						model.edge[i].edgeKnown=false;
						nNeumann++;
						//break;					
					}			
				}
			}
			
			util.pr("Number of edges on Bn0 boundary = "+nDirichlet);
			util.pr("Number of edges on Ht0 boundary = "+nNeumann);

		}
		else if(model.hasBunif){
			model.magMat.setMagBCUniform(model);
		}
		


		for(int i=1;i<=model.numberOfEdges;i++){
			int mp=model.edge[i].map;
			if(mp>0){
				model.edge[i].edgeKnown=model.edge[mp].edgeKnown;
			}
		}


		setMagIndice(model);

	}


	public void setTIndice(Model model){


		int ix=0;                  
		for(int i=1;i<=model.numberOfEdges;i++){

			if(!model.edge[i].hasJ || model.edge[i].edgeKnownT ) continue;

			ix++;

		}

		model.numberOfUnknownT=ix;	
		model.T_unknownIndex=new int[model.numberOfEdges+1];
		int nex=0;
		for(int i=1;i<=model.numberOfEdges;i++){
			if(model.edge[i].hasJ && !model.edge[i].edgeKnownT && model.edge[i].map==0){
				model.T_unknownIndex[i]=++nex;

			}
		}



	}


	public void setMechIndice(Model model){

		int ix=1;
		int jx=0;
		model.U_unknownIndex=new int[model.numberOfNodes+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].isDeformable() && !model.node[i].is_U_known() && model.node[i].getMap()==0){
				model.U_unknownIndex[i]=ix++;
				for(int k=0;k<model.dim;k++){
					if(!model.node[i].is_U_known(k)) jx++;

				}
			}


		for(int i=1;i<=model.numberOfNodes;i++){

			if(model.node[i].getMap()>0){

				model.U_unknownIndex[i]=model.U_unknownIndex[model.node[i].getMap()];

			}
		}


		model.numberOfUnknownU=ix-1;
		model.numberOfUnknownUcomp=jx;

		model.unknownUnumber=new int[model.numberOfUnknownU+1];

		ix=1;
		for(int i=1;i<=model.numberOfNodes;i++){
			if(model.node[i].getMap()>0) continue;
			if(model.U_unknownIndex[i]>0)
				model.unknownUnumber[ix++]=i;
		}


	}


	public void setSeepIndice(Model model){

		int ix=1;
		model.T_unknownIndex=new int[model.numberOfNodes+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].isPhiVar() && !model.node[i].is_T_known() && model.node[i].getMap()==0)
			{


				model.T_unknownIndex[i]=ix++;


			}


		for(int i=1;i<=model.numberOfNodes;i++){

			if(model.node[i].getMap()>0){
				model.T_unknownIndex[i]=model.T_unknownIndex[model.node[i].getMap()];

			}
		}


		model.numberOfUnknownT=ix-1;


	}


	public void mapPBC(Model model){

		boolean[] PBdone=new boolean[model.nBoundary];
		for(int j=0;j<model.nBoundary;j++){
			if(model.BCtype[j]>1 && !PBdone[j]){
				PBdone[model.PBCpair[j]]=true;
				mapPBC(model,j);
			}
		}
	}

	public void mapPBC(Model model, int nb){

		//something wromg about mapping bc next stetps
		int[][] mappedNode=mapBorderNodes(model,nb,model.PBCpair[nb]);
		
		for(int i=0;i<mappedNode.length;i++){

			int nmp=mappedNode[i][1];
			model.node[nmp].setPBC(model.cpb);
			if(model.node[mappedNode[i][0]].getMap()==0){
				model.node[nmp].setMap(mappedNode[i][0]);


			}
			else
				model.node[nmp].setMap(model.node[mappedNode[i][0]].getMap());
		}


	}


	public int[][] mapBorderNodes(Model model,int nb1,int nb2){


		int[][] commonNode=new int[1][1];
		if(model.dim==2){
			commonNode=mapBorderNodes2D(model,nb1,nb2);
		}
		else
			commonNode=mapBorderNodes3D(model,nb1,nb2);


		return commonNode;

	}
	
	public int[] getBorderNodeSorted(Model model,int nb){

		int nNodes=model.numberOfNodes;
		int ix=0;
		for(int i=1;i<=nNodes;i++){

			if(model.node[i].onBound[nb]) {
				ix++;
			}
		}

		int L=ix;

		int[] boundNodeNumb=new int[L+1];
		int[] boundNodeInd=new int[nNodes+1];
		Vect v=new Vect(L+1);
		v.el[0]=-1000;
		ix=0;

		for(int i=1;i<=nNodes;i++)
			if(model.node[i].onBound[nb])  {
				
				Vect r0=model.node[i].getCoord();
				Vect r=r0.v2();

				ix++;
				boundNodeInd[i]=ix;
				boundNodeNumb[ix]=i;
				if(model.motor){
					if(nb<2)
						v.el[ix]=util.getAng(r);
					else 
						if(nb<4) v.el[ix]=r.norm();

				}
				else{
					if(nb<2)
						v.el[ix]=r.el[1];
					else if(nb<4)
						v.el[ix]=r.el[0];
				}
			}

		int[] indSorted=v.bubble();

		if(nb<2)
			for(int i=2;i<=L;i++)
				if(abs(v.el[i]-v.el[i-1])<1e-4) 
					if(model.node[boundNodeNumb[indSorted[i]]].rotor){
						int temp=indSorted[i];
						indSorted[i]=indSorted[i-1];
						indSorted[i-1]=temp;
					}

		int[] boundNodeNumbSorted=new int[L];
		for(int i=0;i<L;i++)
			boundNodeNumbSorted[i]= boundNodeNumb[indSorted[i+1]];


		return boundNodeNumbSorted;
	}


	public void setSliceBounds(Model model){

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			if(model.region[ir].rotor){

				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				{
					int[] vertNumb=model.element[i].getVertNumb();
					for(int j=0;j<model.nElVert;j++)
					{
						model.node[vertNumb[j]].rotor=true;
					}
				}
			}

		double tmax=-10,tmin=10,rmin=1e4,rmax=0,rm=0,h,hmin=1e10,hmax=0;
		double p2=2*PI;
		for(int i=1;i<=model.numberOfNodes;i++)
		{
			Vect z=model.node[i].getCoord();

			double s=new Vect(z.el[0],z.el[1]).norm();
			double t;
			if(s==0) t=0;
			else t=util.getAng(z);
			if(Math.abs(t-p2)<epsAng) t=0;

			if(t>tmax) tmax=t;
			if(t<tmin) tmin=t;
			if(s>rmax) rmax=s;
			if(s<rmin) rmin=s;


			if(model.dim==3) {
				h=z.el[2];
				if(h>hmax) hmax=h;
				if(h<hmin) hmin=h;
			}

			if(model.motor){
				if(model.node[i].rotor) 
					if(s>rm) rm=s;
			}

		}


		model.r1=rmin;
		model.r2=rmax;
		model.h1=hmin;
		model.h2=hmax;


		model.rm=rm;
	

		model.alpha1=tmin;
		model.alpha2=tmax;
		model.spaceBoundary=new double[2*model.dim];
		model.spaceBoundary[0]=model.r1;
		model.spaceBoundary[1]=model.r2;
		model.spaceBoundary[2]=model.alpha1;
		model.spaceBoundary[3]=model.alpha2;

		if(model.dim==3){
			model.spaceBoundary[4]=model.h1;
			model.spaceBoundary[5]=model.h2;
		}

		if(tmax-tmin>3.5) model.fullMotor=true;
		


	}

	public void mapCommonNodes(Model model){

		if(!model.rotateConnect) return;

		model.commonNodes=getCommonNodeSorted(model);
		if( !model.hasTwoNodeNumb) return;

		for(int iz=0;iz<model.commonNodes.length;iz++)
		for(int i=0;i<model.commonNodes[iz][0].length;i++){
			int n1=model.commonNodes[iz][0][i];
			int n2=model.commonNodes[iz][1][i];
			model.node[n1].setMap(n2);

		}

	}
	
	public int[][][] getCommonNodeSorted(Model model){
	

		int nNodes=model.numberOfNodes;
		
		Vect hh=new Vect(model.numberOfNodes);
		
		for(int i=1;i<=nNodes;i++){
			if(model.dim==3)
			hh.el[i-1]=model.node[i].getCoord(2);
		}
		
		int[] index=hh.bubble();

		int nZ=1;
		for(int i=1;i<hh.length;i++){
			if(hh.el[i]-hh.el[i-1]>epsr) nZ++;
		}
		int[][] allNodeNumbsLayered=new int[nZ][nNodes/nZ];

		int ix=0;

		for(int iz=0;iz<nZ;iz++)
			for(int j=0;j<allNodeNumbsLayered[0].length;j++){
				allNodeNumbsLayered[iz][j]=index[ix++]+1;
			}	

		int[][][] commonNodeNumbSorted=new int[nZ][2][];


		boolean rotOut=model.region[model.numberOfRegions].rotor;
				
		double rm=0;
		
		if(rotOut){
			 rm=1e10;
		}
	
	
		for(int iz=0;iz<nZ;iz++)
			for(int j=0;j<allNodeNumbsLayered[iz].length;j++){
				int i=allNodeNumbsLayered[iz][j];
	
			Vect z=model.node[i].getCoord();
			double s=new Vect(z.el[0],z.el[1]).norm();
			if(model.node[i].rotor){
				if(rotOut){
					if(s<rm) rm=s;
				}
				else 
					if(s>rm) rm=s;
			}
				

		}
		model.rm=rm;

			
		for(int iz=0;iz<nZ;iz++)	
		{
			 ix=0;

				for(int j=0;j<allNodeNumbsLayered[iz].length;j++){
					int i=allNodeNumbsLayered[iz][j];

			if(Math.abs(model.node[i].getCoord().v2().norm()-rm)<epsr) 
			{
				ix++;
				model.node[i].common=true;
			}

		}
				int L=ix;
				if(model.hasTwoNodeNumb)
					L=ix/2;
			if(L==0) continue;
		
		commonNodeNumbSorted[iz][0]=new int[L];
		commonNodeNumbSorted[iz][1]=new int[L];

		int[] commonNodeNumb1=new int[L];
		int[] commonNodeInd1=new int[nNodes+1];
		int[] commonNodeNumb2=new int[L];
		int[] commonNodeInd2=new int[nNodes+1];
		Vect angs1=new Vect(L);
		Vect angs2=new Vect(L);
		int ixr=0;
		int ixs=0;
		
		for(int j=0;j<allNodeNumbsLayered[0].length;j++){
		int i=allNodeNumbsLayered[iz][j];
			if(model.node[i].common) {
				if(model.node[i].rotor){

					commonNodeInd1[i]=ixr;
					commonNodeNumb1[ixr]=i;
					angs1.el[ixr]=util.getAng(model.node[i].getCoord().v2());
					ixr++;
				}
				else{

					commonNodeInd2[i]=ixs;
					commonNodeNumb2[ixs]=i;
					angs2.el[ixs]=util.getAng(model.node[i].getCoord().v2());
					ixs++;
				}
			}
		}
			

		int[] indSorted1=angs1.bubble();
		int[] indSorted2=angs2.bubble();

		
		
		for(int i=0;i<L;i++){
			commonNodeNumbSorted[iz][0][i]= commonNodeNumb1[indSorted1[i]];
			commonNodeNumbSorted[iz][1][i]= commonNodeNumb2[indSorted2[i]];
		}
		
			}

		return commonNodeNumbSorted;
	}


	public int[][] mapBorderNodes2D(Model model,int nb1,int nb2){
		
		int nNodes=model.numberOfNodes;

		int[] onb1=new int[nNodes/2];
		int[] onb2=new int[nNodes/2];

		int ix1=0,ix2=0;
		for(int i=1;i<=nNodes;i++){
			
		
			if(model.node[i].onBound[nb1]) {
				onb1[ix1]=i;
				ix1++;

			}
			 if(model.node[i].onBound[nb2]) {
				onb2[ix2]=i;
				ix2++;

			}
		}

		if(ix1!=ix2)
			System.err.println("Mapping boundary nodes failed. " +
					"Boundaries have unequal node numbers;"+ ix1+" and "+ix2);

		Vect v1=new Vect(ix1);
		Vect v2=new Vect(ix1);

		for(int i=0;i<ix1;i++)
		{
			v1.el[i]=model.node[onb1[i]].getR();
			v2.el[i]=model.node[onb2[i]].getR();
		}

		int[] inds1=v1.bubble();
		int[] inds2=v2.bubble();

		int[][] map=new int[ix1][2];
		for(int i=0;i<ix1;i++){
			map[i][0]=onb1[inds1[i]];
			map[i][1]=onb2[inds2[i]];
		}


		return map;

	}


	public int[][] mapBorderNodes3D(Model model,int nb1,int nb2){
		double eps=1e-6;

		int nNodes=model.numberOfNodes;

		int[] onb1=new int[nNodes/2];
		int[] onb2=new int[nNodes/2];

		int ix1=0,ix2=0;
		for(int i=1;i<=nNodes;i++){

			if(model.node[i].onBound[nb1]) {
				onb1[ix1]=i;
				ix1++;

			}
			 if(model.node[i].onBound[nb2]) {
				onb2[ix2]=i;
				ix2++;

			}
		}
		


		if(ix1!=ix2)
			System.err.println("Mapping boundary nodes failed. " +
					"Boundaries have unequal node numbers;"+ ix1+" and "+ix2);



		double theta=model.alpha2-model.alpha1;

		Mat R=util.rotEuler(new Vect(0,0,1),theta);


		int[][] map=new int[ix1][2];

		for(int i=0;i<ix1;i++){
			map[i][0]=onb1[i];
			Vect v1=model.node[onb1[i]].getCoord();
			for(int j=0;j<ix1;j++){
				Vect v2=model.node[onb2[j]].getCoord();
				double d=R.mul(v1).sub(v2).norm();
			//	util.pr(d);
				if(d<eps) {
					map[i][1]=onb2[j];
					break;
				}
			}
		}




		return map;

	}

	public void setNodeOnBound(Model model)
	{
	
		double eps=1e-6;
		if( model.coordCode==1) {setNodeOnMotorBound(model);return;}

		boolean[][] nodeOnBound=new boolean[model.numberOfNodes+1][model.nBoundary];


		for(int i=1;i<=model.numberOfNodes;i++){

			for(int p=0;p<model.dim;p++){
				if(abs(model.node[i].getCoord(p)-model.spaceBoundary[2*p])<eps){
					nodeOnBound[i][2*p]=true;
				

				}
				
				if(abs(model.node[i].getCoord(p)-model.spaceBoundary[2*p+1])<eps){
					nodeOnBound[i][2*p+1]=true;
				}
			}
	}
		



		for(int i=1;i<=model.numberOfNodes;i++){
			model.node[i].onBound=nodeOnBound[i];
		}



	}


	public void setNodeOnMotorBound(Model model){

		boolean[][] nodeOnBound=new boolean[model.numberOfNodes+1][model.nBoundary];
		double r1=model.r1;
		double r2=model.r2;
		double a1=model.alpha1;
		double a2=model.alpha2;
		double h1=model.h1;
		double h2=model.h2;

	

		if(r2-r1==0){
			a1=-1000;
			a2=1000;
			r1=0;
			r2=model.getRmax();
		}

		
		double epsr=this.epsr*model.scaleFactor,p2=PI*2;

		
	
			
		boolean closed=(a2-a1>6.2);

		for(int i=1;i<=1*model.numberOfNodes;i++){	

			Vect v=model.node[i].getCoord();

			double r=v.v2().norm();
			double alpha=0;
			if(!closed) {
				alpha=util.getAng(v.v2());
				if(abs(alpha-p2)<epsAng) alpha=0;
			}

			if(r<epsr) 
			{


				nodeOnBound[i][0]=true;
				if(!closed){
					nodeOnBound[i][2]=true;
				
					nodeOnBound[i][3]=true;
				}
			} 
			else {
				if(abs(r-r1)<epsr) {
			
					nodeOnBound[i][0]=true;

				}

				if (abs(r-r2)<epsr)
					nodeOnBound[i][1]=true;

				if(!closed){
					if(abs(alpha-a1)<epsAng) 
						nodeOnBound[i][2]=true;

					if(abs(alpha-a2)<epsAng) 
						nodeOnBound[i][3]=true;
				}
			}

			if(model.dim==3){

				if( v.el[2]<h1+epsr) nodeOnBound[i][4]=true;
				else if( v.el[2]>h2-epsr) nodeOnBound[i][5]=true;
			}
		}
		


		for(int i=1;i<=model.numberOfNodes;i++){
			model.node[i].onBound=nodeOnBound[i];
		}


	}

	private void detectFarboundaryEdges(Model model){

		if(model.elCode<2)  return;
		
		int nEdge=model.numberOfEdges;

		int ns=6;

		if(model.elCode==2) ns=4;
		if(model.elCode==3) ns=6;
		else if(model.elCode==4) ns=6;
		else if(model.elCode==5) ns=5;
		
		IntVect[] edgeElement=new IntVect[nEdge+1];
		for(int i=1;i<=nEdge;i++)
			edgeElement[i]=new IntVect(ns);


		boolean[] edgeCounted=new boolean[nEdge+1];
		int[] indx=new int[nEdge+1];
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){

		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
			
			int[] edgeNumb=model.element[i].getEdgeNumb();

			for(int j=0;j<model.nElEdge;j++){
				int ne=edgeNumb[j];
				if(indx[ne]==edgeElement[ne].length-1)
					edgeElement[ne].extend(ns);

				edgeElement[ne].el[(indx[ne]++)]=i;
				edgeCounted[ne]=true;
			}			
		}
		}

		
		int nx=0;
		int k;
		boolean[] onSurf=new boolean[nEdge+1];

		for(int i=1;i<=nEdge;i++){
			if(!edgeCounted[i]) continue;
			k=0;
			for(int j=0;j<indx[i];j++)
				if(edgeElement[i].el[j]>0)
					k++;
			int[] nn=new int[2*indx[i]];
			int jx=0;
			int n1=model.edge[i].node[0].id;
			for(int j=0;j<indx[i];j++){
				int[] edgeNumb=model.element[edgeElement[i].el[j]].getEdgeNumb();
				for(int p=0;p<model.nElEdge;p++){
					int ep=edgeNumb[p];
					if(ep==i) continue;
					if(model.edge[ep].node[0].id==n1) {
					nn[jx++]=model.edge[ep].node[1].id;
					}
					else{
						if(model.edge[ep].node[1].id==n1) 
							nn[jx++]=model.edge[ep].node[0].id;
											}
				}

			}
			Arrays.sort(nn);
			int q=0;
			for(int p=1;p<nn.length;p++)
				if(nn[p]!=nn[p-1]) q++;

		if(q==indx[i])
			
			{
			
				onSurf[i]=true;
				nx++;
			}

		}
		
		for(int i=1;i<=nEdge;i++){
			if(onSurf[i]) model.edge[i].setKnownA(0);

		}
		util.pr("Number of edges on far boundary = "+nx);


}	


}
