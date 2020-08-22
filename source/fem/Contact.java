package fem;

import static java.lang.Math.PI;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import io.Loader;

import math.SpMatAsym;

import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan. Created Aug 20, 2012.
 */

public class Contact {
	
	public class MasterEntity{
		
		int[] nodeIds;
		double length;
		
		private MasterEntity(int n){
			nodeIds=new int[n];
		}
	}



	public SpMatAsym[] node_node_mat;
	public SpMatAsym[] constraint_matrix_N;
	public SpMatAsym[] constraint_matrix_N_trp;
	public SpMatAsym[] constraint_matrix_T;
	public SpMatAsym[] constraint_matrix_T_trp;
	//private SpMatAsym Gcft = null;
	public SpMatAsym[] constraint_matrix_STK;
	public SpMatAsym[] constraint_matrix_STK_trp;


	public Vect lamN[], gap[];
	public Vect lamT[], slide[];// ,slide_prev;


	public Vect [] weights;
	public int numContacts;
	public double[] penFactor;
	//public double[] fn_ratio;
	public double[] master_edge_size;
	public double[] fric_coef;
	public int[] type;
	public Node[][] slaveNodes;
	public int[] slaveReg;
	public int[] masterReg;

	public MasterEntity[][]  master_entities;
	public boolean stick[][], landed_stick[][];

	public boolean[][] contacting;

	public Element[][] masterElems;

	public int[][] normalIndex;

	//int[][] u_index;
	public int[] numContactingNodes;
	public int totalnumContactingNodes;


	public Vect[][] normals;
	public  Vect[][] tangentials;
	public boolean frictional=false;


	public Contact(int numCont){



		numContacts = numCont;
		

		slaveNodes = new Node[numCont][];
		
		slaveReg=new int[numCont];
		masterReg=new int[numCont];
		
		master_entities = new MasterEntity[numCont][];

		penFactor = new double[numCont];
		fric_coef = new double[numCont];

		master_edge_size = new double[numCont];
		
		weights = new Vect[numCont];
		
		normalIndex=new int[numCont][];

		gap = new Vect[numCont];
		slide = new Vect[numCont];
		lamN = new Vect[numCont];
		lamT = new Vect[numCont];
		
		normals = new Vect[numContacts][];
		tangentials = new Vect[numContacts][];
		
		type=new int[numContacts];
		
		node_node_mat=new SpMatAsym[numContacts] ;
		constraint_matrix_N=new SpMatAsym[numContacts] ;
		constraint_matrix_N_trp=new SpMatAsym[numContacts] ;
		constraint_matrix_T=new SpMatAsym[numContacts] ;
		constraint_matrix_T_trp=new SpMatAsym[numContacts] ;
		constraint_matrix_STK=new SpMatAsym[numContacts] ;
		constraint_matrix_STK_trp=new SpMatAsym[numContacts] ;


		stick=new  boolean[numContacts][];
		landed_stick=new  boolean[numContacts][];
		landed_stick=new  boolean[numContacts][];
		contacting=new  boolean[numContacts][];
		
		numContactingNodes=new int[numContacts];
	}
	

	public void readContacts(Loader loader, BufferedReader br, Model model) throws IOException {

		model.setEdge();
		

		double minEdgeLenghth = model.minEdgeLength;
		util.pr(" minEdgeLength ------------------------- " + minEdgeLenghth);
		double clearFact = 1e-4;
		
		String line;
		
		double extention_fact = 0.01;

		for (int i = 0; i < numContacts; i++) {
			line = loader.getNextDataLine(br," /* PENALTY FACTOR */");
			penFactor[i] = loader.getScalarData(line);
			line = loader.getNextDataLine(br," /* COEF OF FRICTION */");
			fric_coef[i] = loader.getScalarData(line);
			line = loader.getNextDataLine(br," / * data type */");

			int type = 0;
			if (!line.contains("slav"))
				type = loader.getIntData(line);

			if (type == 0) {
				line = loader.getNextDataLine(br," /* NUM SLAVE NODES */");

				int ns = loader.getIntData(line);

				slaveNodes[i] = new Node[ns];

				for (int k = 0; k < ns; k++) {
					line = br.readLine();
					int sn = loader.getIntData(line);
					slaveNodes[i][k] = model.node[sn];
				}
			} else if (type == 1 || type==2){
				Vect v1 = null;
				Vect v2 = null;
				Vect v3 = null;
				Vect v4 = null;
				
				if (model.dim == 3) {
					//util.pr("This format of setting contact not ready yet.");
					line = loader.getNextDataLine(br," /* REGION ID */");
					int nreg = loader.getIntData(line);
					
					if(type==1){
					line = loader.getNextDataLine(br," /* CORENR 1 */");

					int n1 = loader.getIntData(line);

					line = loader.getNextDataLine(br," /* CORENR 2 */");

					int n2 = loader.getIntData(line);
					
					line = loader.getNextDataLine(br," /* CORENR 3 */");

					int n3 = loader.getIntData(line);
					
					line = loader.getNextDataLine(br,"/* CORENR 4 */");

					int n4 = loader.getIntData(line);

					Node node1 = model.node[n1];
					Node node2 = model.node[n2];
					Node node3 = model.node[n3];
					Node node4 = model.node[n4];
					
			
			

					 v1 = node1.getCoord();
					 v2 = node2.getCoord();
					 v3 = node3.getCoord();
					 v4 = node4.getCoord();
					}
					else{
						line = loader.getNextDataLine(br," /* v1x  v1y v1z */");

						v1=new Vect(loader.getCSV(line));

						line = loader.getNextDataLine(br," /* v2x  v2y v2z */");

						v2=new Vect(loader.getCSV(line));
						
						line = loader.getNextDataLine(br," /* v3x  v3y v3z */");

						v3=new Vect(loader.getCSV(line));

						line = loader.getNextDataLine(br," /* v4x  v4y v4z */");

						v4=new Vect(loader.getCSV(line));
						
					}
		
					Vect v13=v3.sub(v1);
					Vect v24=v4.sub(v2);
					
					v1 = v1.add(v13.times(-extention_fact));
					v3 = v3.add(v13.times(extention_fact));
					v2 = v2.add(v24.times(-extention_fact));
					v4 = v4.add(v24.times(extention_fact));
					
					Vect v12=v2.sub(v1);
					Vect v23=v3.sub(v2);
					Vect v34=v4.sub(v3);
					Vect v41=v1.sub(v4);

						
					int ns = 0;
					int[] nnr = model.getRegNodes(nreg);
					int[] temp = new int[nnr.length];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						Vect v = model.node[n].getCoord();
						Vect v1v = v.sub(v1);
						Vect v2v = v.sub(v2);
						Vect v3v = v.sub(v3);
						Vect v4v = v.sub(v4);
						
				
						
						Vect cross[]=new Vect[4];
						
						cross[0]=v12.cross(v1v);
						cross[1]=v23.cross(v2v);
						cross[2]=v34.cross(v3v);
						cross[3]=v41.cross(v4v);
						
						
						boolean infront=true;
						for(int j=0;j<4;j++){
							for(int m=j+1;m<4;m++){
								double dot=cross[j].dot(cross[m]);
								if(dot<0) {
									infront=false;
									break;
								}
								}
							if(!infront) break;
							}
						
						if(!infront) continue;
						
						Vect normal=v13.cross(v23);
						double dist=Math.abs(v1v.dot(normal));
						
						if(dist>clearFact * minEdgeLenghth) continue;


								temp[ns] = n;
								ns++;
							
						

					}
					
					util.pr("contact "+i+":  num. slave nodes ------------------------- " + ns);
					
					slaveNodes[i] = new Node[ns];

					for (int k = 0; k < ns; k++) {

						int sn = temp[k];
						slaveNodes[i][k] = model.node[sn];
						
						//util.pr(sn);
					}
				//	util.pr(ns);ss


				} else {
					line = loader.getNextDataLine(br," /* REGION ID */");

					int nreg = loader.getIntData(line);
					slaveReg[i] = nreg;
				
					
					line = loader.getNextDataLine(br," / * n1 */");

					int n1 = loader.getIntData(line);

					line = loader.getNextDataLine(br," / * n2 */");
					int n2 = loader.getIntData(line);

					Node node1 = model.node[n1];
					Node node2 = model.node[n2];

					Vect v11 = node1.getCoord();
					Vect v22 = node2.getCoord();
				//	v11.hshow();
				//	v22.hshow();
					Vect edgeDir = v22.sub(v11).normalized();
					int ns = 0;

					int[] nnr = model.getRegNodes(nreg);
					int[] temp = new int[nnr.length];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						Vect v = model.node[n].getCoord();
						Vect vv1 = v.sub(v11);
						Vect vv2 = v.sub(v22);
						double dot = vv1.dot(vv2);

						if (dot <= 0) {

							double proj = vv1.dot(edgeDir);
							double vv1n = vv1.norm();
							if (Math.abs(proj - vv1n) < 1e-4 * minEdgeLenghth) {
								// util.pr(dot);

								temp[ns] = n;
								ns++;
							}
						}

					}
					slaveNodes[i] = new Node[ns];

					for (int k = 0; k < ns; k++) {

						int sn = temp[k];
						slaveNodes[i][k] = model.node[sn];
					}

				}
			}

			line = loader.getNextDataLine(br," / * data type */");
			type = 0;
			if (!line.contains("mast"))
				type = loader.getIntData(line);
			if (type == 0) {
				line = br.readLine();

				int nm = loader.getIntData(line);

		
				master_entities[i] = new MasterEntity[nm];


				for (int k = 0; k < nm; k++) {
					line = br.readLine();

					int[] nn = loader.getCSInt(line);
				//	if (model.dim == 2) {
			
						Node node1=model.node[nn[0]];
						Node node2=model.node[nn[1]];
						master_entities[i][k]=new MasterEntity(2);
						master_entities[i][k].nodeIds[0]=nn[0];
						master_entities[i][k].nodeIds[1]=nn[1];
						master_entities[i][k].length= node1.getCoord().sub(node2.getCoord()).norm();
	/*				} else {

						
						if (nn.length == 3) {
							masterFacets[i][k] = new Element("triangle");
							masterFacets[i][k].setVertNumb(nn);

						} else if (nn.length == 4) {
							masterFacets[i][k] = new Element("quad");
							masterFacets[i][k].setVertNumb(nn);

						}
					}*/
				}
			} else if( type==1 || type==2) {
				if (model.dim == 3) {
				//	util.pr("This format of setting contact not ready yet.");
					
					byte[][] arr_hxa = { { 0, 1,2,3 }, { 4,7,6,5 }, { 0,4,5,1 },{2,6,7,3},{0,3,7,4 }, {1,5,6,2}};
					//byte[][] arr_penta = { { 0, 2,1}, { 3,4,5 }, { 0,1,4,3 }, { 1,2,5,4 },{ 2,0,3,5 } };
					byte[][] arr_penta = { { 0, 2,1}, { 3,4,5 }, { 0,3,4,1 }, { 1,4,5,2 },{ 2,5,3,0 } };


					byte[][] edgeLocalNodes = null;
					if (model.elCode == 4) {
						edgeLocalNodes = arr_hxa;
					//	;
					} else if (model.elCode == 3) {
						edgeLocalNodes = arr_penta;
					}
						
						Vect v1 = null;
						Vect v2 = null;
						Vect v3 = null;
						Vect v4 = null;
						
							//util.pr("This format of setting contact not ready yet.");
							line = loader.getNextDataLine(br," /* REGION ID */");
							int nreg = loader.getIntData(line);
							
							if(type==1){
							line = loader.getNextDataLine(br," /* CORENR 1 */");

							int n1 = loader.getIntData(line);

							line = loader.getNextDataLine(br," /* CORENR 2 */");

							int n2 = loader.getIntData(line);
							
							line = loader.getNextDataLine(br," /* CORENR 3 */");

							int n3 = loader.getIntData(line);
							
							line = loader.getNextDataLine(br,"/* CORENR 4 */");

							int n4 = loader.getIntData(line);

							Node node1 = model.node[n1];
							Node node2 = model.node[n2];
							Node node3 = model.node[n3];
							Node node4 = model.node[n4];
							
					
					

							 v1 = node1.getCoord();
							 v2 = node2.getCoord();
							 v3 = node3.getCoord();
							 v4 = node4.getCoord();
							}
							else{
								line = loader.getNextDataLine(br," /* v1x  v1y v1z */");

								v1=new Vect(loader.getCSV(line));

								line = loader.getNextDataLine(br," /* v2x  v2y v2z */");

								v2=new Vect(loader.getCSV(line));
								
								line = loader.getNextDataLine(br," /* v3x  v3y v3z */");

								v3=new Vect(loader.getCSV(line));

								line = loader.getNextDataLine(br," /* v4x  v4y v4z */");

								v4=new Vect(loader.getCSV(line));
								
							}
				

						
						Vect v13=v3.sub(v1);
						Vect v24=v4.sub(v2);
						
						v1 = v1.add(v13.times(-extention_fact));
						v3 = v3.add(v13.times(extention_fact));
						v2 = v2.add(v24.times(-extention_fact));
						v4 = v4.add(v24.times(extention_fact));
						
						Vect v12=v2.sub(v1);
						Vect v23=v3.sub(v2);
						Vect v34=v4.sub(v3);
						Vect v41=v1.sub(v4);

						

						int[] nnr = model.getRegNodes(nreg);
						boolean[] nc = new boolean[model.numberOfNodes + 1];
						for (int k = 0; k < nnr.length; k++) {
							int n = nnr[k];
							nc[n] = true;
						}
						
						int nm=0;
						
						int[][] tempEd = new int[model.numberOfEdges * edgeLocalNodes.length][4];
						for (int ie = model.region[nreg].getFirstEl(); ie <= model.region[nreg].getLastEl(); ie++) {
							int[] vertNumb = model.element[ie].getVertNumb();
							
							
							
				

							for (int j = 0; j < edgeLocalNodes.length; j++) {
								
								int nv=edgeLocalNodes[j].length;
								Vect cent=new Vect(3);
								for(int k=0;k<nv;k++){
									int nx = vertNumb[edgeLocalNodes[j][k]];
									Vect v = model.node[nx].getCoord();
									cent=cent.add(v);
								}
								
								cent.timesVoid(1./nv);

								Vect v = cent;
								Vect v1v = v.sub(v1);
								Vect v2v = v.sub(v2);
								Vect v3v = v.sub(v3);
								Vect v4v = v.sub(v4);
									
									Vect cross[]=new Vect[4];
									
									cross[0]=v12.cross(v1v);
									cross[1]=v23.cross(v2v);
									cross[2]=v34.cross(v3v);
									cross[3]=v41.cross(v4v);
									
									
									boolean infront=true;
									for(int k=0;k<4;k++){
										for(int m=k+1;m<4;m++){
											double dot=cross[k].dot(cross[m]);
											if(dot<0) {
												infront=false;
												break;
											}
											}
										if(!infront) break;
										}
									
									if(!infront) continue;
									
									Vect normal=v13.cross(v23);
									double dist=Math.abs(v1v.dot(normal));
									
									if(dist>1e-4  * minEdgeLenghth) continue;
					
								for(int p=0;p<4;p++){
									int n12 = vertNumb[edgeLocalNodes[j][p]];
						
										tempEd[nm][p] = n12;
									
									}
							//	util.hshow(tempEd[nm]);
								nm++;
							}

						}

						util.pr("contact "+i+":  num. master facets ------------------------- " + nm);
						master_entities[i] = new MasterEntity[nm];

						for (int k = 0; k < nm; k++) {

							master_entities[i][k] = new MasterEntity(4);
							for(int p=0;p<4;p++){
								master_entities[i][k].nodeIds[p]=tempEd[k][p];
							}
						}
							
						
						
				} else {


					line = loader.getNextDataLine(br," / * Reg ID */");

					int nreg = loader.getIntData(line);
					
					masterReg[i] = nreg;
					
					line = loader.getNextDataLine(br," / * n1 */");
					int n1 = loader.getIntData(line);
					line = loader.getNextDataLine(br," / * n2 */");
					int n2 = loader.getIntData(line);
					Node node1 = model.node[n1];
					Node node2 = model.node[n2];
					Vect v1 = node1.getCoord();
					Vect v2 = node2.getCoord();
					Vect edgeDir = v2.sub(v1).normalized();
					int nm = 0;

					int[] nnr = model.getRegNodes(nreg);
					boolean[] nc = new boolean[model.numberOfNodes + 1];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						nc[n] = true;
					}

					byte[][] arr0 = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
					byte[][] arr1 = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };

					byte[][] edgeLocalNodes = null;
					if (model.elCode == 0) {
						edgeLocalNodes = arr0;
						
					} else if (model.elCode == 1) {
						edgeLocalNodes = arr1;
					}

					int[][] tempEd = new int[model.numberOfEdges * edgeLocalNodes.length][2];
					for (int ie = model.region[nreg].getFirstEl(); ie <= model.region[nreg].getLastEl(); ie++) {
						int[] vertNumb = model.element[ie].getVertNumb();

						for (int j = 0; j < edgeLocalNodes.length; j++) {
							int n11 = vertNumb[edgeLocalNodes[j][0]];
							int n12 = vertNumb[edgeLocalNodes[j][1]];

							Vect v = model.node[n11].getCoord().add(model.node[n12].getCoord()).times(0.5);
							Vect vv1 = v.sub(v1);
							Vect vv2 = v.sub(v2);
							double dot = vv1.dot(vv2);

							if (dot <= 0) {

								double proj = vv1.dot(edgeDir);
								double vv1n = vv1.norm();
								if (Math.abs(proj - vv1n) < clearFact * minEdgeLenghth) {
									// util.pr(dot);

									tempEd[nm][0] = n11;
									tempEd[nm][1] = n12;
									nm++;
								}
							}
						}

					}

					master_entities[i] = new MasterEntity[nm];

					for (int k = 0; k < nm; k++) {

						Node node11 = model.node[tempEd[k][0]];
						Node node12 = model.node[tempEd[k][1]];

						master_entities[i][k]=new MasterEntity(2);
						master_entities[i][k].nodeIds[0]=tempEd[k][0];
						master_entities[i][k].nodeIds[1]=tempEd[k][1];
				
						master_entities[i][k].length= node11.getCoord().sub(node12.getCoord()).norm();

					}
				}

			}

		}
	
	boolean node_duplic=false;
	
	if(node_duplic){
		
		for (int contId = 0; contId < numContacts; contId++) {
			Node[] sns = slaveNodes[contId];
			MasterEntity[]  med = master_entities[contId];
			
			boolean[] coincid=new boolean[model.numberOfNodes+1];

			for (int i = 0; i < coincid.length; i++) 
				coincid[i]=false;
			
			for (int i = 0; i < sns.length; i++) {
				int sn = slaveNodes[contId][i].id;
				for (int j = 0; j < med.length; j++) {
					int n1=med[j].nodeIds[0];
					int n2=med[j].nodeIds[1];
				///	util.pr(sn+" "+n1+"  "+n2);
					if(sn==n1){
						coincid[sn]=true;
						break;
					}
					else if(sn==n2){
						coincid[sn]=true;	
						break;
					}
				
				}
		}
			
			int[] map=new int[model.numberOfNodes+1];
			for (int i = 0; i < map.length; i++)
				map[i]=0;
			
			int extera_nodes=0;
			for (int i = 0; i < coincid.length; i++) {
				if(coincid[i]) {
					extera_nodes++;
					map[i]=model.numberOfNodes+extera_nodes;
				//	util.pr(i+" ---- "+map[i]);
				}
			}

			Vect[] dup_coord=new Vect[extera_nodes];
			int ix=0;

			for (int i = 0; i < map.length; i++) {
				if(map[i]>0) dup_coord[ix++]=model.node[i].getCoord();
			}
			
			int nRegions=model.numberOfRegions;
			int nElements=model.numberOfElements;
			int nNodes=model.numberOfNodes+extera_nodes;
			Model md1=new Model(nRegions,nElements,nNodes,model.elType);
	
			int n1=1,N;
			for(int i=1;i<=md1.numberOfRegions;i++){

				md1.region[i].setFirstEl(model.region[i].getFirstEl());
				md1.region[i].setLastEl(model.region[i].getLastEl());
				md1.region[i].setName(model.region[i].getName());
				md1.region[i].setMaterial(model.region[i].getMaterial());


			}

			for(int i=1;i<=model.numberOfNodes;i++)	
				md1.node[i].setCoord(model.node[i].getCoord());
			
			for(int i=1;i<=extera_nodes;i++)	{
				md1.node[i+model.numberOfNodes].setCoord(dup_coord[i-1]);
			}

			for(int ir=1;ir<=md1.numberOfRegions;ir++)		{

		
				boolean on_slave=slaveReg[contId]==ir;
				
				for(int i=md1.region[ir].getFirstEl();i<=md1.region[ir].getLastEl();i++)	{
					int[] vn=model.element[i].getVertNumb();
						
					for(int k=0;k<md1.nElVert;k++){
						if(map[vn[k]]==0 )
							md1.element[i].setVertNumb(k,vn[k]);
						else{
							if(on_slave)
								md1.element[i].setVertNumb(k,map[vn[k]]);
							else
								md1.element[i].setVertNumb(k,vn[k]);

						}
					}

				
			}
			}



			md1.scaleFactor=model.scaleFactor;
			
			String folder=new File(model.meshFilePath).getParentFile().getPath();
			md1.meshFilePath=model.meshFilePath;
			
			model=md1.deepCopy();
			model.meshFilePath=md1.meshFilePath;
			
			
			for (int i = 0; i <slaveNodes[contId].length; i++) {
				int sn = slaveNodes[contId][i].id;
			
				if(map[sn]>0){ slaveNodes[contId][i]=model.node[map[sn]];
				
				}

		}
			String file = folder + "//duplicated.txt";
			
			for (int i = 1; i <= model.numberOfElements; i++) {
				
				int[] vn=model.element[i].getVertNumb();
				for(int k=0;k<vn.length;k++){
						Vect v=model.node[vn[k]].getCoord();

						v.el[0]+=.01*i;
						//v.el[1]+=.01*i;
						model.node[vn[k]].setCoord(v);
				}
				
			}
			

			model.writeMesh(file);
		

	}
		
		
		
	}

	
	for (int contId = 0; contId < numContacts; contId++) {
		Node[] sns = slaveNodes[contId];
		MasterEntity[]  med = master_entities[contId];
		util.pr("===== Slave Nodes ====  "+sns.length);
		util.pr("===== Master Edges/Facets ====  "+med.length);
		
	}
	}



	public void initialize(Model model) {

		int[][] ne = new int[model.numberOfNodes + 1][20];
		int[] nz = new int[model.numberOfNodes + 1];

		for (int i = 1; i <= model.numberOfElements; i++) {
			int[] vn = model.element[i].getVertNumb();
			for (int j = 0; j < vn.length; j++) {
				ne[vn[j]][nz[vn[j]]++] = i;

			}
		}

		masterElems = new Element[numContacts][];
		for (int contId = 0; contId < numContacts; contId++) {
			
			double mu = fric_coef[contId];
			if(mu!=0) frictional=true;
			
			if (model.dim == 2) {
				masterElems[contId] = new Element[master_entities[contId].length];

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids=master_entities[contId][k].nodeIds;
					
					Node node1 = model.node[nids[0]];
					Node node2 = model.node[nids[1]];

					int ie = 0;
					for (int j = 0; j < nz[node1.id]; j++) {
						for (int p = 0; p < nz[node2.id]; p++) {

							if (ne[node1.id][j] == ne[node2.id][p]) {
								ie = ne[node1.id][j];
								break;
							}
							if (ie > 0)
								break;
						}
					}

					if (ie > 0)
						masterElems[contId][k] = model.element[ie];
					else
						util.pr("master edge ( " + node1.id + ", " + node2.id + " ) belongs to no element.");

				}
			} else {

				masterElems[contId] = new Element[master_entities[contId].length];

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids = master_entities[contId][k].nodeIds;

					int ie = 0;
					for (int j = 0; j < nz[nids[0]]; j++) {
						for (int p = 0; p < nz[nids[1]]; p++) {
							for (int q = 0; q < nz[nids[2]]; q++) {

								if (ne[nids[0]][j] == ne[nids[1]][p] && ne[nids[0]][j] == ne[nids[2]][q]) {
									ie = ne[nids[0]][j];

									break;

								}
								if (ie > 0)
									break;
							}
							if (ie > 0)
								break;
						}
					}

					if (ie > 0) {

						masterElems[contId][k] = model.element[ie];
					} else
						util.pr("master facet ( " + nids[0] + ", " + nids[1] + ", " + nids[2] + ", " + nids[3]
								+ " ) belongs to no element.");

				}
			}

			// for(int k=0;k<masterEdges[contId].length;k++)
			// util.hshow(masterElems[contId][k].getVertNumb());
		}






		for (int contId = 0; contId < numContacts; contId++) {
			
			landed_stick[contId] = new boolean[model.numberOfNodes + 1];

			stick[contId] = new boolean[model.numberOfNodes + 1];
			contacting[contId] = new boolean[model.numberOfNodes + 1];
			
			for (int i = 1; i <= model.numberOfNodes; i++) {
				stick[contId][i] = true;
				landed_stick[contId][i] = true;
			}

			
			
			lamN[contId] = new Vect(model.numberOfNodes + 1);
			lamT[contId] = new Vect(model.numberOfNodes + 1);


			weights[contId] = new Vect(model.numberOfNodes + 1);// .ones(model.Ks.nRow);


			gap[contId] = new Vect(model.numberOfNodes + 1);
			slide[contId] = new Vect(model.numberOfNodes + 1);
			
			type[contId]=0;
			
			int numSn = slaveNodes[contId].length;
			int numMed = 0;
	
			numMed = master_entities[contId].length;
			normalIndex[contId] = new int[numSn];
			normals[contId] = new Vect[numMed];
			tangentials[contId] = new Vect[numMed];
			
			int nnSize = model.numberOfNodes + 1;
			node_node_mat[contId] = new SpMatAsym(nnSize, nnSize);

		}

	
		

	}
	


}

