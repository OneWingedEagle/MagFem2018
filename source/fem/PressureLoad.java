package fem;

import static java.lang.Math.PI;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import io.Loader;

import math.SpMatAsym;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan. Created Aug 20, 2012.
 */

public class PressureLoad {
	
	public class Segment{
		
		int[] nodeIds;
		double length;
		
		private Segment(int n){
			nodeIds=new int[n];
		}
	}

	public Segment[]  segments;
	Vect normal;
	
	double pressure;
	int  press_reg;

	public PressureLoad(){}
	

	public void readPressureLoad(Loader loader, BufferedReader br, Model model) throws IOException {

		model.setEdge();
		

		double minEdgeLenghth = model.minEdgeLength;
		util.pr(" minEdgeLength ------------------------- " + minEdgeLenghth);
		double clearFact = 1e-4;
		
		String line;
		
		double extention_fact = 0.01;

			line = loader.getNextDataLine(br," /* PRESSURE */");
			pressure = 1e6*loader.getScalarData(line);
		
			line = loader.getNextDataLine(br," /* nx  ny nz */");
			
			double[] nrm=loader.getCSV(line);
			
			if(model.dim==2) normal=new Vect(nrm[0],nrm[1]);
			else new Vect(nrm[0],nrm[1],nrm[2]);
			
			normal.normalize();
			
			line = loader.getNextDataLine(br," / * data type */");

			int type =  loader.getIntData(line);


			if (type == 0) {
				line = br.readLine();

				int nm = loader.getIntData(line);

		
				segments = new Segment[nm];


				for (int k = 0; k < nm; k++) {
					line = br.readLine();

					int[] nn = loader.getCSInt(line);
				//	if (model.dim == 2) {
			
						Node node1=model.node[nn[0]];
						Node node2=model.node[nn[1]];
						segments[k]=new Segment(2);
						segments[k].nodeIds[0]=nn[0];
						segments[k].nodeIds[1]=nn[1];
						segments[k].length= node1.getCoord().sub(node2.getCoord()).norm();

				}
			} else if( type==1 || type==2) {
				if (model.dim == 3) {
					
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

						segments = new Segment[nm];

						for (int k = 0; k < nm; k++) {

							segments[k] = new Segment(4);
							for(int p=0;p<4;p++){
								segments[k].nodeIds[p]=tempEd[k][p];
							}
						}
							
						
						
				} else {


					line = loader.getNextDataLine(br," / * Reg ID */");

					int nreg = loader.getIntData(line);
					
					press_reg = nreg;
					
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

					segments = new Segment[nm];

					for (int k = 0; k < nm; k++) {

						Node node11 = model.node[tempEd[k][0]];
						Node node12 = model.node[tempEd[k][1]];

						segments[k]=new Segment(2);
						segments[k].nodeIds[0]=tempEd[k][0];
						segments[k].nodeIds[1]=tempEd[k][1];
				
						segments[k].length= node11.getCoord().sub(node12.getCoord()).norm();

					}
				}

			}

		

	}
	
	public void setPressure(Model model){
		
		for (int k = 0; k < segments.length; k++) {

			int[] nids=segments[k].nodeIds;
			
			Node node1 = model.node[nids[0]];
			Node node2 = model.node[nids[1]];
	

			Vect v1 = node1.getCoord();
			Vect v2 = node2.getCoord();
			Vect v12 = v2.sub(v1);


			v12 = v2.sub(v1);
			double edgeLength = v12.norm();
			
			Vect F=normal.times(0.5*pressure*edgeLength);
			
			if(node1.F==null) node1.F=F.deepCopy();
			else node1.F=node1.F.add(F);
			
			if(node2.F==null) node2.F=F.deepCopy();
			else node2.F=node1.F.add(F);
			

		//	Vect edgeDir = v12.times(1. / edgeLength);


		}

	}
		


}

