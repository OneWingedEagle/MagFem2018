package meshFactory;

import fem.Model;
import fem.ModelGeo;
import math.Vect;

public class Rubiks extends MeshManipulator{
	
	public static void main(String[] args){
		
		Rubiks rb =new Rubiks();

		rb.makeRubix();
	}
	
	
	public void makeRubix(){


		double[][] cube0={{-50,50,-50,50,-50,50}};
		
		

			double scale=1;
		//	scale=1000;

			for(int j=0;j<cube0.length;j++)
				for(int k=0;k<cube0[0].length;k++){
					cube0[j][k]*=scale;
				
				}
			
			Geometry mg0=new Geometry(cube0);


			
							mg0.blockName[0]="cube"+(1);




			
	for(int j=0;j<cube0.length;j++){
				
				for(int k=0;k<cube0[0].length;k++){

					mg0.minMeshRight[j][k]=2;
					mg0.minMeshLeft[j][k]=2;
					mg0.baseLeft[j][k]=2000;
					mg0.baseRight[j][k]=2000;
		
			
				}	
	}


	MeshGeneration mgen=new MeshGeneration();
	
	

			Model model0=mgen.getOrthogMesh(mg0,"",false);
			
			Model model1=new Model(1,27,27*8,"hexahedron");
			
	
			model1.region[1]=model0.region[1].deepCopy();

			int nx=0;
			for(int i=1;i<=27;i++){

				int[] vn=model0.element[i].getVertNumb();
				for(int j=0;j<8;j++){
					nx++;
					model1.node[nx].setCoord(model0.node[vn[j]].getCoord());
					model1.element[i].setVertNumb(j,nx);
				
					
				}
				
			}
			
			model1.scaleFactor=1;
			
			Model[] models=new Model[27];
			
			
			double d=.12;
			int ix=0;
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
					for(int k=0;k<3;k++){
						Vect trans=new Vect((1-i)*d,(1-j)*d,(1-k)*d);
					
						models[ix]=model1.deepCopy();
						this.translate(models[ix],trans,false);
						
						ix++;
						}
			
			Model model=models[0].deepCopy();
			model.region[1].setName("reg"+1);
			for(int j=1;j<27;j++){
				model=this.assemble(model, models[j], false);
				model.region[j+1].setName("reg"+(j+1));
			}


				String bun=System.getProperty("user.dir") + "\\rubix.txt";
				model.writeMesh(bun);	
				
				
				model.setFemCalc();
				
				
				
				for(int ir=1;ir<=model.numberOfRegions;ir++){
					int jx=0;
					
					int[] nn=model.getRegNodes(ir);
					
					Vect v=new Vect(3);
					for(int j=0;j<nn.length;j++)
						v=v.add(model.node[nn[j]].getCoord());
					v=v.times(1.0/nn.length);
				
					
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
						
						int[] vn=model.element[i].getVertNumb();
						
						Vect c=model.getElementCenter(i);
		
					double eps=.001;
					double Ti=0;
						if(model.getElementVolume(i)>1e-5)
						{
							if(c.el[0]<v.el[0]-eps) Ti=1;
							else if(c.el[0]>v.el[0]+eps) Ti=2;
							else if(c.el[1]<v.el[1]-eps) Ti=3;
							else if(c.el[1]>v.el[1]+eps) Ti=4;
							else if(c.el[2]<v.el[2]-eps) Ti=5;
							else if(c.el[2]>v.el[2]+eps) Ti=6;
						}
						
				
							
						for(int j=0;j<8;j++){
						
								model.node[vn[j]].T=Ti;
						}
					}
				}
					String color=System.getProperty("user.dir") + "\\rubixColor.txt";
				model.writeNodalScalar(color);				
				


	}

}
