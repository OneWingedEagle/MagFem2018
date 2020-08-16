package fem;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.BitSet;

import math.Mat;
import math.Vect;
import math.util;

public class Calculator {
	double[][] PW;
	double[][] PWNL;
	double[][] PW3ang;
	double[][] PWtetra;
	boolean centerNu;
	public int dim,elCode,nElVert,nElEdge,numberOfElements,numberOfNodes,numberOfRegions;
	public int struc2D; // 0: plane stress. 1: plane strain

	public Calculator(){	}

	public Calculator(Model model)
	{
		this.nElVert=model.nElVert;
		this.nElEdge=model.nElEdge;
		this.numberOfElements=model.numberOfElements;
		this.numberOfNodes=model.numberOfNodes;
		this.numberOfRegions=model.numberOfRegions;
		this.dim=model.dim;
		this.elCode=model.elCode;
		this.PW=gaussInteg(2);
		this.PWNL=gaussInteg(3);
		this.PW3ang=gaussInteg3(4);
		this.PWtetra=gaussIntegTetra(4);

		this.centerNu=false;
		
		this.struc2D=model.struc2D;
	}


	public Vect getElementB(Model model,int i, Vect[] rotNe){

		Edge[] elEdge=model.elementEdges(i);
		Vect B=new Vect(model.dim);
		for(int j=0;j<model.nElEdge;j++)		{
			B=B.add(rotNe[j].times(elEdge[j].A));
		}

		return B;

	}

	public double getElementAQ(Model model,int i){
		Edge[] elEdge=model.elementEdges(i);
		double[] Ae=new double[4];
		for(int j=0;j<4;j++)
			Ae[j]=elEdge[j].A;
		double A=0;
		Vect zero=new Vect(2);

		double[] Ne=NeQuad(zero);

		for(int j=0;j<model.nElEdge;j++)	{		
			A= A+Ne[j]*Ae[j];
		}
		return  A;	

	}


	public Mat nodalMat(Vect nu,Vect gN1, Vect gN2){
		Vect gNnuN2=nu.times(gN2);
		double p=gN1.dot(gNnuN2);
		Mat M=new Mat(this.dim,this.dim);
		for(int i=0;i<this.dim;i++)
			for(int j=0;j<this.dim;j++)
				if(i==j)
					M.el[i][j]=p-gN1.el[j]*gNnuN2.el[j];
				else
					M.el[i][j]=-gN1.el[j]*gNnuN2.el[j];
		return M;

	}

	public double[][] He(Model model,int nBH,int nLam,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){
		if(model.coupled) return HeCoupled(model,nBH,nLam,ie,nonLinear,eddy,hasJ,hasM);
		if(this.elCode==0 && !model.axiSym) return He3ang(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==0&& model.axiSym) return He3angAxi(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==1&& !model.axiSym) return HeQuad(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==1&& model.axiSym) return HeQuadAxi(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==3) return HePrism(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		//else if(this.elCode==2) return He3ang2nd(model,nBH,ie,nonLinear,eddy,hasJ,hasM);

		Vect nu=new Vect(this.dim),nuVar=new Vect(this.dim);
		Vect B=new Vect();
		Vect  M=new Vect();		


		if(hasM){ M=model.element[ie].getM();}
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		int n;

		if(!nonLinear){
			nu=model.element[ie].getNu();
			n=this.PW[0].length; 
		}
		else{
			n=this.PWNL[0].length; 
			if(this.centerNu){
				double nux,nuv=0;
				B=model.element[ie].getB();
				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				nux=model.BH[nBH].getNu(Bn);
				nu=new Vect(nux,nux,nux);
				if(Bn>0){
					nuv=model.BH[nBH].getdHdB(Bn);	
					nuVar.el[0]=(nuv-nux)/Bn2;
					nuVar.el[1]=nuVar.el[0];
					nuVar.el[2]=nuVar.el[0];

				}
			}
		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];


		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];

		Vect sigma=new Vect();

		if(eddy)
			sigma=model.element[ie].getSigma();

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){


					if(nonLinear){
						localCo.el[0]=this.PWNL[0][p];
						localCo.el[1]=this.PWNL[0][q];
						localCo.el[2]=this.PWNL[0][r];

						if(n!=2)
							ws=this.PWNL[1][p]*this.PWNL[1][q]*this.PWNL[1][r];
						else
							ws=1;
						B=this.getBAt(model,ie,localCo);
						double Bn2=B.norm2();
						double Bn=sqrt(Bn2);

						if(!this.centerNu){
							double nux,nuv=0;


							nux=model.BH[nBH].getNu(Bn);
							nu=new Vect(nux,nux,nux);
							if(Bn>0){
								nuv=model.BH[nBH].getdHdB(Bn);	
								nuVar.el[0]=(nuv-nux)/Bn2;
								nuVar.el[1]=nuVar.el[0];
								nuVar.el[2]=nuVar.el[0];

							}
						}

					}
					else
					{
						localCo.el[0]=this.PW[0][p];
						localCo.el[1]=this.PW[0][q];
						localCo.el[2]=this.PW[0][r];

						if(n!=2)
							ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r];
						else
							ws=1;


					}

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());
				//	if(ie==1) jac.show("%10.5e");
					wsJ=ws*detJac;


					rotNe=rotNe(jac,localCo,edgeDir);
					//if(ie==1) rotNe[0].hshow();
					if(hasJ || eddy)
						Ne=Ne(jac,localCo,edgeDir);

					for(int i=0;i<model.nElEdge;i++){

						if(hasJ){

							Cj[i]=Cj[i].add(Ne[i].times(wsJ));
						}

						if(hasM){
							C[i]+=wsJ*rotNe[i].dot(M);

						}

						for(int j=0;j<=i;j++){

							term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
							H1[i][j]+=term1;
						
			
							if(nonLinear){
								term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
								H2[i][j]+=term2;


							}

							if(eddy){
								term3=wsJ*Ne[i].dot(Ne[j].times(sigma));
								H3[i][j]+=term3;

							}


						}
					}
				}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy){
			lowSym(H3);
		}
		
		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;

/*if(ie==1) {
	int[] nn=model.element[ie].getVertNumb();
	for (int j=0;j<nn.length;j++)
		model.node[nn[j]].getCoord().hshow();
}*/
		//if(ie==1) util.show(H1);

		return H1;


	}


	public Vect elemVectCLN(Model model,int ie, double[] loss){
		
		if(this.elCode==1) return elemVectCLNQuad(model, ie,loss);

		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;

		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();
		
		double sigma=model.element[ie].getSigma().el[0];

		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];

		Mat G;
		Vect Je=new Vect(this.dim);
		Vect[] localGradsN=new Vect[model.nElVert];;
		Vect[] gradN=new Vect[model.nElVert];
		Vect[] gradPhi=new Vect[model.nElVert];;
	

		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	{
			A[j]=elemEdges[j].getA();
		}
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);

					G=jac.inv3();

					localGradsN=localGradN(localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=Ne(jac,localCo,edgeDir);

					for(int j=0;j<model.nElVert;j++){
						gradN[j]=G.mul(localGradsN[j]);
						gradPhi[j]= gradN[j].times(nodePhi[j]);
					}

					
					Je=Je.times(0);
					for(int j=0;j<model.nElEdge;j++){
						Je=Je.add(Ne[j].times(A[j]));
					
					}

					for(int j=0;j<model.nElVert;j++){

						Je=Je.add(gradPhi[j]);
					}

					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i].dot(Je)*wsJ;

						}
					
					heat+=Je.dot(Je)*wsJ;
	
				}
				
		elemV=elemV.times(sigma);

		
		loss[0]=heat*sigma;

		return elemV;


	}

	

	public Vect elemVectCLNQuad(Model model,int ie, double[] loss){
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;

		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect[] localGradsN=new Vect[model.nElVert];;

		double sigma=model.element[ie].getSigma().el[0];

		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		double[] Ne=new double[model.nElEdge];

		Mat G;
		double Jez=0;



		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}
		
		double jezphi=0;
		for(int j=0;j<model.nElVert;j++){
			jezphi+=nodePhi[j]/model.height;
		}

		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	{
			A[j]=elemEdges[j].getA();
		}
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){
			

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);
					
					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=NeQuad(localCo);
				
					Jez=0;
					for(int j=0;j<model.nElEdge;j++){
						Jez+=Ne[j]*A[j];
					
					}
					
					Jez+=jezphi/(n*n);

			

					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i]*Jez*wsJ;

						}
					
					heat+=Jez*Jez*wsJ;

				}
				
		elemV=elemV.times(sigma);

		
		loss[0]=heat*sigma;
		return elemV;


	}

	public double[][] HeCoupled(Model model,int nBH,int nLam,int ie, boolean nonLinear,boolean  eddy,boolean hasJ,boolean hasM){
		if(this.elCode==0) return HeCoupled3ang(model,nBH,nLam,ie,nonLinear,eddy, hasJ,hasM);
		if(this.elCode==1) return HeCoupledQuad(model,nBH,nLam,ie,nonLinear,eddy, hasJ,hasM);


		boolean MS=model.cpvms && model.element[ie].hasMS() ;		

		double E=model.element[ie].getYng().el[0];
		double v=model.element[ie].getPois().el[0];

		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		double G=E/(1+v);


		Vect nu=new Vect(3),nuVar=new Vect(3);
		Vect B=new Vect();
		Vect M=new Vect();		

		if(hasM) M=model.element[ie].getM();

		double  sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];


		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);


		double detJac,ws=1,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];
		Vect[] Ne=new Vect[model.nElEdge];

		Vect localCo=new Vect(3);
		int n;
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					if(nonLinear){

						localCo.el[0]=this.PWNL[0][p];
						localCo.el[1]=this.PWNL[0][q];
						localCo.el[2]=this.PWNL[0][r];

						B=this.getBAt(model, ie,localCo);
						double Bn2=B.norm2();
						double Bn=sqrt(Bn2);

						double sn=model.getStressB(ie, B);

						double nux=model.BHS[nBH].getNu(Bn,sn);
						nu=new Vect(this.dim).add(nux);

						if(Bn>0){
							double nuVx=(model.BHS[nBH].getdHdB(Bn,sn)-nux)/Bn2;
							nuVar=new Vect(this.dim).add(nuVx);
						}


						if(MS &&  Bn>1e-1){

							double lamb=model.lamB[nLam].getLam(Bn);
							double dLambdB=model.lamB[nLam].getdLamdB(Bn);
							double d2LambdB2=model.lamB[nLam].getd2LamdB2(Bn);
							double kms=0,kvar=0;

							kms=G*lamb*dLambdB/Bn;
							kvar=G/Bn*(pow(dLambdB,2)+lamb*d2LambdB2)-kms/Bn;
							kvar=kvar/(2*Bn);
							nu=nu.add(kms);
							nuVar=nuVar.add(kvar);

						}

					}

					else{
						localCo.el[0]=this.PW[0][p];
						localCo.el[1]=this.PW[0][q];
						localCo.el[2]=this.PW[0][r];
						nu=model.element[ie].getNu();
					}


					jac=jacobian(vertexNode, localCo);

					detJac=abs(jac.determinant());

					if(n!=2){
						if(nonLinear){
							ws=this.PWNL[1][p]*this.PWNL[1][q]*this.PWNL[1][r]*detJac;
						}
						else
							ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					}
					else
						ws=detJac;


					rotNe=rotNe(jac,localCo,edgeDir);
					if(hasJ || eddy)
						Ne=Ne(jac,localCo,edgeDir);

					for(int i=0;i<model.nElEdge;i++){

						if(hasJ)
							//C[i]+=ws*Ne[i].dot(J);
							Cj[i]=Cj[i].add(Ne[i].times(ws));

						if(hasM)
							C[i]+=ws*rotNe[i].dot(M);

						for(int j=0;j<=i;j++){

							term1=ws*rotNe[i].times(nu).dot(rotNe[j]);
							H1[i][j]+=term1;

							if(nonLinear){
								term2=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
								H2[i][j]+=term2;
							}

							if(eddy){
								term3=ws*Ne[i].dot(Ne[j].times(sigma));
								H3[i][j]+=term3;
							}


						}
					}
				}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;

		return H1;

	}

	public double[][] HeQuad(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){


		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		

		boolean[] edgeDir=model.element[ie].getEdgeReverse();



		if(hasM){ M=model.element[ie].getM();}


		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double sigma=0;
		if(eddy){
			sigma=model.element[ie].getSigma().el[2];			


		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;

		Vect[] rotNe=new Vect[model.nElEdge];
		double [] Ne=new double[model.nElEdge];

		int n=0;
		Vect localCo=new Vect(2);
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){

				if(nonLinear){
					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];

					if(n!=2)
						ws=this.PWNL[1][p]*this.PWNL[1][q];
					else
						ws=1;
					B=this.getBAt(model,ie,localCo);
					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					if(!this.centerNu){

						double nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux);
						if(Bn>0){
							double nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];


						}
					}


				}
				else
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;

				rotNe=rotNe(jac,localCo,edgeDir);
				if(hasJ || eddy)
					Ne=NeQuad(localCo);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]+=wsJ*Ne[i];
					}

					if(hasM)
						C[i]+=wsJ*rotNe[i].dot(M);

					for(int j=0;j<=i;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;


						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;

						}

						if(eddy){

							term3=wsJ*Ne[i]*Ne[j]*sigma;
							H3[i][j]+=term3;

						}
					}
				}
			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);


		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;



		return H1;

	}

	public double[][] HeQuadAxi(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){


		boolean centerNu=true;

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect  M=new Vect();		


		if(hasM){ M=model.element[ie].getM();}


		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double sigma=0;
		if(eddy){
			sigma=model.element[ie].getSigma().el[2];			


		}


		double rr=0;
		boolean rcent=true;

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];

		double[] rrj=new double[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++){
			double r1=vertexNode[i].getCoord(0)+1e-6;
			rrj[i]=1.0/r1;
		}


		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;

		Vect[] rotNe=new Vect[model.nElEdge];
		double [] Ne=new double[model.nElEdge];

		int n=0;
		Vect localCo=new Vect(2);
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){

				if(nonLinear){
					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];

					if(n!=2)
						ws=this.PWNL[1][p]*this.PWNL[1][q];
					else
						ws=1;
					B=this.getBAt(model,ie,localCo);

					double Bn2=B.norm2();

					double Bn=sqrt(Bn2);

					if(!centerNu){

						double nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux);
						if(Bn>0){
							double nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];


						}
					}


				}
				else
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;

				rotNe=gradN(jac,localCo);

				Ne=NQuad(localCo);


				rr=0;
				if(!rcent){

					for(int i=0;i<model.nElEdge;i++)
						rr+=Ne[i]*rrj[i];
				}
				else{

					for(int i=0;i<model.nElEdge;i++)
						rr+=.25*rrj[i];
				}


				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]+=wsJ*Ne[i];
					}

					if(hasM)
						C[i]+=wsJ*rotNe[i].dot(M);




					for(int j=0;j<4;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;


						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;
						}

						if(eddy){
							term3=wsJ*Ne[i]*Ne[j]*sigma;
							H3[i][j]+=term3;

						}
					}
				}
			}


		for(int i=0;i<model.nElEdge;i++)
			for(int j=0;j<model.nElEdge;j++){
				H1[i][j]*=rr;
				if(nonLinear){
					H2[i][j]*=rr;
				}
				if(eddy)
					H3[i][j]*=rr;

			}

		/*	lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);*/

		//util.show(H2);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;



		return H1;

	}

	public double[][] HePrism(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){

		//return new double[1][1];

		Vect nu=new Vect(this.dim),nuVar=new Vect(this.dim);
		Vect B=new Vect();
		Vect  M=new Vect();		

		if(hasM){ M=model.element[ie].getM();}
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		int n3ang=this.PW3ang.length; 		

		int nGauss;

		if(!nonLinear){
			nu=model.element[ie].getNu();
			nGauss=this.PW[0].length; 
		}
		else{
			nGauss=this.PWNL[0].length; 
			if(this.centerNu){
				double nux,nuv=0;
				B=model.element[ie].getB();
				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				nux=model.BH[nBH].getNu(Bn);
				nu=new Vect(nux,nux,nux);
				if(Bn>0){
					nuv=model.BH[nBH].getdHdB(Bn);	
					nuVar.el[0]=(nuv-nux)/Bn2;
					nuVar.el[1]=nuVar.el[0];
					nuVar.el[2]=nuVar.el[0];

				}
			}
		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];


		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];




		///


		Vect sigma=new Vect();

		if(eddy)
			sigma=model.element[ie].getSigma();

		for(int p=0;p<n3ang;p++)
			for(int q=0;q<nGauss;q++)
			{
				if(nonLinear){
					localCo.el[0]=this.PW3ang[p][0];
					localCo.el[1]=this.PW3ang[p][1];
					localCo.el[2]=this.PWNL[0][q];

					if(nGauss!=2)
						ws=this.PWNL[1][q]*this.PW3ang[p][2];
					else
						ws=this.PW3ang[p][2];

					B=this.getBAt(model,ie,localCo);
					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					if(!this.centerNu){
						double nux,nuv=0;


						nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux,nux);

						if(Bn>0){
							nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];
							nuVar.el[2]=nuVar.el[0];

						}
					}

				}
				else
				{
					localCo.el[0]=this.PW3ang[p][0];
					localCo.el[1]=this.PW3ang[p][1];
					localCo.el[2]=this.PW[0][q];

					if(nGauss!=2)
						ws=this.PW[1][p]*this.PW3ang[p][2];
					else
						ws=this.PW3ang[p][2];


				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;


				rotNe=rotNePrism(jac,localCo,edgeDir);
				if(hasJ || eddy)
					Ne=NePrism(jac,localCo,edgeDir);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]=Cj[i].add(Ne[i].times(wsJ));
					}

					if(hasM){
						C[i]+=wsJ*rotNe[i].dot(M);

					}

					for(int j=0;j<=i;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;

						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;


						}

						if(eddy){
							term3=wsJ*Ne[i].dot(Ne[j].times(sigma));
							H3[i][j]+=term3;

						}


					}
				}
			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy){
			lowSym(H3);
		}

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;


		return H1;


	}

	public double[][] HeCoupledQuad(Model model,int nBH,int nLam,int ie, boolean nonLinear,boolean  eddy,boolean hasJ,boolean hasM){

		boolean MS=!model.calCurve&&model.cpvms && model.element[ie].hasMS() ;	
		double E=model.element[ie].getYng().el[0];
		double v=model.element[ie].getPois().el[0];

		double G=E/(1+v);

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect  M=new Vect();		

		if(hasM) M=model.element[ie].getM();

		double  sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		
		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];

		double detJac,ws=1,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];
		double [] Ne=new double[model.nElEdge];

		Vect localCo=new Vect(2);
		int n;
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){

				if(nonLinear){

					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];

					B=this.getBAt(model, ie,localCo);

					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					double sn=model.getStressB(ie, B);
					double nux=model.BHS[nBH].getNu(Bn,sn);
					nu=new Vect(this.dim).add(nux);

					if(Bn>0){
						double nuVx=(model.BHS[nBH].getdHdB(Bn,sn)-nux)/Bn2;
						nuVar=new Vect(this.dim).add(nuVx);
					}

					if(MS &&  Bn>1e-2){

						double nums=0,numsVar=0;
						double[] numsv;
						numsv=getNumsv(model,nLam,ie,B,Bn,sn,G);					
						nums=numsv[0];
						numsVar=numsv[1];
						nu=nu.add(nums);
						nuVar=nuVar.add(numsVar);

					}
				}
				else{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					nu=model.element[ie].getNu();
				}




				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				if(n!=2){
					if(nonLinear){
						ws=this.PWNL[1][p]*this.PWNL[1][q]*detJac;
					}
					else
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
				}
				else
					ws=detJac;


				rotNe=rotNe(jac,localCo,edgeDir);

				if(hasJ || eddy)
					Ne=NeQuad(localCo);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ)
						Cj[i]+=ws*Ne[i];

					if(hasM)
						C[i]+=ws*rotNe[i].dot(M);

					for(int j=0;j<=i;j++){

						term1=ws*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;

						if(nonLinear){
							term2=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;


						}

						if(eddy){
							term3=ws*Ne[i]*Ne[j]*sigma;
							H3[i][j]+=term3;
						}


					}
				}
			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;

		return H1;


	}

	public double[] getNumsv(Model model,int nLam,int ie,Vect B,double Bn,double sn,double G){

		double[] numsv=new double[2];

		double Bn2=Bn*Bn;
		double dLambdB=model.lamBS[nLam].getdLamdB(Bn,sn);
		double d2LambdB2=model.lamBS[nLam].getd2LamdB2(Bn,sn);

		Vect B3=B.v3();

		Mat S=new Mat(3,3);
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				if(i==j) S.el[i][j]=(3*B3.el[i]*B3.el[j]-Bn2)/(2*Bn2);
				else
					S.el[i][j]=(3*B3.el[i]*B3.el[j])/(2*Bn2);

		Mat sig=model.element[ie].getStressTensor();
		double v=model.element[ie].getPois().el[0];
		//v=0;
		double sz=v*(sig.el[0][0]+sig.el[1][1]);


		Mat sig3=new Mat(3,3);
		for(int i=0;i<this.dim;i++)
			for(int j=0;j<this.dim;j++)
				sig3.el[i][j]=sig.el[i][j];

		sig3.el[2][2]=sz;

		sig3=sig3.times(1e6); // due to MPa to Pa conversion


		double C=0;
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				C+=S.el[i][j]*sig3.el[i][j];


		double dlamdsig=1e-6*model.lamBS[nLam].getdLamds(Bn, sn);

		double pp=G*dlamdsig;
		double nums0=(pp-1)*C*dLambdB/Bn;

		numsv[0]=nums0;	
		numsv[1]=nums0*(d2LambdB2/dLambdB-1.0/Bn);

		return numsv;




	}

	public double[][] He3ang(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);



		if(hasM){ M=model.element[ie].getM();}

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;

			nux=model.BH[nBH].getNu(Bn);

			nu=new Vect(nux,nux);
			if(Bn>0){
				nuv=model.BH[nBH].getdHdB(Bn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}
		}

		else
			nu=model.element[ie].getNu();

		/*		if(model.element[ie].getRegion()==8&& model.getElementCenter(ie).v2().norm()>-.07)
		{

			double  ss=0;
			if(model.element[ie].isDeformable())
				ss=model.element[ie].getStress().norm();
			nu.timesVoid(1+ss/50);
			nuVar.timesVoid(1+ss/50);

		}
		else if(model.element[ie].getRegion()==16){
			double  ss=0;
			if(model.element[ie].isDeformable())
				ss=model.element[ie].getStress().norm();
		//	nu.timesVoid(1-ss/130);
			//nuVar.timesVoid(1-ss/130);
			nu.timesVoid(1-.9);
			nuVar.timesVoid(1-.9);

		}*/

		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[3];
		double[] Cj=new double[3];

		double S=el3angArea(model,ie);
		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		double[] localNe=new double[3];
		Vect[] rotNe=rotNe3ang(model,ie);


		if(hasJ){
			for(int i=0;i<3;i++)
				Cj[i]=S/3;

		}

		if(hasM){
			for(int i=0;i<3;i++)
				C[i]=S*rotNe[i].dot(M);
		}

		for(int i=0;i<this.nElEdge;i++)
			for(int j=0;j<=i;j++)	{
				H1[i][j]=S*rotNe[i].times(nu).dot(rotNe[j]);

				if(nonLinear)
					H2[i][j]=S*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));

				if(eddy){
					int n=this.PW3ang.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PW3ang[p][0];
						lc.el[1]=this.PW3ang[p][1];
						localNe[0]=lc.el[0];
						localNe[1]=lc.el[1];
						localNe[2]=1-lc.el[0]-lc.el[1];
						jac=jacobian3ang(vertexNode,lc);

						detJac=abs(jac.determinant());
						ws=this.PW3ang[p][2]*detJac;
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}
				}


			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}

	public double[][] He3angAxi(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);



		if(hasM){ M=model.element[ie].getM();}

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;

			nux=model.BH[nBH].getNu(Bn);

			nu=new Vect(nux,nux);
			if(Bn>0){
				nuv=model.BH[nBH].getdHdB(Bn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}
		}

		else
			nu=model.element[ie].getNu();



		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];



		double rr=0;

		for(int i=0;i<model.nElEdge;i++){
			double r1=vertexNode[i].getCoord(0);
			rr+=0.33333333/r1;
		}



		double S=el3angArea(model,ie);
		double detJac,ws=1;
		Mat jac;
		double[] localNe=new double[3];
		Vect[] rotNe=rotNe3ang(model,ie);


		if(hasJ){
			for(int i=0;i<3;i++)
				Cj[i]=S/3;

		}

		if(hasM){
			for(int i=0;i<3;i++)
				C[i]=S*rotNe[i].dot(M);
		}

		for(int i=0;i<this.nElEdge;i++)
			for(int j=0;j<3;j++)	{
				H1[i][j]=S*rotNe[i].times(nu).dot(rotNe[j]);

				if(nonLinear)
					H2[i][j]=S*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));

				if(eddy){
					int n=this.PW3ang.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PW3ang[p][0];
						lc.el[1]=this.PW3ang[p][1];
						localNe[0]=lc.el[0];
						localNe[1]=lc.el[1];
						localNe[2]=1-lc.el[0]-lc.el[1];
						jac=jacobian3ang(vertexNode,lc);

						detJac=abs(jac.determinant());
						ws=this.PW3ang[p][2]*detJac;
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}
				}


			}

		for(int i=0;i<model.nElEdge;i++)
			for(int j=0;j<model.nElEdge;j++){
				H1[i][j]*=rr;
				if(nonLinear)
					H2[i][j]*=rr;
				if(eddy)
					H3[i][j]*=rr;

			}


		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}


	public double[][] HeCoupled3ang(Model model,int nBH,int nLam,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		boolean MS=!model.calCurve&&model.cpvms && model.element[ie].hasMS() ;	

		if(nBH==1) MS=false;
		double E=model.element[ie].getYng().el[0];
		double v=model.element[ie].getPois().el[0];

		double G=E/(1+v);
		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);

		if(hasM) M=model.element[ie].getM();

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;
			double sn=model.getStressB(ie,B);

			nux=model.BHS[nBH].getNu(Bn,sn);	
			nu=new Vect(nux,nux);
			if(Bn>0){
				nuv=model.BHS[nBH].getdHdB(Bn,sn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];

			}

			if(MS &&  Bn>1e-2){
				double nums=0,numsVar=0;
				double[] numsv=new double[2];
				numsv=getNumsv(model,nLam,ie,B,Bn,sn,G);	

				nums=numsv[0];
				numsVar=numsv[1];
				nu=nu.add(nums);
				nuVar=nuVar.add(numsVar);

			}

		}


		else
			nu=model.element[ie].getNu();



		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];	
		double[][] H3=new double[model.nElEdge][model.nElEdge];		
		double[] C=new double[3];
		double[] Cj=new double[3];

		double S=el3angArea(model,ie);
		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		double[] localNe=new double[3];
		Vect[] rotNe=rotNe3ang(model,ie);


		if(hasJ){
			for(int i=0;i<3;i++)
				Cj[i]=S/3;
		}

		if(hasM){
			for(int i=0;i<3;i++)
				C[i]=S*rotNe[i].dot(M);
		}
		for(int i=0;i<this.nElEdge;i++)
			for(int j=0;j<=i;j++)	{
				H1[i][j]=S*rotNe[i].times(nu).dot(rotNe[j]);
				if(nonLinear){
					H2[i][j]=S*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
				}

				if(eddy){
					int n=this.PW3ang.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PW3ang[p][0];
						lc.el[1]=this.PW3ang[p][1];
						localNe[0]=lc.el[0];
						localNe[1]=lc.el[1];
						localNe[2]=1-lc.el[0]-lc.el[1];
						jac=jacobian3ang(vertexNode,lc);

						detJac=abs(jac.determinant());
						ws=this.PW3ang[p][2]*detJac;
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}

				}
			}


		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;

		return H1;

	}

	public double[][] He3ang1st(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);

		if(hasM) M=model.element[ie].getM();

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[3];
		double[] Cj=new double[3];

		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		double[] localNe;
		Vect[] rotNe;


		int n=this.PW3ang.length; 					
		for(int p=0;p<n;p++){

			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];

			localNe=Ne3ang1st(lc);


			jac=jacobian3ang(vertexNode,lc);

			rotNe=rotNe3ang1st(jac,ie);

			detJac=abs(.5*jac.determinant());

			ws=this.PW3ang[p][2]*detJac;

			if(nonLinear){

				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				if(!this.centerNu){

					double nux=model.BH[nBH].getNu(Bn);
					nu=new Vect(nux,nux);
					if(Bn>0){
						double nuv=model.BH[nBH].getdHdB(Bn);	
						nuVar.el[0]=(nuv-nux)/Bn2;
						nuVar.el[1]=nuVar.el[0];


					}
				}

			}


			if(hasJ){
				for(int i=0;i<3;i++)
					Cj[i]+=ws*localNe[i];
			}

			if(hasM){
				for(int i=0;i<3;i++)
					C[i]+=ws*rotNe[i].dot(M);
			}

			for(int i=0;i<this.nElEdge;i++)
				for(int j=0;j<=i;j++)	{
					H1[i][j]+=ws*rotNe[i].times(nu).dot(rotNe[j]);
					if(nonLinear){
						H2[i][j]+=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
					}



					if(eddy){
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}

				}
		}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}

	/*	public double[][] HeTetra(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(3),nuVar=new Vect(3);
		Vect B=new Vect();
		Vect J=new Vect(), M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(3);

		if( hasJ)	{ J=model.element[ie].getJ().times(model.region[model.element[ie].getRegion()].getFactJ());}
		if(hasM) M=model.element[ie].getM();

		Vect sigma=new Vect(3);
		if(eddy)
			sigma=model.element[ie].getSigma().times(sigma);

		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
				nuVar.el[2]=nuVar.el[0];
			}

		}

		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];

		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		Vect[] localNe;
		Vect[] rotNe;


		int n=this.PWtetra.length; 					
		for(int p=0;p<n;p++){

			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];

			localNe=NeTetra(lc);

			jac=jacobianTetra(vertexNode,lc);

			rotNe=rotNeTetra(jac,ie);

			detJac=abs(.5*jac.determinant());

			ws=this.PWtetra[p][3]*detJac;

			if(nonLinear){

				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				if(!this.centerNu){

					double nux=model.BH[nBH].getNu(Bn);
					nu=new Vect(nux,nux);
					if(Bn>0){
						double nuv=model.BH[nBH].getdHdB(Bn);	
						nuVar.el[0]=(nuv-nux)/Bn2;
						nuVar.el[1]=nuVar.el[0];
						nuVar.el[2]=nuVar.el[0];

					}
				}

			}


			if(hasJ){
				for(int i=0;i<model.nElEdge;i++)
					C[i]+=localNe[i].dot(J)*ws;
			}

			if(hasM){
				for(int i=0;i<model.nElEdge;i++)
					C[i]+=ws*rotNe[i].dot(M);
			}

			for(int i=0;i<this.nElEdge;i++)
				for(int j=0;j<=i;j++)	{
					H1[i][j]+=ws*rotNe[i].times(nu).dot(rotNe[j]);
					if(nonLinear){
						H2[i][j]+=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
					}



					if(eddy){
						H3[i][j]+=ws*localNe[i].dot(localNe[j].times(sigma));	
					}

				}
		}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;


		return H1;
	}
	 */

	public Vect[] gradN(Mat jac,Vect localCo){

		if(this.elCode==3) return gradNPrism(jac,localCo);

		Mat invJac=jac.inv();
		Vect[] gradN=new Vect[this.nElVert];
		Vect[] localGradN=localGradN(localCo);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=invJac.mul(localGradN[i]);
		}

		return gradN;
	}

	public Vect[] gradNPrism(Mat jac,Vect localCo){


		Mat invJac=jac.inv();
		Vect[] gradN=new Vect[this.nElVert];
		Vect[] localGradN=localGradNPrism(localCo);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=invJac.mul(localGradN[i]);
		}

		return gradN;
	}

	

	public  Vect[] rotNe(Mat jac, Vect localCo, boolean[] edgeReverse){

		if(this.elCode==1) return rotNeQuad(jac,localCo);
		if(this.elCode==3) return rotNePrism(jac,localCo,edgeReverse);
		Vect[] rotNe=new Vect[this.nElEdge];

		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		rotNe[0]= grad[1].times(-(1-c)).add(grad[2].times(-(1-b))).times(0.125).cross(grad[0]); 
		rotNe[1]= grad[1].times((1-c)).add(grad[2].times(-(1+b))).times(0.125).cross(grad[0]);
		rotNe[2]= grad[1].times(-(1+c)).add(grad[2].times(+(1-b))).times(0.125).cross(grad[0]);
		rotNe[3]= grad[1].times((1+c)).add(grad[2].times(+(1+b))).times(0.125).cross(grad[0]);
		rotNe[4]= grad[0].times(-(1-c)).add(grad[2].times(-(1-a))).times(0.125).cross(grad[1]);
		rotNe[5]= grad[0].times((1-c)).add(grad[2].times(-(1+a))).times(0.125).cross(grad[1]);
		rotNe[6]= grad[0].times(-(1+c)).add(grad[2].times(+(1-a))).times(0.125).cross(grad[1]);
		rotNe[7]= grad[0].times((1+c)).add(grad[2].times(+(1+a))).times(0.125).cross(grad[1]);
		rotNe[8]= grad[0].times(-(1-b)).add(grad[1].times(-(1-a))).times(0.125).cross(grad[2]);
		rotNe[9]= grad[0].times((1-b)).add(grad[1].times(-(1+a))).times(0.125).cross(grad[2]);
		rotNe[10]= grad[0].times(-(1+b)).add(grad[1].times(+(1-a))).times(0.125).cross(grad[2]);
		rotNe[11]= grad[0].times((1+b)).add(grad[1].times(+(1+a))).times(0.125).cross(grad[2]);

	
		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}


	public double[] NQuad(Vect localCo){
		double a,b;
		a=localCo.el[0];b=localCo.el[1];;
		double[] N=new double[this.nElVert];
		double c=.25;
		N[0]=(1+a)*(1+b)*c;
		N[1]=(1-a)*(1+b)*c;
		N[2]=(1-a)*(1-b)*c;
		N[3]=(1+a)*(1-b)*c;

		return N;
	}
	public Vect[] rotNeQuad(Mat jac, Vect localCo){


		double a=localCo.el[0];
		double b=localCo.el[1];

		Vect[] rotNe=new Vect[this.nElVert];

		Mat invJac=jac.inv();

		Vect[] grad=new Vect[2];
		grad[0]=invJac.getColVect(0);
		grad[1]=invJac.getColVect(1);

		Vect v;
		v=grad[0].times(1+b).add(grad[1].times(1+a)).times(.25);
		rotNe[0]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(-(1+b)).add(grad[1].times(1-a)).times(.25);
		rotNe[1]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(b-1).add(grad[1].times(a-1)).times(.25);
		rotNe[2]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(1-b).add(grad[1].times(-(1+a))).times(.25);
		rotNe[3]= new Vect(v.el[1],-v.el[0]);

		return rotNe;
	}




	public Vect[] rotNePrism(Mat jac,Vect localCo,boolean[] edgeReverse){

		double c=localCo.el[2];


		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		Vect[] rotNe=new Vect[9];


		rotNe[0]=grad[0].cross(grad[1]).times(2*c);
		rotNe[1]=rotNe[0];
		rotNe[2]=rotNe[0];

		rotNe[3]=grad[0].cross(grad[1]).times(2*(1-c));
		rotNe[4]=rotNe[3];
		rotNe[5]=rotNe[3];

		rotNe[6]=grad[0].cross(grad[2]);
		rotNe[7]=grad[1].cross(grad[2]);
		rotNe[8]=grad[0].add(grad[1]).times(-1).cross(grad[2]);

		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}


	public double[] NeQuad(Vect localCo){

		double a=localCo.el[0];
		double b=localCo.el[1];

		double[] Ne=new double[this.nElVert];

		Ne[0]= (1+a)*(1+b)*0.25; 
		Ne[1]= (1-a)*(1+b)*0.25; 
		Ne[2]= (1-a)*(1-b)*0.25; 
		Ne[3]= (1+a)*(1-b)*0.25; 


		return Ne;
	}


	Vect[] localGradN(Vect localCo){
		if(this.elCode==0) return localGradN3ang();
		else if(this.elCode==1) return localGradNQuad(localCo);
		else if(this.elCode==3) return localGradNPrism(localCo);
		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Vect[] gradN=new Vect[this.nElVert];

		gradN[0]=new Vect((1+b)*(1+c),(1+a)*(1+c),(1+a)*(1+b));
		gradN[1]=new Vect(-(1+b)*(1+c),(1-a)*(1+c),(1-a)*(1+b));
		gradN[2]=new Vect(-(1-b)*(1+c),-(1-a)*(1+c),(1-a)*(1-b));
		gradN[3]=new Vect((1-b)*(1+c),-(1+a)*(1+c),(1+a)*(1-b));
		gradN[4]=new Vect((1+b)*(1-c),(1+a)*(1-c),-(1+a)*(1+b));
		gradN[5]=new Vect(-(1+b)*(1-c),(1-a)*(1-c),-(1-a)*(1+b));
		gradN[6]=new Vect(-(1-b)*(1-c),-(1-a)*(1-c),-(1-a)*(1-b));
		gradN[7]=new Vect((1-b)*(1-c),-(1+a)*(1-c),-(1+a)*(1-b));


		for(int i=0;i<this.nElVert;i++)
			gradN[i].timesVoid(0.125);
		return gradN;
	}

	public Vect[] localGradN3ang(){


		Vect[] gradN=new Vect[3];

		gradN[0]=new Vect(1,0,0);
		gradN[1]=new Vect(0,1,0);
		gradN[2]=new Vect(0,0,1);

		return gradN;
	}

	public Vect[] localGradNtet(){


		Vect[] gradN=new Vect[4];

		gradN[0]=new Vect(1,0,0);
		gradN[1]=new Vect(0,1,0);
		gradN[2]=new Vect(0,0,1);
		gradN[3]=new Vect(-1,-1,-1);

		return gradN;
	}



	public Vect[] localGradNPrism(Vect localCo){

		double a=localCo.el[0]/2;
		double b=localCo.el[1]/2;
		double c=.5-a-b;

		double h1=(1-localCo.el[2])/2;
		double h2=(1+localCo.el[2])/2;
		Vect[] gradN=new Vect[6];

		gradN[0]=new Vect(h2,0,a);
		gradN[1]=new Vect(0,h2,b);
		gradN[2]=new Vect(-h2,-h2,c);
		gradN[3]=new Vect(h1,0,-a);
		gradN[4]=new Vect(0,h1,-b);
		gradN[5]=new Vect(-h1,-h1,-c);

		return gradN;
	}

	public Vect[] localGradNQuad(Vect localCo){
		double a=localCo.el[0];
		double b=localCo.el[1];

		Vect[] gradN=new Vect[this.nElVert];

		gradN[0]=new Vect((1+b),(1+a));
		gradN[1]=new Vect(-(1+b),(1-a));
		gradN[2]=new Vect(-(1-b),-(1-a));
		gradN[3]=new Vect((1-b),-(1+a));


		for(int i=0;i<this.nElVert;i++)
			gradN[i].timesVoid(0.25);
		return gradN;
	}

	Vect[] gradN3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);

		Vect v1=vertexNode[0].getCoord();
		Vect v2=vertexNode[1].getCoord();
		Vect v3=vertexNode[2].getCoord();

		double rS=.5/el3angArea(model,ie);

		Vect[] gradN=new Vect[3];

		gradN[0]=new Vect(v2.el[1]-v3.el[1],v3.el[0]-v2.el[0]).times(rS);
		gradN[1]=new Vect(v3.el[1]-v1.el[1],v1.el[0]-v3.el[0]).times(rS);
		gradN[2]=new Vect(v1.el[1]-v2.el[1],v2.el[0]-v1.el[0]).times(rS);

		return gradN;
	}

	Vect[] gradNTet(Mat jac){

		Vect[] locgN=this.localGradNtet();
		Mat grads=new Mat(3,4) ;
		for(int i=0;i<4;i++)
			grads.setCol(locgN[i], i);

		grads=jac.inv().transp().mul(grads);


		Vect[] gradN=new Vect[4];

		for(int i=0;i<4;i++)
			gradN[i]=grads.getColVect(i);
		return gradN;
	}

	public double getElementArea(Model model,int i){
		if(this.elCode==0) return el3angArea(model,i);
		else 	if(this.elCode==1) return elQuadArea(model,i);
		else throw new NullPointerException(" Element is not 2D. ");

	}

	public double el3angArea(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=abs(v1.cross(v2).norm())/2;
		return S;
	}

	public double elPrismVol(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;
		double h=abs(vertexNode[0].getCoord(2)-vertexNode[3].getCoord(2));
		return S*h;
	}

	public double elPrismBase(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;

		return S;
	}

	public double elQuadArea(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[3].getCoord().sub(vertexNode[0].getCoord());
		Vect v3=vertexNode[1].getCoord().sub(vertexNode[2].getCoord());
		Vect v4=vertexNode[3].getCoord().sub(vertexNode[2].getCoord());
		double S=(v1.cross(v2).norm()+v4.cross(v3).norm())/2;
		return S;
	}

	double[] N(Vect localCo){
		if(elCode==0) return N3ang1st(localCo);
		else if(elCode==1) return NQuad(localCo);
		double a,b,c;
		a=localCo.el[0];b=localCo.el[1];c=localCo.el[2];
		double[] N=new double[this.nElVert];
		N[0]=(1+a)*(1+b)*(1+c)*0.125;
		N[1]=(1-a)*(1+b)*(1+c)*0.125;
		N[2]=(1-a)*(1-b)*(1+c)*0.125;
		N[3]=(1+a)*(1-b)*(1+c)*0.125;
		N[4]=(1+a)*(1+b)*(1-c)*0.125;
		N[5]=(1-a)*(1+b)*(1-c)*0.125;
		N[6]=(1-a)*(1-b)*(1-c)*0.125;
		N[7]=(1+a)*(1-b)*(1-c)*0.125;

		return N;
	}

	public  Vect[] Ne(Mat jac,Vect localCo, boolean[] edgeReverse){
		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];
		Vect[] Ne=new Vect[this.nElEdge];
		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		Ne[0]= grad[0].times((1-b)*(1-c)*0.125); 
		Ne[1]= grad[0].times((1+b)*(1-c)*0.125); 
		Ne[2]= grad[0].times((1-b)*(1+c)*0.125); 
		Ne[3]= grad[0].times((1+b)*(1+c)*0.125); 
		Ne[4]= grad[1].times((1-a)*(1-c)*0.125); 
		Ne[5]= grad[1].times((1+a)*(1-c)*0.125); 
		Ne[6]= grad[1].times((1-a)*(1+c)*0.125); 
		Ne[7]= grad[1].times((1+a)*(1+c)*0.125); 
		Ne[8]= grad[2].times((1-a)*(1-b)*0.125); 
		Ne[9]= grad[2].times((1+a)*(1-b)*0.125); 
		Ne[10]= grad[2].times((1-a)*(1+b)*0.125); 
		Ne[11]= grad[2].times((1+a)*(1+b)*0.125); 
		
		for(int k=0;k<Ne.length;k++)
			if(edgeReverse[k])
				Ne[k].timesVoid(-1);


		return Ne;

	}

	public  Vect[] NePrism(Mat jac,Vect localCo, boolean[] edgeReverse)
	{

		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Vect[] Ne=new Vect[this.nElEdge];

		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);


		Ne[0]= grad[1].times(a).sub( grad[0].times(b)).times(c);
		Ne[1]= Ne[0].sub( grad[1]).times(c);
		Ne[2]= Ne[0].add( grad[0]).times(c);


		Ne[3]= grad[1].times(a).sub( grad[0].times(b)).times(1-c);
		Ne[4]= Ne[3].sub( grad[1]).times(1-c);
		Ne[5]= Ne[3].add( grad[0]).times(1-c);

		Ne[6]= grad[2].times(a);
		Ne[7]= grad[2].times(b);
		Ne[8]= grad[2].times(1-a-b);
		
		for(int k=0;k<Ne.length;k++)
			if(edgeReverse[k])
				Ne[k].timesVoid(-1);

		return Ne;

	}

	double[] localNPrims(Vect lc){

		double[] N=new double[6];


		double c1=(1+lc.el[2])/2;
		double c2=(1-lc.el[2])/2;

		N[0]=lc.el[0]*c1;
		N[1]=lc.el[1]*c1;
		N[2]=(1-lc.el[0]-lc.el[1])*c1;
		N[3]=lc.el[0]*c2;
		N[4]=lc.el[1]*c2;
		N[5]=(1-lc.el[0]-lc.el[1])*c2;

		return N;


	}

	double[] NPrims(Node[] vertexNode, Vect globalCo){

		Vect P=new Vect(globalCo.el[0],globalCo.el[1]);
		Vect v3D1=vertexNode[3].getCoord();
		Vect v3D2=vertexNode[4].getCoord();
		Vect v3D3=vertexNode[5].getCoord();
		Vect v3D4=vertexNode[0].getCoord();

		double H=v3D1.el[2]-v3D4.el[2];
		double h=globalCo.el[2]-v3D4.el[2];

		Vect v1=new Vect(v3D1.el[0],v3D1.el[1]).sub(P);
		Vect v2=new Vect(v3D2.el[0],v3D2.el[1]).sub(P);
		Vect v3=new Vect(v3D3.el[0],v3D3.el[1]).sub(P);

		double c1=h/H;
		double c2=(1-h/H);
		double[] N=new double[6];

		double[] s=new double[3];
		s[0]=v2.cross(v3).norm();
		s[1]=v3.cross(v1).norm();
		s[2]=v1.cross(v2).norm();

		for(int i=0;i<6;i++)
			if(i<3)
				N[i]=s[i]*c1/2;
			else
				N[i]=s[i-3]*c2/2;

		return N;


	}



	//============

	public double[] N3ang(Node[] vertexNode, Vect globalCo){

		Vect v1=vertexNode[0].getCoord().sub(globalCo);
		Vect v2=vertexNode[1].getCoord().sub(globalCo);
		Vect v3=vertexNode[2].getCoord().sub(globalCo);

		double[] N=new double[3];

		N[0]=v2.cross(v3).norm();
		N[1]=v3.cross(v1).norm();
		N[2]=v1.cross(v2).norm();

		double S=N[0]+N[1]+N[2];

		for(int i=0;i<3;i++)
			N[i]/=S;

		return N;


	}

	public double[] N3ang(Model model,int ie, Vect globalCo){
		Node[] vertexNode=model.elementNodes(ie);

		return N3ang(vertexNode,globalCo);

	}

	public double[] N3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		Vect globalCo=vertexNode[0].getCoord().add(vertexNode[1].getCoord()).add(vertexNode[2].getCoord()).times(1.0/3);

		return N3ang(vertexNode,globalCo);


	}

	public double[] N3ang(Node[] vertexNode){

		Vect globalCo=vertexNode[0].getCoord().add(vertexNode[1].getCoord()).add(vertexNode[2].getCoord()).times(1.0/3);
		return N3ang(vertexNode,globalCo);


	}


	public double[] Ne3ang(Node[] vertexNode){
		return N3ang(vertexNode);

	}

	double[] Ne3ang(Node[] vertexNode, Vect globalCo){

		return N3ang(vertexNode, globalCo);
	}

	double[] N3ang1st(Vect localCo){
		double[] N=new double[3];
		N[0]=localCo.el[0];
		N[1]=localCo.el[1];
		N[2]=1-N[0]-N[1];
		return N;
	}

	public double [] Ntetra(Vect localCo){
		double[] N=new double[3];
		N[0]=localCo.el[0];
		N[1]=localCo.el[1];
		N[0]=localCo.el[2];
		N[2]=1-N[0]-N[1]-N[2];
		return N;
	}


	public double[] Ne3ang1st(Vect localCo){
		return N3ang1st(localCo);
	}


	double[] Ne3ang(Model model,int ie, Vect globalCo){
		return Ne3ang(model.elementNodes(ie), globalCo);
	}


	public Vect[] rotNe3ang(Node[] vertexNode){

		Vect v1=vertexNode[0].getCoord();
		Vect v2=vertexNode[1].getCoord();
		Vect v3=vertexNode[2].getCoord();

		double rS=1.0/v2.sub(v1).cross(v3.sub(v1)).norm();

		Vect[] rotNe=new Vect[3];

		rotNe[0]=new Vect(v3.el[0]-v2.el[0], v3.el[1]-v2.el[1]).times(rS);
		rotNe[1]=new Vect(v1.el[0]-v3.el[0], v1.el[1]-v3.el[1]).times(rS);
		rotNe[2]=new Vect(v2.el[0]-v1.el[0],v2.el[1]-v1.el[1]).times(rS);

		return rotNe;
	}

	private Vect[] rotNe3ang1st(Mat jac,int ie){


		Vect[] rotNe=new Vect[this.nElEdge];


		double Jr=1.0/(jac.determinant());
		Mat P=new Mat(3,2);
		P.el[0][0]=jac.el[2][1]-jac.el[2][2];
		P.el[0][1]=jac.el[1][2]-jac.el[1][1];
		P.el[1][0]=jac.el[2][2]-jac.el[2][0];
		P.el[1][1]=jac.el[1][0]-jac.el[1][2];
		P.el[2][0]=jac.el[2][0]-jac.el[2][1];
		P.el[2][1]=jac.el[1][1]-jac.el[1][0];

		P=P.times(Jr);

		Vect[] grad=new Vect[3];

		for(int i=0;i<3;i++){
			grad[i]=P.rowVect(i);

			rotNe[i]=new Vect(grad[i].el[1], -grad[i].el[0]);

		}

		return rotNe;


	}

	private Vect[] rotNe3ang2nd(Mat jac, Vect lc){

		double[] tu=lc.el;

		Vect[] rotNe=new Vect[6];


		double Jr=1.0/(jac.determinant());
		Mat P=new Mat(3,2);
		P.el[0][0]=jac.el[2][1]-jac.el[2][2];
		P.el[0][1]=jac.el[1][2]-jac.el[1][1];
		P.el[1][0]=jac.el[2][2]-jac.el[2][0];
		P.el[1][1]=jac.el[1][0]-jac.el[1][2];
		P.el[2][0]=jac.el[2][0]-jac.el[2][1];
		P.el[2][1]=jac.el[1][1]-jac.el[1][0];

		P=P.times(Jr);


		Vect[] grad=new Vect[3];

		for(int i=0;i<3;i++){
			grad[i]=P.rowVect(i).times(4*tu[i]-1);
		}

		grad[3]=P.rowVect(0).times(4*tu[1]).add(P.rowVect(1).times(4*tu[0]));
		grad[4]=P.rowVect(1).times(4*tu[2]).add(P.rowVect(2).times(4*tu[1]));
		grad[5]=P.rowVect(0).times(4*tu[2]).add(P.rowVect(2).times(4*tu[0]));



		for(int i=0;i<3;i++){
			rotNe[i]=new Vect(grad[i].el[1], -grad[i].el[0]);

		}


		return rotNe;


	}


	private Vect[] rotNe3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		return rotNe3ang(vertexNode);

	}

	private double[] Ne3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		return Ne3ang(vertexNode);

	}

	public double[] N3ang2nd(Vect localCo){

		double[] N=new double[6];
		double[] tu=localCo.el;

		for(int i=0;i<3;i++)
			N[i]=tu[i]*(2*tu[i]-1);

		N[3]=4*tu[0]*tu[1];
		N[4]=4*tu[1]*tu[2];
		N[5]=4*tu[0]*tu[2];

		return N;
	}

	public double[] Ne3ang2nd(Vect localCo){
		return N3ang2nd(localCo);
	}

	public Vect[] rotNe3ang2nd(Node[] vertexNode,Vect globalCo){

		Vect v1=vertexNode[0].getCoord().sub(globalCo);
		Vect v2=vertexNode[1].getCoord().sub(globalCo);
		Vect v3=vertexNode[2].getCoord().sub(globalCo);

		double[] N=new double[6];

		double[] su=new double [3];
		su[0]=v2.cross(v3).norm();
		su[1]=v3.cross(v1).norm();
		su[2]=v1.cross(v2).norm();

		double S=su[0]+su[1]+su[2];
		double rS=1.0/S;
		double[] tu=new double [3];
		for(int i=0;i<3;i++)
			tu[i]=su[i]*rS;

		Vect[] rotNe=new Vect[6];

		Vect[] rn=new Vect[3];
		rn[0]=new Vect(v3.el[0]-v2.el[0], v3.el[1]-v2.el[1]).times(rS);
		rn[1]=new Vect(v1.el[0]-v3.el[0], v1.el[1]-v3.el[1]).times(rS);
		rn[2]=new Vect(v2.el[0]-v1.el[0],v2.el[1]-v1.el[1]).times(rS);
		rotNe[0]=rn[0].times(4*tu[0]-1);
		rotNe[1]=rn[1].times(4*tu[1]-1);
		rotNe[2]=rn[2].times(4*tu[2]-1);
		rotNe[3]=rn[0].times(4*tu[1]).add(rn[1].times(4*tu[0]));
		rotNe[4]=rn[1].times(4*tu[2]).add(rn[2].times(4*tu[1]));
		rotNe[5]=rn[2].times(4*tu[0]).add(rn[0].times(4*tu[2]));


		return rotNe;
	}


	public Mat jacobian(Node[] vertexNode,Vect localCo){
		if(this.elCode==0) return jacobian3ang(vertexNode,localCo);
		if(this.elCode==2) return jacobianTetra(vertexNode,localCo);
		if(this.elCode==3) return jacobianPrism(vertexNode,localCo);



		Mat J=new Mat(this.dim,this.dim);
		Vect[] gN=new Vect[this.nElVert];

		gN=localGradN(localCo);

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}

		return J;
	}



	public Mat jacobian3ang(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);
		Vect[] gN=new Vect[3];

		gN=localGradN3ang();

		for(int i=0;i<3;i++) 
			J.el[0][i]=1;

		for(int i=0;i<3;i++) {
			for(int j=1;j<3;j++)
				for(int k=0;k<3;k++)
					J.el[j][k]+=gN[i].el[k]* vertexNode[i].getCoord(j-1);
		}


		return J;
	}

	public Mat jacobianTetra(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);

		/*		Vect[] gN=this.localGradNtet();

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}*/

		/*		for(int i=0;i<this.nElVert;i++) 
			for(int j=0;j<this.nElVert;j++)
				if(i==0)
					J.el[i][j]=1;
				else
				J.el[i][j]=vertexNode[j].getCoord(i-1);*/


		J.el[0][0]=vertexNode[0].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][0]=vertexNode[0].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][0]=vertexNode[0].getCoord(2)-vertexNode[3].getCoord(2);
		J.el[0][1]=vertexNode[1].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][1]=vertexNode[1].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][1]=vertexNode[1].getCoord(2)-vertexNode[3].getCoord(2);		
		J.el[0][2]=vertexNode[2].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][2]=vertexNode[2].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][2]=vertexNode[2].getCoord(2)-vertexNode[3].getCoord(2);



		return J;
	}

	public Mat jacobianPrism(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);

		Vect[] gN=new Vect[this.nElVert];

		gN=localGradNPrism(localCo);

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}

		return J;
	}


	public double[][] Qe(Model model,int ie){

		if(model.elCode==0) return Qe3ang(model,ie);

		Node[] vertexNode=model.elementNodes(ie);
		double[][] M=new double[this.nElVert][this.nElVert];
		Vect dtSigma=model.element[ie].getSigma();

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(3);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M[i][j]+=ws*gradN[i].dot(gradN[j].times(dtSigma));
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
			{
				M[j][i]=M[i][j];
			}

		return M;

	}

	public Mat elemPhiMat(Model model,int ie){

		if(model.elCode==1) return elemPhiMatQuad(model, ie);
			
		Node[] vertexNode=model.elementNodes(ie);
		Mat M=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(3);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M.el[i][j]+=ws*gradN[i].dot(gradN[j]);
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
			{
				M.el[j][i]=M.el[i][j];
			}

		return M;

	}
	
	public Mat elemPhiMatQuad(Model model,int ie){

			
		Node[] vertexNode=model.elementNodes(ie);
		Mat M=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(model.dim);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M.el[i][j]+=ws*gradN[i].dot(gradN[j]);
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
			{
				M.el[j][i]=M.el[i][j];
			}

		return M;

	}
	

	public Vect elemPhiVect(Model model,int ie){

		if(model.elCode==1) return elemPhiVectQuad(model, ie);
		
		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect elemVec=new Vect(this.nElVert);
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		
		Vect localCo=new Vect(3);
		Vect[] Ne=new Vect[this.nElEdge];
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);

					Ne=Ne(jac,localCo,edgeDir);
					
					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
							elemVec.el[i]+=ws*gradN[i].dot(Ne[j].times(A[j]));
				}	


		return elemVec;

	}
	
	public Vect elemPhiVectQuad(Model model,int ie){

		
		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect elemVec=new Vect(this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect localCo=new Vect(model.dim);
		double[] N=new double[this.nElEdge];
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();
		
		double gradN=1./model.height;
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
					else
						ws=detJac;


					N=N(localCo);
					

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
									elemVec.el[i]+=ws*gradN*N[j]*A[j];
				}	
		

		return elemVec;

	}

	public double[][] Qe3ang(Model model, int ie){
		double rdt=1.0/model.dt;

		double[][] Qe=new double[this.nElVert][this.nElEdge];

		double S=el3angArea(model,ie);
		double sigmaZ=model.element[ie].getSigma().el[2]*rdt;


		Vect[] gradN=gradN3ang(model,ie);


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElEdge;j++)	{
				Qe[i][j]=S*gradN[i].dot(gradN[j].times(sigmaZ));

			}



		return Qe;


	}


	public double[][] Pe(Model model, int ie){

		Vect sigma=model.element[ie].getSigma();

		Node[] vertexNode=model.elementNodes(ie);
		double[][] M=new double[this.nElVert][this.nElEdge];

		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		double detJac,ws=1;

		Mat jac;

		Vect[] Ne=new Vect[this.nElEdge];
		Vect[] gradN=new Vect[this.nElVert];

		Vect localCo=new Vect(3);
		int n=this.PW[0].length; 
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;

					Ne=Ne(jac,localCo,edgeDir);


					gradN=gradN(jac,localCo);

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
							M[i][j]+=ws*gradN[i].dot(Ne[j].times(sigma));

				}	

		return M;
	}

	public Vect gradPhi(Node[] vertexNode,double[] nodePhi){

		Vect gradPhi=new Vect(3);
		Vect zero=new Vect(3);
		Vect[] gradN=new Vect[this.nElVert];
		Mat G=jacobian(vertexNode,zero).inv3();
		Vect[] localGradN=localGradN(zero);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=G.mul(localGradN[i]);
			gradPhi= gradPhi.add(gradN[i].times(nodePhi[i]));
		}

		return  gradPhi;	

	}

	public Mat[][] Ke(Model model,int ie){

		if(this.elCode==0) return Ke3ang(model,ie);
		else if(this.elCode==1) return KeQuad(model,ie);
		else if(elCode==3) {return KePrism(model,ie);}
		Node[] vertexNode=model.elementNodes(ie);
		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(3,3);

		boolean isotElast=model.region[model.element[ie].getRegion()].isotElast;
		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		Vect v=model.element[ie].getPois();
		Vect E=model.element[ie].getYng();
		Vect G=model.element[ie].getShear();




		Mat BtDB;
		Vect[] gradN;
		Vect localCo=new Vect(3);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());


					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;

					gradN=gradN(jac,localCo);

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElVert;j++){

							if(isotElast)
								BtDB=BtDB(v.el[0],G.el[0],gradN[i],gradN[j]).times(ws);
							else
								BtDB=BtDB(E,v,G,gradN[i],gradN[j]).times(ws);
							Ke[i][j]=Ke[i][j].add(BtDB);

						}


				}	



		return Ke;
	}
	
	
	
	public Mat[][] Ke_tang(Model model,int ie,boolean[] yield_state){

		if(this.elCode==0) return null;
		else if(this.elCode==1) return KeQuad_tang(model,ie, yield_state);
		else return null;
	}
	

	public Mat[][] KeTet(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);


		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(3,3);

		boolean isotElast=model.region[model.element[ie].getRegion()].isotElast;

		Vect v=model.element[ie].getPois();

		Vect E=model.element[ie].getYng();
		Vect G=new Vect(E.length);


		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));
		Vect[] gradN=this.gradNTet(jac);

		double ws=abs(jac.determinant())/6;

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				if(isotElast)
					Ke[i][j]=BtDB(v.el[0],G.el[0],gradN[i],gradN[j]).times(ws);

				else
					Ke[i][j]=BtDB(E,v,G,gradN[i],gradN[j]).times(ws);

				//Ke[i][j].show();
			}




		return Ke;
	}

	public Mat[][] Ke(Model model,int ie, Mat[][] Me){

		if(elCode==0) { fillMe3ang(model,ie,Me); return Ke3ang(model,ie);}
		else if(elCode==1)   { fillMeQuad(model,ie,Me); return KeQuad(model,ie);}
		else if(elCode==2) { fillMeTet(model,ie,Me); return KeTet(model,ie);}
		else if(elCode==3) { fillMePrism(model,ie,Me); return KePrism(model,ie);}


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me[i][j]=new Mat(3,3);


		Node[] vertexNode=model.elementNodes(ie);
		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		Mat M1=new Mat(this.nElVert,this.nElVert);
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(3,3);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		boolean isotElast=model.region[model.element[ie].getRegion()].isotElast;

		Vect v=model.element[ie].getPois();
		Vect E=model.element[ie].getYng();
		Vect G=model.element[ie].getShear();

		double ro=model.element[ie].getRo();
		Mat BtDB;
		Vect[] gradN;
		double[] N;
		Mat Nv=new Mat(1,this.nElVert);
		Vect localCo=new Vect(3);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());


					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;

					gradN=gradN(jac,localCo);

					N=N(localCo);

					for(int j=0;j<this.nElVert;j++)
						Nv.el[0][j]=N[j];

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElVert;j++){

							if(isotElast)
								BtDB=BtDB(v.el[0],G.el[0],gradN[i],gradN[j]).times(ws);
							else
								BtDB=BtDB(E,v,G,gradN[i],gradN[j]).times(ws);
							Ke[i][j]=Ke[i][j].add(BtDB);
						}



					M1=M1.add(Nv.transp().mul(Nv).times(ws));

				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				for(int ir=0;ir<3;ir++)
				{
					Me[i][j].el[ir][ir]=ro*M1.el[i][j];

				}

			}


		return Ke;
	}	


	public Mat[][] KeQuad(Model model,int ie){
		
		if(struc2D==2)  return KeQuadAxi(model,ie);


		Node[] vertexNode=model.elementNodes(ie);

		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(2,2);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];
		double G=0;
		if(struc2D==0)
			G=E/(1-v*v);
		else
			G=E/((1+v)*(1-2*v));

		Mat BtDB;
		Vect[] gradN;
		Vect localCo=new Vect(3);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				
				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());


				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;

				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++)
					for(int j=0;j<this.nElVert;j++){

						BtDB=BtDB2D(v,gradN[i],gradN[j]).times(ws);
						Ke[i][j]=Ke[i][j].add(BtDB);

					}
			}	

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				Ke[i][j]=Ke[i][j].times(G);

			}
		return Ke;
	}
	
	
public Mat[][] KeQuad_tang(Model model,int ie, boolean[] yield_state){
		

		Node[] vertexNode=model.elementNodes(ie);

		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(2,2);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		int nr=model.element[ie].getRegion();

		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];
		
		double Et=model.region[nr].getTangYoung();
		
		double Gi=0;
		if(struc2D==0)
			Gi=E/(1-v*v);
		else
			Gi=E/((1+v)*(1-2*v));
		
		double Gt=Gi*Et/E;

		Mat BtDB;
		Vect[] gradN;
		Vect localCo=new Vect(3);
		
		int count=0;
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				
				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				
				boolean y_state=yield_state[count++];
			
				
				double G=Gi;
				if(y_state) G=Gt;
				

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;

				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++)
					for(int j=0;j<this.nElVert;j++){

						BtDB=BtDB2D(v,gradN[i],gradN[j]).times(ws*G);
						Ke[i][j]=Ke[i][j].add(BtDB);

					}
			}	


		return Ke;
	}
	
	public Mat[][] KeQuadAxi(Model model,int ie){


		Node[] vertexNode=model.elementNodes(ie);
		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(2,2);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		boolean rcent=false;
		
		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];
		double G=E/((1+v)*(1-2*v));

		double b=v;
		double a=(1-v);
		double c=(0.5-v);
		Mat D=new Mat(4,4);
		D.el[0][0]=a;
		D.el[0][1]=b;
		D.el[1][0]=b;
		D.el[1][1]=a;
		D.el[2][2]=a;
		D.el[2][0]=b;
		D.el[2][1]=b;
		D.el[0][2]=b;
		D.el[1][2]=b;
		D.el[3][3]=c;
		
		
		Mat BtDB;
		Vect[] gradN;
		double [] Ne=new double[model.nElVert];
		
		double[] ri=new double[model.nElVert];
		for(int i=0;i<model.nElVert;i++){
			double r1=vertexNode[i].getCoord(0);
			ri[i]=r1;
		}

		
		Vect localCo=new Vect(3);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
	
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());




				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;

				gradN=gradN(jac,localCo);
				
				Ne=NQuad(localCo);


				double r=0;
				
			if(!rcent){

					for(int i=0;i<model.nElEdge;i++){
						r+=Ne[i]*ri[i];
					}
				}
				else{

					for(int i=0;i<model.nElEdge;i++){
						r+=.25*ri[i];
					}
				}
	
				
				ws*=r;
				
				for(int i=0;i<this.nElVert;i++)
					for(int j=0;j<this.nElVert;j++){

							
							Mat Bj=new Mat(4,2);
							
							Bj.el[0][0]=gradN[j].el[0];
							Bj.el[1][1]=gradN[j].el[1];
							
							Bj.el[2][0]=Ne[j]/r;
							Bj.el[3][0]=gradN[j].el[1];
							Bj.el[3][1]=gradN[j].el[0];
	
								Mat Bi=new Mat(4,2);
								
								Bi.el[0][0]=gradN[i].el[0];
								Bi.el[1][1]=gradN[i].el[1];
								
								Bi.el[2][0]=Ne[i]/r;
								Bi.el[3][0]=gradN[i].el[1];
								Bi.el[3][1]=gradN[i].el[0];
								
						BtDB=Bi.transp().mul(D.mul(Bj)).times(ws);

						Ke[i][j]=Ke[i][j].add(BtDB);

					}
			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				Ke[i][j]=Ke[i][j].times(2*Math.PI*G);
			}
		

		return Ke;
	}





	public Mat[][] Ke3ang(Model model,int ie){
		
	if(struc2D==2) return Ke3angAxi(model,ie);

		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(2,2);

		double S=el3angArea(model,ie);

		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];

		double G=0;
		if(struc2D==0)
			G=E/(1-v*v);
		else
			G=E/((1+v)*(1-2*v));

		double GS=G*S;

		Vect[] gradN=gradN3ang(model,ie);

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){

				Ke[i][j]=BtDB2D(v,gradN[i],gradN[j]);

			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				Ke[i][j]=Ke[i][j].times(GS);

			}

		return Ke;
	}
	
	public Mat[][] Ke3angAxi(Model model,int ie){

		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(2,2);

		double S=el3angArea(model,ie);

		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];

		double G=E/((1+v)*(1-2*v));

		double b=v;
		double a=(1-v);
		double c=(0.5-v);
		Mat D=new Mat(4,4);
		D.el[0][0]=a;
		D.el[0][1]=b;
		D.el[1][0]=b;
		D.el[1][1]=a;
		D.el[2][2]=a;
		D.el[2][0]=b;
		D.el[2][1]=b;
		D.el[0][2]=b;
		D.el[1][2]=b;
		D.el[3][3]=c;

		double GS=G*S;
		
		Mat BtDB;

		Vect[] gradN=gradN3ang(model,ie);
		
		double[] N=N3ang(model,ie);

		Node[] vertexNode=model.elementNodes(ie);

		double r=0;
		

		for(int i=0;i<model.nElVert;i++){
			double r1=vertexNode[i].getCoord(0);
			r+=0.33333333*r1;
		}


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				
				Mat Bj=new Mat(4,2);
				
				Bj.el[0][0]=gradN[j].el[0];
				Bj.el[1][1]=gradN[j].el[1];
				
				Bj.el[2][0]=N[j]/r;
				Bj.el[3][0]=gradN[j].el[1];
				Bj.el[3][1]=gradN[j].el[0];

					Mat Bi=new Mat(4,2);
					
					Bi.el[0][0]=gradN[i].el[0];
					Bi.el[1][1]=gradN[i].el[1];
					
					Bi.el[2][0]=N[i]/r;
					Bi.el[3][0]=gradN[i].el[1];
					Bi.el[3][1]=gradN[i].el[0];
					
					BtDB=Bi.transp().mul(D.mul(Bj)).times(r);
					
					Ke[i][j]=BtDB.times(GS*2*Math.PI);

			}


		return Ke;
	}


	private void fillMeQuad(Model model,int ie,Mat[][] Me){

		Node[] vertexNode=model.elementNodes(ie);


		Mat M1=new Mat(this.nElVert,this.nElVert);
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me[i][j]=new Mat(2,2);

		double detJac,ws=1;
		Mat jac;
		double ro=model.element[ie].getRo();
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;
				M1=M1.add(Nv.transp().mul(Nv).times(ws));


			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				for(int ir=0;ir<dim;ir++)
				{
					Me[i][j].el[ir][ir]=ro*M1.el[i][j];

				}

			}

	}

	private void fillMe3ang(Model model,int ie,Mat[][] Me){

		Node[] vertexNode=model.elementNodes(ie);
		Mat M1=new Mat(this.nElVert,this.nElVert);
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me[i][j]=new Mat(2,2);

		double detJac,ws=1;

		Mat jac;
		double ro=model.element[ie].getRo();
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];
			localNe[0]=lc.el[0];
			localNe[1]=lc.el[1];
			localNe[2]=1-lc.el[0]-lc.el[1];


			for(int j=0;j<this.nElVert;j++)
				Nv.el[0][j]=localNe[j];

			jac=jacobian3ang(vertexNode,lc);

			detJac=abs(jac.determinant());
			ws=this.PW3ang[p][2]*detJac;
			M1=M1.add(Nv.transp().mul(Nv).times(ws));
		}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				for(int ir=0;ir<2;ir++)
				{
					Me[i][j].el[ir][ir]=ro*M1.el[i][j];

				}

			}

	}

	private void fillMeTet(Model model,int ie,Mat[][] Me){

		Node[] vertexNode=model.elementNodes(ie);
		Mat M1=new Mat(this.nElVert,this.nElVert);
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me[i][j]=new Mat(3,3);

		double detJac,ws=1;

		Mat jac;
		double ro=model.element[ie].getRo();
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PWtetra.length;	
		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];


			for(int j=0;j<this.nElVert-1;j++)
				Nv.el[0][j]=lc.el[j];

			Nv.el[0][3]=1-lc.el[0]-lc.el[1]-lc.el[2];

			jac=jacobianTetra(vertexNode,lc);

			detJac=abs(jac.determinant());
			ws=this.PWtetra[p][3]*detJac;
			M1=M1.add(Nv.transp().mul(Nv).times(ws));
		}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				for(int ir=0;ir<3;ir++)
				{
					Me[i][j].el[ir][ir]=ro*M1.el[i][j];

				}

			}

	}


	private void fillMePrism(Model model,int ie,Mat[][] Me){

		Node[] vertexNode=model.elementNodes(ie);

		Mat M1=new Mat(this.nElVert,this.nElVert);
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me[i][j]=new Mat(3,3);

		double detJac,ws=1;
		Mat jac;
		double ro=model.element[ie].getRo();

		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		int nr=this.PW[0].length;	
		for(int r=0;r<nr;r++)
			for(int p=0;p<n;p++)
			{

				lc.el[0]=this.PW3ang[p][0];
				lc.el[1]=this.PW3ang[p][1];
				lc.el[2]=this.PW[0][r];

				localNe=localNPrims(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobianPrism(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW3ang[p][2]*this.PW[1][r]*detJac;
				M1=M1.add(Nv.transp().mul(Nv).times(ws));


			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
				for(int ir=0;ir<dim;ir++)
				{
					Me[i][j].el[ir][ir]=ro*M1.el[i][j];

				}

			}

	}


	public Mat[][] KePrism(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);

		Mat[][] Ke=new Mat[this.nElVert][this.nElVert];
		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Ke[i][j]=new Mat(3,3);

		Vect v=model.element[ie].getPois();
		Vect E=model.element[ie].getYng();
		Vect G=model.element[ie].getShear();

		boolean isotElast=model.region[model.element[ie].getRegion()].isotElast;




		Vect[] gradN;

		int nr=this.PW[0].length;
		int n=this.PW3ang.length;
		Vect lc=new Vect(3);
		Mat jac;
		double detJac,ws;
		for(int r=0;r<nr;r++)
			for(int p=0;p<n;p++)
			{

				lc.el[0]=this.PW3ang[p][0];
				lc.el[1]=this.PW3ang[p][1];				
				lc.el[2]=this.PW[0][r];

				jac=jacobianPrism(vertexNode,lc);

				detJac=abs(jac.determinant());
				ws=this.PW3ang[p][2]*this.PW[1][r]*detJac;


				gradN=gradN(jac,lc);

				for(int i=0;i<6;i++)
					for(int j=0;j<6;j++){


						if(isotElast)
							Ke[i][j]=Ke[i][j].add(BtDB(v.el[0],G.el[0],gradN[i],gradN[j]).times(ws));
						else
							Ke[i][j]=Ke[i][j].add(BtDB(E,v,G,gradN[i],gradN[j]).times(ws));
					}
			}


		return Ke;
	}




	public Mat Se(Model model,int ie){

		if(this.elCode==0) return Se3ang(model,ie);
		else if(this.elCode==1) return SeQuad(model,ie);
		else return null;
	}

	public Mat Se(Model model,int ie,Mat Me){

		if(this.elCode==0) { fillMe3angScalar(model,ie,Me); return Se3ang(model,ie);}

		if(this.elCode==1) { fillMeQuadScalar(model,ie,Me); return SeQuad(model,ie);}
		else return null;
	}



	public Mat SeQuad(Model model,int ie){


		Node[] vertexNode=model.elementNodes(ie);
		Mat Ke=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		Vect kt=model.element[ie].getSigma();

		Vect[] gradN;
		Vect localCo=new Vect(dim);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(kt.times(gradN[j]))*ws;
						Ke.el[i][j]+=BtkB;

					}
				}

			}	


		return Ke;
	}	



	public Mat gradNi_gradNjQ(Model model,int ie){


		Node[] vertexNode=model.elementNodes(ie);
		Mat Ke=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		Vect[] gradN;
		Vect localCo=new Vect(dim);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(gradN[j])*ws;
						Ke.el[i][j]+=BtkB;

					}
				}

			}	


		return Ke;
	}

	public Mat NiNjQ(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		Mat Me=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		Mat jac;
		Mat M1=new Mat(this.nElVert,this.nElVert);
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;

				M1=M1.add(Nv.transp().mul(Nv).times(ws));



			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me.el[i][j]=M1.el[i][j];

		return Me;
	}



	public Mat Se3ang(Model model,int ie){

		Mat Ke=new Mat(this.nElVert,this.nElVert);	

		Vect kt=model.element[ie].getSigma();
		double S=el3angArea(model,ie);


		Vect[] gradN=gradN3ang(model,ie);

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){

				Ke.el[i][j]=gradN[i].dot(kt.times(gradN[j]))*S;

			}


		return Ke;
	}


	private void fillMe3angScalar(Model model,int ie,Mat Me){

		Node[] vertexNode=model.elementNodes(ie);

		double detJac,ws=1;
		Mat jac;
		double gw=model.gw;
		double he=model.getElementT(ie);
		double lam=gw*Math.exp(-.05*he);
		Mat M1=new Mat(this.nElVert,this.nElVert);

		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];
			localNe[0]=lc.el[0];
			localNe[1]=lc.el[1];
			localNe[2]=1-lc.el[0]-lc.el[1];


			for(int j=0;j<this.nElVert;j++)
				Nv.el[0][j]=localNe[j];

			jac=jacobian3ang(vertexNode,lc);

			detJac=abs(jac.determinant());
			ws=this.PW3ang[p][2]*detJac;
			M1=M1.add(Nv.transp().mul(Nv).times(ws));
		}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){

				Me.el[i][j]=lam*M1.el[i][j];

			}

	}

	private void fillMeQuadScalar(Model model,int ie,Mat Me){

		Node[] vertexNode=model.elementNodes(ie);


		double detJac,ws=1;
		Mat jac;
		double gw=model.gw;
		double he=model.getElementT(ie);
		double lam=gw*Math.exp(-.05*he);
		Mat M1=new Mat(this.nElVert,this.nElVert);
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;

				M1=M1.add(Nv.transp().mul(Nv).times(ws));



			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me.el[i][j]=lam*M1.el[i][j];



	}

	Mat BtDB(Vect yng,Vect pois, Vect shear,Vect gNi,Vect gNj){


		Mat D=hook( yng, pois,  shear);



		Mat Kn=new Mat(3,3);	

		double pipj=gNi.el[0]*gNj.el[0];
		double qiqj=gNi.el[1]*gNj.el[1];
		double rirj=gNi.el[2]*gNj.el[2];

		double piqj=gNi.el[0]*gNj.el[1];
		double qipj=gNi.el[1]*gNj.el[0];

		double pirj=gNi.el[0]*gNj.el[2];
		double ripj=gNi.el[2]*gNj.el[0];

		double qirj=gNi.el[1]*gNj.el[2];
		double riqj=gNi.el[2]*gNj.el[1];

		Kn.el[0][0]=D.el[0][0]*pipj+D.el[3][3]*qiqj+D.el[5][5]*rirj;
		Kn.el[0][1]=D.el[0][1]*piqj+D.el[3][3]*qipj;
		Kn.el[0][2]=D.el[0][2]*pirj+D.el[5][5]*ripj;
		Kn.el[1][0]=D.el[0][1]*qipj+D.el[3][3]*piqj;
		Kn.el[1][1]=D.el[1][1]*qiqj+D.el[3][3]*pipj+D.el[4][4]*rirj;
		Kn.el[1][2]=D.el[1][2]*qirj+D.el[4][4]*riqj;
		Kn.el[2][0]=D.el[0][2]*ripj+D.el[5][5]*pirj;
		Kn.el[2][1]=D.el[1][2]*riqj+D.el[4][4]*qirj;
		Kn.el[2][2]=D.el[2][2]*rirj+D.el[4][4]*pipj+D.el[5][5]*qiqj;

		//util.pr(D.el[5][5]);

		return Kn;
	}



	Mat BtDB(double v,double G,Vect gNi,Vect gNj){

		Mat Kn=new Mat(3,3);	
		double a=2*(1-v)/(1-2*v)*G;
		double b=2*v/(1-2*v)*G;;
		double c=G;

		double pipj=gNi.el[0]*gNj.el[0];
		double qiqj=gNi.el[1]*gNj.el[1];
		double rirj=gNi.el[2]*gNj.el[2];

		double piqj=gNi.el[0]*gNj.el[1];
		double qipj=gNi.el[1]*gNj.el[0];

		double pirj=gNi.el[0]*gNj.el[2];
		double ripj=gNi.el[2]*gNj.el[0];

		double qirj=gNi.el[1]*gNj.el[2];
		double riqj=gNi.el[2]*gNj.el[1];

		Kn.el[0][0]=a*pipj+c*(qiqj+rirj);
		Kn.el[0][1]=b*piqj+c*qipj;
		Kn.el[0][2]=b*pirj+c*ripj;
		Kn.el[1][0]=b*qipj+c*piqj;
		Kn.el[1][1]=a*qiqj+c*(pipj+rirj);
		Kn.el[1][2]=b*qirj+c*riqj;
		Kn.el[2][0]=b*ripj+c*pirj;
		Kn.el[2][1]=b*riqj+c*qirj;
		Kn.el[2][2]=a*rirj+c*(pipj+qiqj);

		//util.pr(c);

		return Kn;
	}



	public Mat BtDB2D(double v,Vect gNi,Vect gNj){

		Mat Kn=new Mat(2,2);

		double a,b,c;
		if(struc2D==0){
			b=v;
			a=1;
			c=(1.0-v)/2;
		}
		else {
			b=v;
			a=(1-v);
			c=(0.5-v);
		}


		double pipj=gNi.el[0]*gNj.el[0];
		double qiqj=gNi.el[1]*gNj.el[1];

		double piqj=gNi.el[0]*gNj.el[1];
		double qipj=gNi.el[1]*gNj.el[0];


		Kn.el[0][0]=a*pipj+c*qiqj;
		Kn.el[0][1]=b*piqj+c*qipj;
		Kn.el[1][0]=b*qipj+c*piqj;
		Kn.el[1][1]=a*qiqj+c*pipj;

		return Kn;
	}



	public Mat hookNormal(double E,double v){

		Mat D=new Mat(this.dim,this.dim);
		for(int i=0;i<this.dim;i++)
			for(int j=0;j<this.dim;j++)
				if(i==j) D.el[i][j]=1-v;
				else
					D.el[i][j]=v;

		return D.times(E/((1-2*v)*(1+v)));


	}



	public Mat hook3D(Model model, int ie){
		Vect yng=model.element[ie].getYng();
		Vect pois=model.element[ie].getPois();
		Vect shear=model.element[ie].getShear();
		return  hook( yng, pois,  shear);
	}

	public Mat hook(Vect yng,Vect pois, Vect shear){
		if(yng.length==1 || yng.length==2)
			return hookIsot(yng.el[0],pois.el[0]);

		double vxy=pois.el[0];
		double vyz=pois.el[1];
		double vxz=pois.el[2];





		double vyx=yng.el[1]/yng.el[0]*vxy;
		double vzy=yng.el[2]/yng.el[1]*vyz;
		double vzx=yng.el[2]/yng.el[0]*vxz;




		double delta=(1-vxy*vyx-vyz*vzy-vxz*vzx-2*vxy*vyz*vzx)/(yng.el[0]*yng.el[1]*yng.el[2]);

		Mat D=new Mat(6,6);

		double denum=(delta*yng.el[1]*yng.el[2]);
		D.el[0][0]=(1-vyz*vzy)/denum;
		D.el[0][1]=(vyx+vzx*vyz)/denum;
		D.el[0][2]=(vzx+vyx*vzy)/denum;

		denum=(delta*yng.el[0]*yng.el[2]);
		D.el[1][1]=(1-vzx*vxz)/denum;
		D.el[1][2]=(vzy+vzx*vxy)/denum;

		denum=(delta*yng.el[0]*yng.el[1]);
		D.el[2][2]=(1-vxy*vyx)/denum;

		D.el[3][3]=shear.el[0];
		D.el[4][4]=shear.el[1];
		D.el[5][5]=shear.el[2];

		return D;

	}
	public Mat hookIsot(double E,double v){



		if(this.dim==3 || struc2D==1){
		int m=3*(this.dim-1);
		Mat D=new Mat(m,m);
		for(int i=0;i<this.dim;i++)
			for(int j=0;j<this.dim;j++)
				if(i==j) D.el[i][j]=1-v;
				else
					D.el[i][j]=v;

		for(int i=this.dim;i<m;i++)
			D.el[i][i]=.5-v;
		
		return D.times(E/((1-2*v)*(1+v)));
		
		}else{
			double	G=E/(1-v*v);
		
			Mat D=new Mat(3,3);;
			
			D.el[0][0]=1;
			D.el[1][1]=1;
			D.el[2][2]=(1.-v)/2;
			D.el[0][1]=v;
			D.el[1][0]=v;
			
			D=D.times(G);

			return D;
		}

		


	}

	public Mat hook(Model model,int ie){

		return hook(model.element[ie].getYng(),model.element[ie].getPois(),model.element[ie].getShear() );

		//	return hook(model.element[ie].getYng().el[0],model.element[ie].getPois().el[0] );

	}

	public Mat BtDdBdA(int k,double v,Vect gNi,Vect gNj){

		Mat Kn=new Mat(3,3);
		double b=v;
		double a=(1-v);
		double c=(0.5-v);

		double pipj=gNi.el[0]*gNj.el[0];
		double qiqj=gNi.el[1]*gNj.el[1];
		double rirj=gNi.el[2]*gNj.el[2];

		double piqj=gNi.el[0]*gNj.el[1];
		double qipj=gNi.el[1]*gNj.el[0];

		double pirj=gNi.el[0]*gNj.el[2];
		double ripj=gNi.el[2]*gNj.el[0];

		double qirj=gNi.el[1]*gNj.el[2];
		double riqj=gNi.el[2]*gNj.el[1];

		Kn.el[0][0]=a*pipj+c*(qiqj+rirj);
		Kn.el[0][1]=b*piqj+c*qipj;
		Kn.el[0][2]=b*pirj+c*ripj;
		Kn.el[1][0]=b*qipj+c*piqj;
		Kn.el[1][1]=a*qiqj+c*(pipj+rirj);
		Kn.el[1][2]=b*qirj+c*riqj;
		Kn.el[2][0]=b*ripj+c*pirj;
		Kn.el[2][1]=b*riqj+c*qirj;
		Kn.el[2][2]=a*rirj+c*(pipj+qiqj);

		return Kn;
	}


	public Vect getStrain(Model model,int ie){
		Vect lc=new Vect(this.dim);

		return getStrain(model,ie,lc);
	}

	public Vect getStrain(Model model,int ie,Vect lc){
		if(this.elCode==0) return getEl3angStrain(model,ie);
		else if(this.elCode==1) return getElQuadStrain(model,ie,lc);

		Node[] vertexNode=model.elementNodes(ie);
		Mat jac=jacobian(vertexNode,lc);

		Vect[] gradN=gradN(jac,lc);
		Vect strain=new Vect(6);
		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<this.nElVert;j++){
			Vect u=model.node[vertNumb[j]].getU();


			if(model.coordCode==1){
				Mat R2D=util.rotMat2D(util.getAng(model.node[vertNumb[j]].getCoord().v2()));
				Mat R=new Mat(dim,dim);
				for(int m=0;m<2;m++)
					for(int n=0;n<2;n++)
						R.el[m][n]=R2D.el[m][n];

				R.el[2][2]=1;
				u=R.mul(u);
			}

			strain.el[0]+=u.el[0]*gradN[j].el[0];
			strain.el[1]+=u.el[1]*gradN[j].el[1];
			strain.el[2]+=u.el[2]*gradN[j].el[2];
			strain.el[3]+=u.el[0]*gradN[j].el[1]+u.el[1]*gradN[j].el[0];
			strain.el[4]+=u.el[1]*gradN[j].el[2]+u.el[2]*gradN[j].el[1];
			strain.el[5]+=u.el[0]*gradN[j].el[2]+u.el[2]*gradN[j].el[0];


		}


		return strain;
	}

	
	
	
	public Vect[] getGpStrainQuad(Model model,int ie){
		
		
		int n=this.PW[0].length;	
		
		 Vect[] strain=new Vect[n*n];
		Vect localCo=new Vect(2);
		int count=0;
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];
				
				strain[count++]=getElQuadStrain(model,ie,localCo);

			}
		
		return strain;

	}


	public Vect getElQuadStrain(Model model,int ie){

		Vect lc=new Vect(2);

		return getElQuadStrain(model,ie,lc);
	}


	public Vect getElQuadStrain(Model model,int ie,Vect lc){

		Node[] vertexNode=model.elementNodes(ie);

		Mat jac=jacobian(vertexNode,lc);

		Vect[] gradN=gradN(jac,lc);
		Vect strain=new Vect(3);
		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<this.nElVert;j++){
			Vect u=model.node[vertNumb[j]].getU();

			if(model.coordCode==1)
				u=util.rotMat2D(util.getAng(model.node[vertNumb[j]].getCoord())).mul(u);

			strain.el[0]+=u.el[0]*gradN[j].el[0];
			strain.el[1]+=u.el[1]*gradN[j].el[1];
			strain.el[2]+=u.el[0]*gradN[j].el[1]+u.el[1]*gradN[j].el[0];
		}


		return strain;
	}


	public Vect getEl3angStrain(Model model,int ie){
		/*Vect zero=new Vect(2);
		Mat jac=jacobian3ang(zero,vertexNode);*/
		Vect[] gradN=gradN3ang(model,ie);
		Vect strain=new Vect(3);
		int[] vertNumb=model.element[ie].getVertNumb();


		for(int j=0;j<this.nElVert;j++){
			Vect u=model.node[vertNumb[j]].getU();

			if(model.coordCode==1)
				u=util.rotMat2D(util.getAng(model.node[vertNumb[j]].getCoord())).mul(u);

			strain.el[0]+=u.el[0]*gradN[j].el[0];
			strain.el[1]+=u.el[1]*gradN[j].el[1];
			strain.el[2]+=u.el[0]*gradN[j].el[1]+u.el[1]*gradN[j].el[0];

		}

		return strain;
	}


	public Vect getElementCenter(Model model,int ie){


		Vect center=new Vect(this.dim);
		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<this.nElVert;j++)
			center=center.add(model.node[vertNumb[j]].getCoord());


		return center.times(1.0/this.nElVert);
	}


	public Vect getElementMS(Model model,int ie){

		int nLam=model.region[model.element[ie].getRegion()].lamBNumber;
		double lamb1=0;

		Vect B=model.element[ie].getB();
		/*		
		if(model.coupled){
		double sn=model.getStressB(ie, B);
			lamb1=model.lamBS[nLam].getLam(B.norm(),sn);
		}
		else*/
		lamb1=model.lamB[nLam].getLam(B.norm());
		Vect lamb=new Vect(this.dim);

		double vms=.5;

		lamb.el[0]=lamb1;
		lamb.el[1]=-vms*lamb.el[0];
		if(this.dim==3)
			lamb.el[2]=-vms*lamb.el[0];
		else{
			double v=model.element[ie].getPois().el[0];
			double lambx=v*lamb.el[1];
			lamb=lamb.add(lambx);
		}



		return lamb;
	}


	public Vect globalCo(Node[] vertex,double[]  localN){
		Vect v=new Vect(this.dim);
		for(int i=0;i<this.nElVert;i++) 
			v=v.add( vertex[i].getCoord().times(localN[i]));
		return v;
	}

	public  Vect localCo(Model model,int[] m,Vect globalCo){

		if(model.elCode==1) return localCoQuad(model,m[0],globalCo);
		Vect lc=new Vect(3);
		Vect dlc=new Vect(3);
		Vect gc=new Vect(3);
		double[] localN;

		Node[] vertex=model.elementNodes(m[0]);
		double resNorm=1;
		double error=1e-6;

		if(m[1]==0) lc.el[0]=-1;
		if(m[2]==0) lc.el[0]=1;
		if(m[3]==0) lc.el[1]=-1;
		if(m[4]==0) lc.el[1]=1;
		if(m[5]==0) lc.el[2]=-1;
		if(m[6]==0) lc.el[2]=1;

		boolean[] known=new boolean[3];
		for(int j=0;j<3;j++)
			if(lc.el[j]!=0) known[j]=true;

		for(int i=0;(i<1 && resNorm>error);i++){
			localN=N(lc);
			gc=globalCo(vertex,localN);

			Vect res=globalCo.sub(gc);

			resNorm=res.norm();
			for(int j=0;j<3;j++)
				if(!known[j])
					lc.el[j]+=dlc.el[j];
		}

		return lc;
	}

	private Vect localCoQuad(Model model,int ie,Vect globalCo){
		Vect lc=new Vect(2);
		Vect dlc=new Vect(2);
		Vect gc=new Vect(2);
		double[] localN;

		Node[] vertex=model.elementNodes(ie);
		double resNorm=1;
		double error=1e-12*model.scaleFactor;

		Vect[] v=new Vect[4];


		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<4;j++)
			v[j]=model.node[vertNumb[j]].getCoord();


		if(v[0].sub(globalCo).cross(v[1].sub(globalCo)).norm()<error) lc.el[1]=1;
		if(v[1].sub(globalCo).cross(v[2].sub(globalCo)).norm()<error) lc.el[0]=-1;
		if(v[2].sub(globalCo).cross(v[3].sub(globalCo)).norm()<error) lc.el[1]=-1;
		if(v[3].sub(globalCo).cross(v[0].sub(globalCo)).norm()<error) lc.el[0]=1;

		boolean[] known=new boolean[2];
		for(int j=0;j<2;j++)
			if(lc.el[j]!=0) known[j]=true;

		for(int i=0;(i<1 && resNorm>error);i++){
			localN=NQuad(lc);
			gc=globalCo(vertex,localN);

			Vect res=globalCo.sub(gc);

			resNorm=res.norm();
			for(int j=0;j<2;j++)
				if(!known[j])
					lc.el[j]+=dlc.el[j];
		}

		return lc;
	}

	public double[] getApAnAt(Model model,Vect globalCo){
		if(dim==3){
			throw new IllegalArgumentException(" for 3D not set up yet.");
		}

		if(model.elCode==0){
			throw new IllegalArgumentException(" for triangular element not set up yet.");
		}

		double na[]=new double[2];
		na[0]=-1e10;
		na[1]=-1e10;

		int[] m=model.getContainingElement(globalCo);
		if(m[0]<=0) return na;


		Vect	lc=localCoQuad(model,m[0],globalCo);


		return getApAnAt(model, m[0], lc);

	}

	public Vect getBAt(Model model,Vect globalCo){
		if(this.elCode==3) return new Vect(dim);

		Vect na=new Vect(this.dim);
		na.el[0]=1e10;
		na.el[1]=-1e10;
		int[] m=model.getContainingElement(globalCo);
		if(m[0]<=0) return na;

		if(this.elCode==0) return model.element[m[0]].getB();
		Vect lc;
		if(this.dim==3){
			lc=localCo(model,m,globalCo);
		}
		else{
			lc=localCoQuad(model,m[0],globalCo);

		}
		return getBAt(model, m[0], lc);

	}

	public double[] getApAnAt(Model model,int i,Vect lc){

		Edge[] elEdges=model.elementEdges(i);

		double[] Ne=NeQuad(lc);


		double[] ApAn=new double[3];

		for(int j=0;j<this.nElEdge;j++)	{	
			ApAn[0]= ApAn[0]+Ne[j]*elEdges[j].Ap;
			ApAn[1]= ApAn[1]+Ne[j]*elEdges[j].A;
		}


		return 	ApAn;
	}


	public Vect getBAt(Model model,int i,Vect lc){
		if(model.fluxLoaded || model.elCode==0) return model.element[i].getB();
		else if(this.elCode==3) model.element[i].getB();

		boolean[] edgeDir=model.element[i].getEdgeReverse();

		Vect B=new Vect(this.dim);
		Mat jac=new Mat(this.dim,this.dim);
		Node[] vertex=model.elementNodes(i);
		Edge[] elEdges=model.elementEdges(i);
		Vect[] rotNe=new Vect[this.nElEdge];
		if(this.dim==3){
			jac=jacobian(vertex,lc);
			rotNe=rotNe(jac,lc,edgeDir);
		}
		else if(this.elCode==1){
			jac=jacobian(vertex,lc);

			rotNe=rotNeQuad(jac,lc);

		}


		for(int j=0;j<this.nElEdge;j++)	{	
			B= B.add(rotNe[j].times(elEdges[j].A));
		}

		if(model.axiSym){
			double rr=0;
			int L=vertex.length;
			for(int j=0;j<L;j++){
				double r1=vertex[j].getCoord(0);
				rr+=1.0/L/r1;
			}
			B=B.times(rr);
			//B.show();
		}

		return B;
	}

	public double[][] Te(Model model,int ie){



		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		int n=this.PWNL[0].length; 


		Node[] vertexNode=model.elementNodes(ie);
		double[][] H=new double[model.nElEdge][model.nElEdge];


		double detJac,ws=1,wsJ=0,term1;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];
		Vect localCo=new Vect(model.dim);


		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];
					localCo.el[2]=this.PWNL[0][r];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;


					rotNe=rotNe(jac,localCo,edgeDir);


					for(int i=0;i<model.nElEdge;i++){

						for(int j=0;j<=i;j++){

							term1=wsJ*rotNe[i].dot(rotNe[j]);
							H[i][j]+=term1;


						}
					}
				}

		lowSym(H);


		return H;


	}




	public double[] nodalMass(Model model,int ie){
		if(this.elCode==0) return nodalMass3ang(model,ie);
		else if(this.elCode==1) return nodalMassQuad(model,ie);
		else if(this.elCode==2) return nodalMassTet(model,ie);
		else if(elCode==3) {return nodalMassPrism(model,ie);}
		Node[] vertexNode=model.elementNodes(ie);
		double[] nm=new double[this.nElVert];

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		double rho=model.element[ie].getRo();

		Vect localCo=new Vect(3);


		double[] N;

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());


					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					N=N(localCo);


					for(int i=0;i<this.nElVert;i++)
						nm[i]+=ws*N[i];
				}


		for(int i=0;i<this.nElVert;i++)
			nm[i]*=rho;


		return nm;
	}


	private double[] nodalMass3ang(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		double[] nm=new double[this.nElVert];
		double rho=model.element[ie].getRo();

		double detJac,ws=1;
		Mat jac;
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		double[] localNe=new double[3];

		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];
			localNe[0]=lc.el[0];
			localNe[1]=lc.el[1];
			localNe[2]=1-lc.el[0]-lc.el[1];


			jac=jacobian3ang(vertexNode,lc);

			detJac=abs(jac.determinant());

			ws=this.PW3ang[p][2]*detJac;
			for(int i=0;i<this.nElVert;i++)

				nm[i]+=localNe[2]*ws;
		}

		for(int i=0;i<this.nElVert;i++)		
			nm[i]*=rho;

		return nm;

	}

	private double[] nodalMassTet(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		double[] nm=new double[this.nElVert];
		double rho=model.element[ie].getRo();

		double detJac,ws=1;
		Mat jac;
		Vect lc=new Vect(dim);		
		int n=this.PWtetra.length;	
		double[] localNe=new double[4];

		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];
			localNe[0]=lc.el[0];
			localNe[1]=lc.el[1];
			localNe[2]=lc.el[2];
			localNe[3]=1-lc.el[0]-lc.el[1]-lc.el[2];

			jac=jacobianTetra(vertexNode,lc);

			detJac=abs(jac.determinant());

			ws=this.PWtetra[p][3]*detJac;
			for(int i=0;i<this.nElVert;i++)			
				nm[i]+=localNe[i]*ws;
		}


		for(int i=0;i<this.nElVert;i++)			
			nm[i]*=rho;



		return nm;

	}
	
	public Vect[] BtSigQuad(Model model,int ie,Vect[] gp_stress){


		Node[] vertexNode=model.elementNodes(ie);
		

		
		Vect[] F=new Vect[this.nElVert];
		
		for(int i=0;i<this.nElVert;i++)
		F[i]=new Vect(model.dim);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		double v=model.element[ie].getPois().el[0];
		double E=model.element[ie].getYng().el[0];
		
		double G=0;
		if(struc2D==0)
			G=E/(1-v*v);
		else
			G=E/((1+v)*(1-2*v));
		
		Mat Bi=new Mat(3,2);;
		//Mat Bj=new Mat(3,2);;
		
		Mat D=new Mat(3,3);;
		
		D.el[0][0]=1;
		D.el[1][1]=1;
		D.el[2][2]=(1.-v)/2;
		D.el[0][1]=v;
		D.el[1][0]=v;
		
		D=D.times(G);


		Vect[] gradN;
		Vect localCo=new Vect(3);
		
		int count=0;
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				
				Vect gp_sig=gp_stress[count++];
				
			
				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());


				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;

				gradN=gradN(jac,localCo);

				


				for(int i=0;i<this.nElVert;i++){
					
					Bi.el[0][0]=gradN[i].el[0];
					Bi.el[1][1]=gradN[i].el[1];
					Bi.el[2][0]=gradN[i].el[1];
					Bi.el[2][1]=gradN[i].el[0];
					
					Mat Bit=Bi.transp();


					Vect BitSig=Bit.mul(gp_sig);
					
		
					F[i]=F[i].add(BitSig.times(ws));

					

					}
			}	


		return F;
	}

	public double[] nodalMassQuad(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		double[] nm=new double[this.nElVert];

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		double rho=model.element[ie].getRo();


		Vect localCo=new Vect(3);


		double[] N;

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());


					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
					else
						ws=detJac;


					N=N(localCo);


					for(int i=0;i<this.nElVert;i++)
						nm[i]+=ws*N[i];
				}


		for(int i=0;i<this.nElVert;i++)
			nm[i]*=rho;


		return nm;
	}


	private double[] nodalMassPrism(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		double[] nm=new double[this.nElVert];
		double rho=model.element[ie].getRo();

		double detJac,ws=1;
		Mat jac;

		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		int nr=this.PW[0].length;	
		for(int r=0;r<nr;r++)
			for(int p=0;p<n;p++)
			{

				lc.el[0]=this.PW3ang[p][0];
				lc.el[1]=this.PW3ang[p][1];
				lc.el[2]=this.PW[0][r];

				localNe=localNPrims(lc);



				jac=jacobianPrism(vertexNode,lc);




				detJac=abs(jac.determinant());



				ws=this.PW3ang[p][2]*this.PW[1][r]*detJac;

				for(int j=0;j<this.nElVert;j++)
					nm[j]+=localNe[j]*ws;

			}

		for(int i=0;i<this.nElVert;i++)
			nm[i]*=rho;

		return nm;

	}




	public Mat SePCQuad(Model model,int ie,Mat Te,Mat Pe,double kw){


		///
		/*double detJac,ws=1;
		Mat jac;
		double gw=model.gw;
		double he=model.getElementT(ie);
		double lam=gw*Math.exp(-.05*he);
		Mat M1=new Mat(this.nElVert,this.nElVert);
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;

				M1=M1.add(Nv.transp().mul(Nv).times(ws));



			}
		 */

		///

		Node[] vertexNode=model.elementNodes(ie);
		Mat Se=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN;
		double[] Nshp;
		Vect lc=new Vect(dim);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,lc);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,lc);
				Nshp=NQuad(lc);

				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(gradN[j])*ws;
						Se.el[i][j]+=BtkB;
						Te.el[i][j]+=Nshp[i]*Nshp[j]*ws;

					}
				}

			}	


		return Se;
	}





	private void lowSym(double[][] H) {
		for(int i=0;i<H.length;i++)
			for(int j=i+1;j<H[0].length;j++)
				H[i][j]= H[j][i];
	}


	public double[][] gaussInteg(int n){
		double[] pg;
		double[] wg;
		if(n==1) {pg=new double[1]; wg=new double[1]; wg[0]=2;}
		else if(n==2) {double a=1.0/sqrt(3.0); pg=new double[2]; pg[0]=-a;pg[1]=a; wg=new double[2]; wg[0]=1;wg[1]=1;}
		else if(n==3) {double a=sqrt(3.0/5.0); pg=new double[3]; pg[0]=-a;pg[1]=0; pg[2]=a;wg=new double[3]; wg[0]=5.0/9;wg[1]=8.0/9;wg[2]=5.0/9;}
		else if(n==4) {double a1=sqrt((3.0-2*sqrt(6.0/5))/7.0); double a2=sqrt((3.0+2*sqrt(6.0/5))/7.0);
		pg=new double[4]; pg[0]=-a2;pg[1]=-a1; pg[2]=a1; pg[3]=a2;wg=new double[4]; double w1=(18-sqrt(30))/36; double w2=(18+sqrt(30))/36;
		wg[0]=w1;wg[1]=w2;wg[2]=w2;wg[3]=w1;}
		else if(n==5) {double a1=sqrt(5.0-2*sqrt(10.0/7))/3.0; double a2=sqrt(5.0+2*sqrt(10.0/7))/3.0;
		pg=new double[5]; pg[0]=-a2;pg[1]=-a1; pg[2]=0;pg[3]=a1; pg[4]=a2;wg=new double[5]; double w1=(322-13*sqrt(70))/900; double w2=(322+13*sqrt(70))/900;
		wg[0]=w1;wg[1]=w2; wg[2]=128.0/225; wg[3]=w2;wg[4]=w1;}
		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		

		double[][] WP=new double[2][wg.length];
		WP[0]=pg;
		WP[1]=wg;

		return WP;
	}

	public double[][] gaussInteg3(int n){

		double[][] PW=new double[n][3];

		if(n==1) {
			PW[0][0]=1.0/3;
			PW[0][1]=1.0/3;
			PW[0][2]=.5;
		}
		else if(n==3) { 
			PW[0][0]=.5; 
			PW[0][1]=.5;
			PW[1][0]=.5; 
			PW[1][1]=0; 
			PW[2][0]=0; 
			PW[2][1]=0.5; 
			double a=1.0/6;
			for(int i=0;i<3;i++)
				PW[i][2]=a;
		}

		else if(n==4) { 
			PW[0][0]=1./3; 
			PW[0][1]=1./3;
			PW[1][0]=.2; 
			PW[1][1]=.6; 
			PW[2][0]=.2; 
			PW[2][1]=0.2; 
			PW[3][0]=.6; 
			PW[3][1]=0.2;

			PW[0][2]=-27./96;

			double a=25./96;
			for(int i=1;i<4;i++)
				PW[i][2]=a;
		}

		else if(n==7) { 
			double c1=1.0/3;
			double c2=1.0/15;
			double c3=1.0/40;

			PW[0][0]=c1; 
			PW[0][1]=c1;
			PW[1][0]=0.5; 
			PW[1][1]=0; 
			PW[2][0]=0.5; 
			PW[2][1]=0.5; 
			PW[3][0]=0; 
			PW[3][1]=0.5;
			PW[4][0]=1; 
			PW[4][1]=0; 
			PW[5][0]=0; 
			PW[5][1]=1; 
			PW[6][0]=0; 
			PW[6][1]=0;

			PW[0][2]=9*c3;


			for(int i=1;i<4;i++)
				PW[i][2]=c2;
			for(int i=4;i<7;i++)
				PW[i][2]=c3;
		}
		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		


		return PW;
	}

	public double[][] gaussIntegTetra(int n){

		double[][] PW=new double[n][4];

		if(n==1) {
			PW[0][0]=1.0/4;
			PW[0][1]=PW[0][0];
			PW[0][2]=PW[0][0];
			PW[0][3]=1.0/6;
		}
		else if(n==4) { 
			double a=(5+3*sqrt(5))/20;
			double b=(5-3*sqrt(5))/20;
			double w=1.0/24;

			PW[0][0]=a; 
			PW[0][1]=b;
			PW[0][2]=b;

			PW[1][0]=b; 
			PW[1][1]=a;
			PW[1][2]=b;

			PW[2][0]=b; 
			PW[2][1]=b;
			PW[2][2]=a;

			PW[3][0]=b; 
			PW[3][1]=b;
			PW[3][2]=b;

			for(int k=0;k<4;k++)
				PW[k][3]=w;

		}

		else if(n==5) { 

			double a=.25;
			double b=1.0/6;
			double w=9.0/120;

			PW[0][0]=a; 
			PW[0][1]=a;
			PW[0][2]=a;

			PW[1][0]=.5; 
			PW[1][1]=b;
			PW[1][2]=b;

			PW[2][0]=b; 
			PW[2][1]=.5;
			PW[2][2]=b;

			PW[3][0]=b; 
			PW[3][1]=b;
			PW[3][2]=.5;

			PW[3][0]=b; 
			PW[3][1]=b;
			PW[3][2]=b;

			for(int k=0;k<5;k++)
				if(k==0)
					PW[k][3]=-4.0/30;
				else
					PW[k][3]=w;

		}

		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		


		return PW;
	}

}
