package materialData;

import io.Loader;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;


import jxl.Workbook;
import jxl.write.Label;
import jxl.write.Number;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import jxl.write.biff.RowsExceededException;
import fem.Model;
import math.Mat;
import math.Vect;
import math.util;
import static java.lang.Math.abs;

public class BHSCurve {

	public BHCurve[] BH;
	public double[] stress;
	public int nstr,nZeroStress;

	public BHSCurve(){}
	
	public BHSCurve(int L){
		
		nstr=L;
		BH=new BHCurve[L];
		stress=new double[L];
		
	}
	
	public BHSCurve(double[][] HBS,double[] str) throws Exception{
		
		nstr=str.length;
		stress=str;
		BH=new BHCurve[nstr];
		int N=HBS.length;
		for(int i=0;i<nstr;i++){
			double[][] BHx=new double[N][2];
			for(int j=0;j<N;j++){
				BHx[j][0]=HBS[j][0];
				BHx[j][1]=HBS[j][i+1];
			}
				
			BH[i]=new BHCurve(BHx);
		}
		
	}
	

	
	public BHSCurve(String mateName, boolean vms) throws Exception{

	Loader loader=new Loader();
	loader.loadBHS(this,mateName,vms);
	

	}
	
	
	
	
	public BHSCurve(Model model, int ir, double[] stress,int nLamBS) throws Exception{
	util.hshow(stress);
		BHSCurve BHS1=new BHSCurve(model.region[ir].getMaterial(),true);
		int nz=BHS1.nZeroStress;
		BHCurve BH1=BHS1.BH[nz];

		double E=model.region[ir].getYng().el[0];
		double v=model.region[ir].getPois().el[0];
		LamBSCurve lamBS=model.lamBS[nLamBS];
		double s;
		double dlds,dldB,H0,B0,nums;
		double G=E/(1+v);
		int L=stress.length;
		this.stress=stress;
		this.nstr=L;
		BH=new BHCurve[nstr];
			
		BHSCurve BHS2=new BHSCurve(L);
		
		int K=BH1.length;
		
		int p=2;
		int N=p*(K-1)+1;
	
		double[][] BHx=new double[N][2];
		double[] B2=new double[N];
		
		for(int j=0;j<K-1;j++){		
			double dB=(BH1.BH[j+1][1]-BH1.BH[j][1])/p;		
			for(int k=0;k<p;k++){
				 B2[j*p+k]=BH1.BH[j][1]+k*dB;
			}
				
		}
		 B2[N-1]=BH1.BH[K-1][1];

		for(int j=0;j<L;j++){
			
			s=stress[j];
				for(int i=0;i<N;i++){
					B0=B2[i];
					H0=BH1.getH(B0);					
					BHx[i][1]=B0;

			dldB=lamBS.getdLamdB(B0, s);			

			dlds=1e-6*lamBS.getdLamds(B0, s);	
			nums=(G*dlds-1)*1e6*s*dldB/(B0+1e-6);
	
			if(B0<.001) {
				BHx[i][0]=H0;
				}
				else{					
					BHx[i][0]=H0+nums*B0;
		}	
			
		}
		
				BHS2.BH[j]=new BHCurve(BHx);
	}
		
		
	
		
		N=p*(K-1)+1;
		BHx=new double[N][2];
		
		double[] H2=new double[N];
		
			for(int j=0;j<K-1;j++){
				
				double dB=(BH1.BH[j+1][1]-BH1.BH[j][1])/p;
				
				for(int k=0;k<p;k++){
					double B3=BH1.BH[j][1]+k*dB;
					H2[j*p+k]=BH1.getH(B3);
				
				}
			}
			H2[N-1]=BH1.BH[K-1][0];

			for(int i=0;i<nstr;i++){
		
			for(int j=0;j<N;j++){
				BHx[j][0]=H2[j];			
				BHx[j][1]=BHS2.BH[i].getB(H2[j]);
			}
		
			BH[i]=new BHCurve(BHx);
		}

		
		for(int i=0;i<nstr;i++)
			if(abs(stress[i])<1e-5) {nZeroStress=i; break;}
	
	/*	String bhs = System.getProperty("user.dir") + "\\belahBHS.xls";
		this.writexls(bhs);*/
		//BHSCurve BHS2=new BHSCurve("oitaBHS",false);
	
		
		/*Curve cv=new Curve(this,600,600);
		cv.show(true);*/
	}

	
	
	public  Vect getH(Vect B,Vect  strs){
		return new Vect(getH(B.el[0],strs.el[0]),getH(B.el[1],strs.el[1]),getH(B.el[2],strs.el[2]));
	}
	
	public  double getH(double B,double strs){

		if(strs>=this.stress[nstr-1])
			return BH[nstr-1].getH(B);
		else if(strs<=this.stress[0])
			return BH[0].getH(B);
		else{
		
			int j=getj(strs);
		double H1= BH[j].getH(B);
		double H2= BH[j+1].getH(B);

		return (H2*(strs-this.stress[j])+H1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	}
	}
	
		public  double getPem(double H1,double H2,double strs){
			if(this.nstr==1)
			return BH[0].getPem(H1,H2);
			
			if(strs>=this.stress[nstr-1]){
				return BH[nstr-1].getPem(H1, H2);
			}
		
			int j=getj(strs);
	
			double pem1=BH[j].getPem(H1,H2);
			double pem2=BH[j+1].getPem(H1,H2);
		
			return (pem2*(strs-this.stress[j])+pem1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
			

	}
	
		public  double getPem(Vect H,Vect strs){
			
			return  getPem(new Vect(3), H,strs);
		}
		
		public double getPem(Vect H1, Vect H2, Vect strs){
			double pem=0;
			
			pem=getPem(H1.el[0],H2.el[0],strs.el[0])+getPem(H1.el[1],H2.el[1],strs.el[1])+getPem(H1.el[2],H2.el[2],strs.el[2]);
			
			return pem;
			
		}

		public  double getPemB(double B1,double B2,double strs){
			if(this.nstr==1)
			return BH[0].getPemB(B1,B2);
			
			if(strs>=this.stress[nstr-1]){
				return BH[nstr-1].getPemB(B1, B2);
			}
		
			int j=getj(strs);
	
			double pem1=BH[j].getPemB(B1,B2);
			double pem2=BH[j+1].getPemB(B1,B2);
		
			return (pem2*(strs-this.stress[j])+pem1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
			

	}
	
	public  double getB(double H,double strs){
		if(strs>=this.stress[nstr-1])
			return BH[nstr-1].getB(H);
		else if(strs<=this.stress[0])
			return BH[0].getB(H);

		
		int j=getj(strs);

		double B1= BH[j].getB(H);
		double B2= BH[j+1].getB(H);
		return (B2*(strs-this.stress[j])+B1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	
	}
	
public  double getdHdB(double B,double strs){

		
		if(strs>=this.stress[nstr-1])
			return BH[nstr-1].getdHdB(B);
		else if(strs<=this.stress[0])
			return BH[0].getdHdB(B);

		int j=getj(strs);
		

		double dHdB1= BH[j].getdHdB(B);
		double dHdB2= BH[j+1].getdHdB(B);
		return (dHdB2*(strs-this.stress[j])+dHdB1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	}
	
	public  double getNu(double B,double strs){

		if(strs>=this.stress[nstr-1])
			return BH[nstr-1].getNu(B);
		else if(strs<=this.stress[0])
			return BH[0].getNu(B);
		
		int j=getj(strs);
		
		double nu1= BH[j].getNu(B);
		double nu2= BH[j+1].getNu(B);
		double nu=(nu2*(strs-this.stress[j])+nu1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
		
		return nu;
	}
	
	public  Vect getNu(Vect B,Vect stress){
		Vect nu=new Vect(this.getNu(abs(B.el[0]),stress.el[0]),this.getNu(abs(B.el[1]),stress.el[1]),this.getNu(abs(B.el[2]),stress.el[2]));
		return nu;
	}
		
	public  double getNuVar(double B,double strs){
		if(strs>=this.stress[nstr-1])
			return BH[nstr-1].getNuVar(B);
		else if(strs<=this.stress[0])
			return BH[0].getNuVar(B);

		int j=getj(strs);
		
		double nuVar1= BH[j].getNuVar(B);
		double nuVar2= BH[j+1].getNuVar(B);
	
		return (nuVar2*(strs-this.stress[j])+nuVar1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	}
	
	public  double getdNuds(double B,double strs){
		if(this.nstr==1)return 0;
		if(strs>=this.stress[nstr-1])
			return (BH[nstr-1].getNu(B)-BH[nstr-2].getNu(B))/(this.stress[nstr-1]-this.stress[nstr-2]);
		else if(strs<=this.stress[0])
			return (BH[1].getNu(B)-BH[0].getNu(B))/(this.stress[1]-this.stress[0]);
		else{
		
		int j=getj(strs);

		double nu0,nu1,nu2,s0,s1,s2,m1,m2;
		 nu0= this.BH[j].getNu(B);
		 nu1= this.BH[j+1].getNu(B);
		 s0=this.stress[j];
		s1=this.stress[j+1];
		m1=(nu1-nu0)/(s1-s0);
		if(j==nstr-2){
			return m1;
		}

		 nu2= this.BH[j+2].getNu(B);		
		 s2=this.stress[j+2];
	
		 m2=(nu2-nu1)/(s2-s1);
		return (m2*(strs-s0)+m1*(s2-strs))/(s2-s0);
					
	}
	}
	
	public  double getdPemBds(double B1,double B2,double strs){
		//if(strs<1e1) return 0;
		if(this.nstr==1)return 0;
		if(strs>=this.stress[nstr-1])
			return (BH[nstr-1].getPemB(B1,B2)-BH[nstr-2].getPemB(B1,B2))/(this.stress[nstr-1]-this.stress[nstr-2]);
		else if(strs<=this.stress[0])
			return (BH[1].getPemB(B1,B2)-BH[0].getPemB(B1,B2))/(this.stress[1]-this.stress[0]);
		else{
		
		int j=getj(strs);

		double pemB0,pemB1,pemB2,s0,s1,s2,m1,m2;
		 pemB0= this.BH[j].getPemB(B1,B2);
		 pemB1= this.BH[j+1].getPemB(B1,B2);
		 s0=this.stress[j];
		s1=this.stress[j+1];
		m1=(pemB1-pemB0)/(s1-s0);
		
		if(j==nstr-2){
			return m1;
		}

		 pemB2= this.BH[j+2].getPemB(B1,B2);
		 s2=this.stress[j+2];
		 m2=(pemB2-pemB1)/(s2-s1);
		 double m=(m2*(strs-s0)+m1*(s2-strs))/(s2-s0);
		return m;
					
	}
	}
	public void writexls(String file) throws RowsExceededException, WriteException, IOException{
			
		WritableWorkbook workbook = Workbook.createWorkbook(new File(file)); 
		WritableSheet sheet = workbook.createSheet("First Sheet",0);

		int L=this.BH[0].length;

		int N=this.BH.length;
		
	
		
		Number[][] number=new Number[L][N+1];
		
		int c0=6,r0=2;
		Label[] stress=new Label[L];
		
		DecimalFormat df=new DecimalFormat("0.00");
		for(int p=0;p<N;p++){
			stress[p] = new Label(r0+p+1, 5, df.format(this.stress[p])+" MPa");
			sheet.addCell(stress[p]);
		}
		
		for(int i=0;i<L;i++)
			number[i][0]=new Number(r0, c0+i, this.BH[0].BH[i][0]);
		for(int i=0;i<L;i++)
			for(int p=0;p<N;p++){
				number[i][p+1] = new Number(r0+p+1, c0+i,this.BH[p].BH[i][1]);
		
		}
		
		Label label = new Label(2,2,"BH-stress data");
		sheet.addCell(label);

	for(int i=0;i<L;i++)
		for(int p=0;p<=N;p++){
			sheet.addCell(number[i][p]);
		}
			

		workbook.write();
		workbook.close();
		
		}
	
	
	public int getj(double strs){
		
		int j=0;
		if(this.stress.length==1) return j;
		while(this.stress[j+1]<strs){j++;}
		return j;
	}
	
	public void showCurve(){
		Curve cv=new Curve(this,800,600);
		cv.show(true);

	}

	public static void main(String[] args) throws Exception {
	
		//BHSCurve BHS=new BHSCurve("NipponSteel");
		//HSCurve BHS=new BHSCurve("steelLinear");
		BHSCurve BHS=new BHSCurve("H350",false);
//for(int i=0;i<8;i++)
util.pr(BHS.getdPemBds(0,1,0));	
Curve cv=new Curve(BHS,800,600);
cv.show(true);
/*
String bhs = System.getProperty("user.dir") + "\\ttt.xls";
BHS.writexls(bhs);*/


	}




}
