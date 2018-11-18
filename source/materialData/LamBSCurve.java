package materialData;

import static java.lang.Math.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;

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


public class LamBSCurve {

	public LamBCurve[] lamB;
	public double[] stress;
	public int nstr,nZeroStress;

	public LamBSCurve() throws Exception{

		int code=0;
		
		if(code==0){

			double[] ss1={-4.8832,-3.4008,-2.616,0,1.26,5.66,6.45};
			double[] ss2={-6.1,-1.7,0,3.9};
		
			double[] cc1={1,.7,.5,.3,-.3,-.65,-1};
			double[] cc2={1.2,.9,.3,-.27};
		
			double kk=3;
			double jj=1;
			
			int mm=1;
			
	double[][] ss;
	if(mm==0){
		ss=new double[ss1.length][2];
		for(int i=0;i<ss.length;i++){
			ss[i][0]=ss1[i];
			ss[i][1]=cc1[i];
		}
	}
	else{
		ss=new double[ss2.length][2];
		for(int i=0;i<ss.length;i++){
			ss[i][0]=ss2[i];
			ss[i][1]=cc2[i];
		}
	}

			//double[] cc={1,.7,.5};

	
		BHCurve BH=new BHCurve("H350");

		//double[] ss={-6,-4,-2,0,2,4,6};
		double a1,a2,a3,a4,a5,a6;
		a1=1;
		a2=-.2;
		a3=4.5;;
		a4=0.09;
		a5=1.5;	
		double st=0;
		nstr=ss.length;
		this.lamB=new LamBCurve[nstr];

		stress=new double[nstr];
		
		for(int i=0;i<nstr;i++){
			stress[i]=jj*ss[i][0];
			}


	nZeroStress=-1;
		for(int i=0;i<nstr;i++)
			if(abs(stress[i])<1e-6) {nZeroStress=i; break;}


		this.lamB[nZeroStress]=new LamBCurve("H350");

		double alpha=2,p=1.2, b=.15,c=2,tt=1.2;
		
		int N=181;
		double dB=.01;
		for(int i=0;i<nstr;i++){
			
			//if(i==nposStress) continue;
			double[][] lamB1=new double[N][2];
		
			for(int j=0;j<N;j++){
				double B=j*dB;
				lamB1[j][1]=B;
				
				
				double s;
				s=stress[i];

			
		

			s=s+st;
		
				//lamB1[j][0]=cc[i]*((1.2369)*pow(B,2)+0.914*pow(B,4)-0.3122*pow(B,6));
				lamB1[j][0]=kk*ss[i][1]*((1.0369)*pow(B,2)+0.714*pow(B,4)-0.2122*pow(B,6));
			//lamB1[j][0]=1.4*cc[i]*((1.0369)*pow(B,2)+0.714*pow(B,4)-0.2122*pow(B,6));
				
				//lamB1[j][0]=2*B*B;
		
			}

			lamB[i]=new LamBCurve(lamB1);
			
		}
		
		//makeBHS(stress);
		

		}
		else if(code==1){
			

			double[] ss={-15,1,10,20};

			nstr=ss.length;
			this.lamB=new LamBCurve[nstr];
			double stressMax=2e1;
			
			stress=new double[nstr];
			
		/*	for(int i=0;i<nstr;i++){
				stress[i]=2*(i-.5*(nstr-1))*1.0/(nstr-1)*stressMax;
				}*/
			

		nZeroStress=-1;
			for(int i=0;i<nstr;i++)
				if(abs(stress[i])<1e-6) {nZeroStress=i; break;}

			this.lamB[nZeroStress]=new LamBCurve();

			double alpha=2.0,p=1.2, b=.15,c=2,tt=1.2;
			
			stress=ss;
			nstr=stress.length;
			
			for(int i=0;i<nstr;i++){
				
				double[][] lamB1=new double[this.lamB[nZeroStress].lamB.length][2];
			
				for(int j=0;j<lamB1.length;j++){
					lamB1[j][1]=lamB[nZeroStress].lamB[j][1];
					
					double s;
					s=stress[i];
					double B=lamB[nZeroStress].lamB[j][1];

					//lamB1[j][0]=alpha*(pow(tt-s/stressMax,p)*pow(B,2)-b*(c-s/stressMax)*pow(B,4));
					if(s>18){
					lamB1[j][0]=alpha*(pow(1-s/stressMax,1)*pow(B,2)-.11*pow(B,4));
					}
					else
					if(s>6)
						lamB1[j][0]=.4*alpha*(pow(B,2)-.35*pow(B,4));
						else if(s>=0)
						lamB1[j][0]=1.1*alpha*(pow(B,2)-.2*pow(B,4));
						else if(s<0){
							//if(B<2.4)
							lamB1[j][0]=3*alpha*(pow(B,2)-.18*pow(B,4));
						/*	else
								lamB1[j][0]=3*alpha*(pow(1.4,2)-.5*pow(1.4,4));*/
							//lamB1[j][0]=3*alpha*(1-cos(B*2));
						}
					
					//lamB1[j][0]=lamB[nZeroStress].lamB[j][0]*(1-.0*s)*(1+1./(1+pow(B-.93,2)));
				}
				lamB[i]=new LamBCurve(lamB1);
				
			}
			
			
			//makeBHS(ss);
			
			
		}
		
		if(code==2){
		double[] ss={-50,-30,-10,0,10,30,50};

		BHCurve BH=new BHCurve("belah");
		
	
		double st=0;
		nstr=ss.length;
		this.lamB=new LamBCurve[nstr];
	
		stress=new double[nstr];
		
		for(int i=0;i<nstr;i++){
			stress[i]=ss[i];
			}
		

	nZeroStress=-1;
		for(int i=0;i<nstr;i++)
			if(abs(stress[i])<1e-6) {nZeroStress=i; break;}


		this.lamB[nZeroStress]=new LamBCurve();

		
		int N=201;
		double dB=.01;
		for(int i=0;i<nstr;i++){
			
			//if(i==nposStress) continue;
			double[][] lamB1=new double[N][2];
		
			for(int j=0;j<N;j++){
				double B=j*dB;
				lamB1[j][1]=B;
				
				
				double s;
				s=stress[i];

			
		

			s=s+st;
		
				//lamB1[j][0]=cc[i]*((1.2369)*pow(B,2)+0.914*pow(B,4)-0.3122*pow(B,6));
				lamB1[j][0]=ss[i]*((1.0369)*pow(B,2)+0.714*pow(B,4)-0.2122*pow(B,6));
			//lamB1[j][0]=1.4*cc[i]*((1.0369)*pow(B,2)+0.714*pow(B,4)-0.2122*pow(B,6));
				
				lamB1[j][0]=3*pow((55-ss[i])/50,2)*((.8369)*pow(B,1.6)+0.84*pow(B,4));
		
			}
			lamB[i]=new LamBCurve(lamB1);
			
		}
		
		makeBHS(ss);

		}
		
/*		BHCurve BH=new BHCurve("oita+15");
		
		for(int j=0;j<BH.length;j++)
		{
		util.pr(BH.getH(BH.BH[j][1]));
	
		}*/
		

		
		}
	

	public void makeBHS(double[] ss) throws Exception{
		//BHCurve BH=new BHCurve("belah");
		BHCurve BH=new BHCurve("H350");


		double[] str=new double[ss.length];
		double s;
		double dlds,dldB;
		double E=2.1e11,v=.3,nums,H0,B0;
		double Ex=E/(1+v);
		int L=ss.length;
	
		int N=BH.length;
		N=200;
		double dB=.01;
		double Bst=.01;
		int nbg=1;

		double[][] BHSH=new double[N-nbg][L+1];
		double[][] BHSB=new double[N-nbg][L+1];
		for(int j=0;j<N-nbg;j++){
			
			double Bx=j*dB+Bst;
		
					H0=BH.getH(Bx);
					B0=Bx;
				BHSH[j][0]=B0;
				BHSB[j][0]=H0;
			
		}
			
		for(int j=0;j<L;j++){
			
			
			s=ss[j];
			str[j]=s;
				for(int i=0;i<N-nbg;i++){
					double Bx=i*dB+Bst;
			/*		
					B0=BH.BH[i+nbg][1];
					H0=BH.BH[i+nbg][0];*/
					H0=BH.getH(Bx);
					B0=Bx;
		
				
			dldB=1e6*this.getdLamdB(B0, s);
			dlds=1e-6*this.getdLamds(B0, s);
		
			
			nums=(Ex*dlds-1)*s*dldB/B0;
			//util.pr(nums);
			//if(s<-1 && B0<1) util.pr(nums+ " "+dldB) ;
			
			if(B0<.01 || abs(nums)<.01) {
			BHSB[i][j+1]=B0;
				BHSH[i][j+1]=H0;
			}
			else{
				
				double nu0=H0/(B0+1e-6);
				
				BHSB[i][j+1]=H0/(nu0+nums);
				//util.pr(H0);
			
			BHSH[i][j+1]=H0+B0*nums;

			}
			
		}
		}
		

			
		
			BHSCurve BHSc=new BHSCurve(BHSB,str);
			Curve cv=new Curve(BHSc,600,600);
			cv.show(true);
			
		/*	 file = System.getProperty("user.dir") + "\\BHSoita.xls";
			BHSc.writexls(file);*/
	}

	public  Vect getLam(Vect B,Vect  strs){
		return new Vect(getLam(B.el[0],strs.el[0]),getLam(B.el[1],strs.el[1]),getLam(B.el[2],strs.el[2]));
	}
	
	public  double getLam(double B,double strs){

		if(strs>=this.stress[nstr-1])
			return lamB[nstr-1].getLam(B);
		else if(strs<=this.stress[0])
			return lamB[0].getLam(B);
		else{
		
			int j=getj(strs);
		

		double lam1= this.lamB[j].getLam(B);
		double lam2= this.lamB[j+1].getLam(B);

		return (lam2*(strs-this.stress[j])+lam1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	}
	}
	
	public  double getdLamdB(double B,double strs){
		if(strs>=this.stress[nstr-1])
			return lamB[nstr-1].getdLamdB(B);
		else if(strs<=this.stress[0])
			return lamB[0].getdLamdB(B);
		else{
		
		int j=getj(strs);

		double lam1= this.lamB[j].getdLamdB(B);
		double lam2= this.lamB[j+1].getdLamdB(B);
		double dldB=(lam2*(strs-this.stress[j])+lam1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
		return dldB;
		
		
				
	}
	}
	
	public  double getd2LamdB2(double B,double strs){
		if(strs>=this.stress[nstr-1])
			return lamB[nstr-1].getd2LamdB2(B);
		else if(strs<=this.stress[0])
			return lamB[0].getd2LamdB2(B);
		else{
		
		int j=getj(strs);

		double lam1= this.lamB[j].getd2LamdB2(B);
		double lam2= this.lamB[j+1].getd2LamdB2(B);
		return (lam2*(strs-this.stress[j])+lam1*(this.stress[j+1]-strs))/(this.stress[j+1]-this.stress[j]);
				
	}
	}

	public  double getdLamds(double B,double strs){
		if(nstr==1) return 0;
		if(strs>=this.stress[nstr-1])
			return (lamB[nstr-1].getLam(B)-lamB[nstr-2].getLam(B))/(this.stress[nstr-1]-this.stress[nstr-2]);
		else if(strs<=this.stress[0])
			return (lamB[1].getLam(B)-lamB[0].getLam(B))/(this.stress[1]-this.stress[0]);
		else{
		
		int j=getj(strs);

		double lam0,lam1,lam2,s0,s1,s2,m1,m2;
		 lam0= this.lamB[j].getLam(B);
		 lam1= this.lamB[j+1].getLam(B);
		 s0=this.stress[j];
		s1=this.stress[j+1];
		m1=(lam1-lam0)/(s1-s0);
		if(j==nstr-2){
			return m1;
		}

		 lam2= this.lamB[j+2].getLam(B);		
		 s2=this.stress[j+2];
	
		 m2=(lam2-lam1)/(s2-s1);
		return (m2*(strs-s0)+m1*(s2-strs))/(s2-s0);
				
	}
	}


	public  double[][] getLamBSArray(){
		
	double[][] ar=new double[this.lamB[0].length][this.lamB.length+1];
	
	for(int i=0;i<ar.length;i++){
		ar[i][0]=this.lamB[0].lamB[i][1];
	
	}

	
	
	for(int i=1;i<ar[0].length;i++)
		for(int j=0;j<ar.length;j++){
		ar[j][i]=this.lamB[i-1].lamB[j][0];
		
	}

	
	return ar;
	
}
	
	public int getj(double strs){

		int j=0;
		while(this.stress[j+1]<strs){  j++;}
		
		
		return j;
	}
	
	public void writexls(String file) throws RowsExceededException, WriteException, IOException{
		
		WritableWorkbook workbook = Workbook.createWorkbook(new File(file)); 
		WritableSheet sheet = workbook.createSheet("First Sheet",0);

		int L=this.lamB[0].length;

		int N=this.lamB.length;
		Number[][] number=new Number[L][N+1];
		int c0=6,r0=2;
	Label[] stress=new Label[L];
		
		DecimalFormat df=new DecimalFormat("0.00");
		for(int p=0;p<N;p++){
			stress[p] = new Label(r0+p+1, 5, df.format(this.stress[p])+" MPa");
			sheet.addCell(stress[p]);
		}
		for(int i=0;i<L;i++)
			number[i][0]=new Number(r0, c0+i, this.lamB[0].lamB[i][1]);
		for(int i=0;i<L;i++)
			for(int p=0;p<N;p++){
				number[i][p+1] = new Number(r0+p+1, c0+i,this.lamB[p].lamB[i][0]);
		
		}
		
		Label label = new Label(2,2,"Lambda-stress data");
		sheet.addCell(label);

	for(int i=0;i<L;i++)
		for(int p=0;p<=N;p++){
			sheet.addCell(number[i][p]);
		}
			

		workbook.write();
		workbook.close();
		
		}
	
	public void showCurve(){
		Curve cv=new Curve(this,800,600);
		cv.show(true);

	}


	public static void main(String[] args) throws Exception {
		
		int L=50;
/*		double[][] xy=new double[L][2];
		for(int i=0;i<L;i++)
			xy[i][0]=i*.2*1e6;	*/
		String file = System.getProperty("user.dir") + "\\BH\\"+"lambs.txt";
		LamBSCurve lamBS=new LamBSCurve();
	//	System.out.println(lamBS.getLam(1,-2));
		double s=-5;
		double B=1.13;
	
		System.out.println(lamBS.getLam(B,s));
		System.out.println(lamBS.getdLamds(B,s));
		
		double C;
		if(s>0) C=4;
		else
			C=10;
double mm;

mm=-C/PI*(3e-7)*lamBS.getLam(B, 0)/(1+pow(s*3e-7,2));

		System.out.println(mm);
	//	System.out.println(lamBS.getdLamdB(1.1,1e7));
		
	/*	 file = System.getProperty("user.dir") + "\\LambS.xls";
		lamBS.writexls(file);*/
		
Curve cv=new Curve(lamBS,800,600);
cv.show(true);
	}




}
