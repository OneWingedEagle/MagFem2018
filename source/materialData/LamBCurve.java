package materialData;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.Scanner;

import jxl.Workbook;
import jxl.write.Label;
import jxl.write.Number;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import jxl.write.biff.RowsExceededException;
import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;



public class LamBCurve {
	
		public double[][] lamB;
		public double[][] gradLamB;
		public double[][] grad2LamB;
		public int length;
		public double scale=1e-6;

		public LamBCurve(){
			
			this.lamB=new double[200][2];
			this.length=this.lamB.length;
			for(int i=0;i<this.lamB.length;i++){
				this.lamB[i][1]=.01*i;
				this.lamB[i][0]=1*(Math.pow(this.lamB[i][1],2)-1*.2*Math.pow(this.lamB[i][1],4));
				
			}
			 setGradLamB();
			 setGrad2LamB();
		}
		
		public LamBCurve(String file){

			this.lamB=new double[20][2];
			this.length=this.lamB.length;
			for(int i=0;i<this.lamB.length;i++){

				this.lamB[i][1]=.1*i;
				//this.lamB[i][0]=10*(Math.pow(this.lamB[i][1],2)-1*.41*Math.pow(this.lamB[i][1],4)); // Hameyer
				this.lamB[i][0]=.8*(Math.pow(this.lamB[i][1],2)-1*.2*Math.pow(this.lamB[i][1],4))+.0*util.triangWave(1.,5*this.lamB[i][1]); 
				
				
			}
			
			
			 setGradLamB();
			 setGrad2LamB();
			
		}
		

	
	public LamBCurve(double[][] lamB){
		this.length=lamB.length;
		this.lamB=new double[this.length][2];
		for(int i=0;i<this.length;i++)
			for(int j=0;j<2;j++)
			this.lamB[i][j]=lamB[i][j];
		 setGradLamB();
		setGrad2LamB();
	}
			
	public void setGradLamB(){

		 this.gradLamB=new double[this.lamB.length][2];
		 for(int k=0;k<this.gradLamB.length-1;k++){

			 this.gradLamB[k][0]=(this.lamB[k+1][1]-this.lamB[k][1]);
			
			 this.gradLamB[k][1]=(this.lamB[k+1][0]-this.lamB[k][0])/this.gradLamB[k][0];
			
		 }
		 this.gradLamB[this.gradLamB.length-1][0]=1000;
		 this.gradLamB[this.gradLamB.length-1][1]=this.gradLamB[this.gradLamB.length-2][1];


		
	}
	
	public void setGrad2LamB(){

		 
		 this.grad2LamB=new double[this.lamB.length][2];
		 for(int k=0;k<this.grad2LamB.length-1;k++){
			 this.grad2LamB[k][0]= this.gradLamB[k][0];
			 this.grad2LamB[k][1]=(this.gradLamB[k+1][1]-this.gradLamB[k][1])/this.grad2LamB[k][0];
			
		 }
		
		 this.grad2LamB[this.grad2LamB.length-1][0]=10;
		 this.grad2LamB[this.grad2LamB.length-1][1]=this.grad2LamB[this.grad2LamB.length-2][1];
		 this.grad2LamB[0][1]=this.grad2LamB[1][1];
	
	}
	
public void writexls(String file) throws RowsExceededException, WriteException, IOException{
		
		WritableWorkbook workbook = Workbook.createWorkbook(new File(file)); 
		WritableSheet sheet = workbook.createSheet("Lambda curve",0);


		int L=this.lamB.length;
		Number[][] number=new Number[L][2];
		int c0=6,r0=2;
		for(int i=0;i<L;i++){
			number[i][0]=new Number(r0, c0+i, this.lamB[i][0]);
			number[i][1]=new Number(r0+1, c0+i, this.lamB[i][1]);
		}
		
		Label label = new Label(2,2,"Lambda data");
		sheet.addCell(label);

	for(int i=0;i<L;i++)
	{
			sheet.addCell(number[i][1]);
			sheet.addCell(number[i][0]);
		}
			

		workbook.write();
		workbook.close();
		
		}

		
	
	public static void main(String[] args) throws Exception{
		String file = System.getProperty("user.dir") + "\\LamB\\50H400.txt";
		
		LamBCurve lamB=new LamBCurve(file);
		
		util.show(lamB.lamB);
	
		Curve cv=new Curve(lamB,800,600);

cv.show(true);
		
		 file = System.getProperty("user.dir") + "\\Lamb1.xls";
			lamB.writexls(file);

	
	}
	
	public  double getLam(double B){
	

		if(B>=this.lamB[this.length-1][1])
			return this.scale*(this.lamB[this.length-1][0]+(B-this.lamB[this.length-1][1])*this.gradLamB[this.length-1][1]);
		int j=0;
		while(this.lamB[j+1][1]<B){j++;}
	double lam=(this.lamB[j][0]*(this.lamB[j+1][1]-B)+this.lamB[j+1][0]*(B-this.lamB[j][1]))/this.gradLamB[j][0];
	
		return this.scale*lam;
				
	}
	

	
	public  double getdLamdB(double B){
		
		
		if(B>=this.lamB[this.length-1][1])
			return this.scale*this.gradLamB[this.length-1][1];
		int j=0;
		while(this.lamB[j+1][1]<B){j++;}

		if(j>=this.length-1)
			return	this.scale*this.gradLamB[this.length-2][1];
		double  interp=((B-this.lamB[j][1])*this.gradLamB[j+1][1]+(this.lamB[j+1][1]-B)*this.gradLamB[j][1])/this.gradLamB[j][0];
	
		return this.scale*interp;
				
	}
	
	public  double getd2LamdB2(double B){
		
		
		if(B>=this.lamB[this.length-1][1])
			return this.scale*this.grad2LamB[this.length-1][1];
		int j=0;
		while(this.lamB[j+1][1]<B){j++;}

		if(j>=this.length-1)
			return	this.scale*this.grad2LamB[this.length-2][1];
		double  interp=((B-this.lamB[j][1])*this.grad2LamB[j+1][1]+(this.lamB[j+1][1]-B)*this.grad2LamB[j][1])/this.grad2LamB[j][0];
		

		return this.scale*interp;
				
	}
	
	public  double[][] getLamBArray(){
		
		double[][] lamB1=new double[this.lamB.length][this.lamB[0].length];
		for(int i=0;i<this.length;i++)
		for(int j=0;j<2;j++)
			lamB1[i][j]=this.lamB[i][1-j];

		return lamB1;
	
				
	}
	
	public void showCurve(){
		Curve cv=new Curve(this,800,600);
		cv.show(true);

	}
	
	
}
