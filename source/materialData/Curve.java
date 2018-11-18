package materialData;

import java.awt.Color;
import java.awt.Font;
import java.text.DecimalFormat;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.math.plot.Plot2DPanel;

import materialData.BHCurve;
import materialData.BHSCurve;
import materialData.LamBCurve;
import materialData.LamBSCurve;
import math.Vect;
import math.util;
public class Curve extends JPanel{

	JFrame frame;

	public Curve(BHCurve BH,int width, int height){
		DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		plot.setFont(new Font("Arialw Roman", 2, 19));

		double[] H=new double[BH.length];
		double[] B=new double[H.length];
		for(int i=0;i<BH.length;i++)
			H[i]=1e-3*BH.BH[i][0];

	
			for(int i=0;i<B.length;i++)
			B[i]=BH.BH[i][1];
			plot.addLinePlot("B-H ", H, B);
		
					
		plot.setAxisLabel(0,"H (kA/m)");
		plot.setAxisLabel(1,"B (T)");
		setFrame("B-H curves ",width,height);
		frame.add(plot);

	}

	
	public Curve(LamBCurve lamB,int width, int height){

		
		Plot2DPanel plot = new Plot2DPanel();
		double[] B=new double[lamB.length];
		double[] lam=new double[B.length];
		for(int i=0;i<lamB.lamB.length;i++)
			B[i]=lamB.lamB[i][1];

	
			for(int i=0;i<B.length;i++)
			lam[i]=lamB.lamB[i][0];
	
			DecimalFormat df=new DecimalFormat("00.0");	
		
			plot.addLinePlot("lamB ", B, lam);
		
		plot.setAxisLabel(0,"B");
		plot.setAxisLabel(1,"lamda");
		setFrame("Magnetostriction Curve",width,height);
		frame.add(plot);

	
	
	}
	
	public Curve(CurrentWaveForm Ix,int width, int height){

		
		Plot2DPanel plot = new Plot2DPanel();
		double[] t=new double[Ix.length];
		double[] ix=new double[t.length];
		for(int i=0;i<Ix.TI.length;i++)
			t[i]=Ix.TI[i][0];

	
			for(int i=0;i<t.length;i++)
			ix[i]=Ix.TI[i][1];
		
			plot.addLinePlot("Ix ", t, ix);
		
		plot.setAxisLabel(0,"t");
		plot.setAxisLabel(1,"I");
		setFrame("Current wave form",width,height);
		frame.add(plot);

	
	
	}
	
	public Curve(BHSCurve BHS,int width, int height){
		DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		double[] H=new double[BHS.BH[0].BH.length];
		double[] B=new double[H.length];
		for(int i=0;i<BHS.BH[0].BH.length;i++)
			H[i]=BHS.BH[0].BH[i][0];

		for(int j=0;j<BHS.nstr;j++){
			for(int i=0;i<B.length;i++)
			B[i]=BHS.BH[j].BH[i][1];
			plot.addLinePlot(df.format(BHS.stress[j])+" MPa ", H, B);
		
				}
		
		plot.setAxisLabel(0,"H");
		plot.setAxisLabel(1,"B");
		plot.addLegend("EAST");
		setFrame("B-H curves for various stress levels",width,height);
		frame.add(plot);

	}
	

	public Curve(LamBSCurve lamBS,int width, int height){
	
		Plot2DPanel plot = new Plot2DPanel();
		double[] B=new double[lamBS.lamB[0].lamB.length];
		double[] lam=new double[B.length];
		for(int i=0;i<lamBS.lamB[0].lamB.length;i++)
			B[i]=lamBS.lamB[0].lamB[i][1];

		for(int j=0;j<lamBS.nstr;j++){
			for(int i=0;i<B.length;i++)
			lam[i]=lamBS.lamB[j].lamB[i][0];
	
			DecimalFormat df=new DecimalFormat("00.0");	
			plot.addLinePlot(df.format(lamBS.stress[j])+" MPa ", B, lam);
		
		
		}
		
		plot.setAxisLabel(0,"B");
		plot.setAxisLabel(1,"lamda");
		plot.addLegend("EAST");
		setFrame("lambda curve for various stress levels",width,height);
		frame.add(plot);

	}
	

	private void setFrame(String name,int width, int height){
		frame=new JFrame(name);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(width,height);
	}

	public void show(boolean b){
		frame.setVisible(b);
	}

}